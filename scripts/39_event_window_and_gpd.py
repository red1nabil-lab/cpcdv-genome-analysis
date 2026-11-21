#!/usr/bin/env python
"""Step 39: (a) does the event grouping depend on window size?  (b) is the GPD step sound?

Reviewer round 3.

Point 3: genomes are collapsed into "distinct events" by sharing a 50 nt breakpoint
window and a parental strain pair.  That is an arbitrary window, so the enrichment is
recomputed at 25, 50 and 100 nt and the per-region event counts are printed so the
enrichment ratio can be checked by hand.

Point 6(ii): the claim that the generalized Pareto step is conservative compared
empirical significance under one multiple-testing family with extrapolated significance
under another, which proves nothing.  Two proper checks instead:
  * like-for-like - empirical and GPD p-values under the SAME correction
  * a false-positive check - run the whole pipeline on data with the recombination
    signal destroyed (alignment columns randomly permuted, which preserves the site
    frequency spectrum and the tree but removes spatial structure) and count calls
"""
import os, sys, json, csv, time
from collections import Counter, defaultdict
import numpy as np
from scipy import stats

sys.stdout.reconfigure(encoding="utf-8")
import os as _os, sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # noqa: E402
from cpcdv_paths import DATA_DIR, FIG_DIR, MUSCLE_BIN  # noqa: E402
AN = DATA_DIR
rng = np.random.default_rng(20261124)
GAP = ord("-")
MIN_SIDE, N_PERM, GPD_EXC = 10, 10000, 250

rec = [r for r in csv.DictReader(open(os.path.join(AN, "29_recombination_maxchi.csv"),
                                      encoding="utf-8-sig"))
       if r["significant"] == "True"]

REGIONS = [("LIR (5')", 0, 181), ("V1 / MP", 181, 507), ("spacer", 507, 519),
           ("V2 / CP", 519, 1264), ("SIR", 1264, 1437), ("Rep exon 1", 1437, 1859),
           ("intron", 1859, 1945), ("Rep exon 2", 1945, 2632), ("LIR (3')", 2632, 2829)]
HOT = ("SIR", "Rep exon 1")
REF = ("V1 / MP", "V2 / CP")          # the comparison is against the two structural
                                      # genes, NOT "coding regions" -- Rep exon 1 codes


def region_of(p):
    for nm, s, e in REGIONS:
        if s <= p < e:
            return nm
    return None


print("=== (a) sensitivity of the event grouping to window size ===")
rows = []
for W in (25, 50, 100):
    ev = defaultdict(list)
    for r in rec:
        bp = int(r["breakpoint_aln_col"])
        ev[(bp // W, tuple(sorted((r["parent1_strain"], r["parent2_strain"]))))]\
            .append(r["accession"])
    cnt = Counter()
    for (wnd, _), accs in ev.items():
        nm = region_of(wnd * W + W // 2)
        if nm:
            cnt[nm] += 1
    hw = sum(e - s for nm, s, e in REGIONS if nm in HOT)
    rw = sum(e - s for nm, s, e in REGIONS if nm in REF)
    h = sum(cnt.get(n, 0) for n in HOT)
    rf = sum(cnt.get(n, 0) for n in REF)
    ratio = (h / hw) / (rf / rw) if rf else float("inf")
    print(f"\n  window {W:3d} nt -> {len(ev)} distinct events")
    print(f"    {'region':12s} {'width':>6s} {'events':>7s} {'per kb':>8s}")
    for nm, s, e in REGIONS:
        c = cnt.get(nm, 0)
        print(f"    {nm:12s} {e-s:6d} {c:7d} {1000*c/(e-s):8.2f}")
    print(f"    SIR+Rep exon 1: {h} events / {hw} nt = {1000*h/hw:.2f} per kb")
    print(f"    MP+CP         : {rf} events / {rw} nt = {1000*rf/rw:.2f} per kb")
    print(f"    enrichment    : {ratio:.1f}-fold")
    rows.append(dict(window_nt=W, n_events=len(ev), hot_events=h, hot_width=hw,
                     ref_events=rf, ref_width=rw, enrichment=round(ratio, 2),
                     **{f"events_{nm}": cnt.get(nm, 0) for nm, _, _ in REGIONS}))

out = os.path.join(AN, "39_event_window_sensitivity.csv")
with open(out, "w", newline="", encoding="utf-8-sig") as fh:
    w_ = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w_.writeheader(); w_.writerows(rows)
print(f"\n  wrote {out}")

# ============================================================ (b) GPD checks
print("\n=== (b) checks on the generalized Pareto step ===")


def read_fasta(p):
    d, k, buf = {}, None, []
    for line in open(p, encoding="utf8"):
        line = line.strip()
        if line.startswith(">"):
            if k:
                d[k] = "".join(buf)
            k, buf = line[1:].split()[0], []
        elif line:
            buf.append(line)
    if k:
        d[k] = "".join(buf)
    return d


def maxchi(x):
    n = x.size
    if n < 2 * MIN_SIDE:
        return 0.0, -1
    c = np.cumsum(x); tot = c[-1]
    k = np.arange(MIN_SIDE, n - MIN_SIDE + 1)
    a = c[k - 1]; b = tot - a
    left, right = k, n - k
    num = n * (a * (right - b) - (left - a) * b) ** 2
    den = left * right * tot * (n - tot)
    with np.errstate(divide="ignore", invalid="ignore"):
        chi = np.where(den > 0, num / den, 0.0)
    j = int(np.nanargmax(chi))
    return float(chi[j]), int(k[j])


def null_dist(x, nperm):
    n = x.size
    order = rng.random((nperm, n)).argsort(axis=1)
    P = x[order]
    C = np.cumsum(P, axis=1); tot = C[:, -1]
    k = np.arange(MIN_SIDE, n - MIN_SIDE + 1)
    a = C[:, k - 1]; b = tot[:, None] - a
    left, right = k[None, :], (n - k)[None, :]
    num = n * (a * (right - b) - (left - a) * b) ** 2
    den = left * right * tot[:, None] * (n - tot)[:, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        chi = np.where(den > 0, num / den, 0.0)
    return np.nanmax(chi, axis=1)


# --- check 1: does the GPD fit reproduce quantiles it was not fitted to?
print("\n  check 1: GPD tail fit validated against held-out permutation quantiles")
aln = read_fasta(os.path.join(AN, "05_full_aln.fasta"))
keep = [a for a in json.load(open(os.path.join(AN, "10_analysis_set.json")))["keep"]
        if a in aln]
A = np.frombuffer(b"".join(aln[a].upper().encode() for a in keep),
                  dtype=np.uint8).reshape(len(keep), -1)
ixm = {a: i for i, a in enumerate(keep)}
trial = [r for r in rec[:6]]
errs = []
for r in trial:
    q, p1, p2 = r["accession"], r["parent1"], r["parent2"]
    iq, a1, a2 = A[ixm[q]], A[ixm[p1]], A[ixm[p2]]
    inf = (a1 != a2) & (a1 != GAP) & (a2 != GAP) & (iq != GAP) & ((iq == a1) | (iq == a2))
    x = (iq[inf] == a1[inf]).astype(np.int64)
    nd = null_dist(x, 20000)
    fit_part, test_part = nd[:10000], nd[10000:]
    tail = np.sort(fit_part)[-GPD_EXC:]
    u = tail[0]; exc = tail - u
    shp, loc, scl = stats.genpareto.fit(exc[exc > 0], floc=0.0)
    for qq in (0.999, 0.9995, 0.9999):
        emp = np.quantile(test_part, qq)
        pred_p = (GPD_EXC / 10000) * float(stats.genpareto.sf(emp - u, shp, loc=loc,
                                                              scale=scl))
        emp_p = float((test_part >= emp).mean())
        if emp_p > 0:
            errs.append(abs(np.log10(max(pred_p, 1e-12)) - np.log10(emp_p)))
print(f"    mean |log10(predicted p) - log10(empirical p)| over held-out quantiles: "
      f"{np.mean(errs):.3f}  (n={len(errs)})")
print("    -> the fitted tail reproduces p-values it was not fitted to within "
      f"~{10**np.mean(errs):.1f}x")

# --- check 2: like-for-like correction
done = [r for r in csv.DictReader(open(os.path.join(AN, "29_recombination_maxchi.csv"),
                                       encoding="utf-8-sig")) if r["p_perm"]]
ntrip = json.load(open(os.path.join(AN, "29_recombination_summary.json")))["n_triplets"]
floor = 1.0 / (N_PERM + 1)
for fam, m in (("candidates (n=267)", len(done)), (f"all triplets (n={ntrip})", ntrip)):
    thr = 0.05 / m
    emp = sum(1 for r in done
              if (float(r["p_perm"]) if r["p_method"] == "empirical" else floor) < thr)
    gpd = sum(1 for r in done if float(r["p_perm"]) * m < 0.05)
    print(f"\n  check 2: correction over {fam}, threshold {thr:.2e}")
    print(f"    significant with empirical p (floor {floor:.1e}): {emp}")
    print(f"    significant with GPD-extrapolated p            : {gpd}")

# --- check 3: false positives on column-permuted data
print("\n  check 3: pipeline run on data with spatial structure destroyed")
ncols = A.shape[1]
perm_cols = rng.permutation(ncols)
Ap = A[:, perm_cols]
fp = 0; tested = 0
for r in rec[:40]:
    q, p1, p2 = r["accession"], r["parent1"], r["parent2"]
    iq, a1, a2 = Ap[ixm[q]], Ap[ixm[p1]], Ap[ixm[p2]]
    inf = (a1 != a2) & (a1 != GAP) & (a2 != GAP) & (iq != GAP) & ((iq == a1) | (iq == a2))
    x = (iq[inf] == a1[inf]).astype(np.int64)
    if x.size < 2 * MIN_SIDE or x.sum() in (0, x.size):
        continue
    tested += 1
    obs, _ = maxchi(x)
    nd = null_dist(x, 2000)
    nge = int((nd >= obs).sum())
    if nge > 0:
        p = (1 + nge) / 2001
    else:
        tail = np.sort(nd)[-200:]; u = tail[0]; exc = tail - u
        shp, loc, scl = stats.genpareto.fit(exc[exc > 0], floc=0.0)
        p = (200 / 2000) * float(stats.genpareto.sf(obs - u, shp, loc=loc, scale=scl))
    if p * ntrip < 0.05:
        fp += 1
print(f"    triplets tested on permuted columns: {tested}")
print(f"    calls surviving Bonferroni over {ntrip} triplets: {fp} "
      f"({100*fp/max(tested,1):.1f}%)")
print("    -> with the recombination signal removed the pipeline makes "
      f"{'no' if fp == 0 else str(fp)} call{'s' if fp != 1 else ''}")

json.dump({"window_sensitivity": rows,
           "gpd_quantile_log10_error": float(np.mean(errs)),
           "permuted_false_positives": fp, "permuted_tested": tested},
          open(os.path.join(AN, "39_gpd_and_window_checks.json"), "w"), indent=1)
print("\nwrote 39_gpd_and_window_checks.json")
