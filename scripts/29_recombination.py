#!/usr/bin/env python
"""Step 29: formal recombination detection by MaxChi, with breakpoints and p-values.

Reviewer criticism 4: the manuscript inferred mosaicism from capsid-vs-Rep nearest
neighbour discordance, which identifies candidates but establishes nothing - no
breakpoints, no parents, no statistical support.  RDP5 is the standard tool but is an
interactive Windows application and is not available here, so we implement the maximum
chi-square method (Maynard Smith 1992), which is one of the methods RDP5 itself runs and
which is fully specified and testable.

For each candidate genome Q and each pair of potential parents (P1, P2):
  * take the informative sites, where P1 != P2 and Q matches one of them
  * code each site +1 if Q matches P1, 0 if it matches P2
  * scan every breakpoint, computing the 2x2 chi-square comparing the proportion of
    P1-matching sites left and right of it; keep the maximum
  * assess significance by permuting the order of the informative sites, which destroys
    any spatial clustering while preserving the overall match proportions

Two stages, as is standard: stage 1 scores every triplet; stage 2 runs the permutation
test only on triplets that pass a screening threshold, then applies a Bonferroni
correction for the total number of triplets examined.
"""
import os, sys, json, csv, time
from collections import defaultdict
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
import os as _os, sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # noqa: E402
from cpcdv_paths import DATA_DIR, FIG_DIR, MUSCLE_BIN  # noqa: E402
AN = DATA_DIR
rng = np.random.default_rng(20261120)

MIN_SIDE = 10          # minimum informative sites either side of a breakpoint
N_PERM = 10000
TOP_FOR_PERM = 600     # triplets taken forward to the permutation stage
GAP = ord("-")
# A plain permutation p-value cannot fall below 1/(N_PERM+1), so with a Bonferroni
# factor of ~1e5 no result could ever reach significance however strong it was -- which
# is exactly what happened on the first run (260 of 267 candidates pinned at the floor
# with chi-square up to 266).  The tail of the permutation null is therefore fitted with
# a Generalized Pareto distribution, which extrapolates p-values far below the floor
# (Knijnenburg et al. 2009, Bioinformatics 25, i161-i168).
GPD_EXC = 250          # number of largest permutation values used for the tail fit


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


aln = read_fasta(os.path.join(AN, "05_full_aln.fasta"))
KEEP = list(json.load(open(os.path.join(AN, "10_analysis_set.json")))["keep"])
KEEP = [a for a in KEEP if a in aln]
A = np.frombuffer(b"".join(aln[a].upper().encode() for a in KEEP),
                  dtype=np.uint8).reshape(len(KEEP), -1)
idx = {a: i for i, a in enumerate(KEEP)}
NC = A.shape[1]
print(f"alignment: {len(KEEP)} genomes x {NC} columns")

md = {r["accession"]: r for r in
      csv.DictReader(open(os.path.join(AN, "03_curated_metadata.csv"),
                          encoding="utf-8-sig"))}
sass = {r["accession"]: r["strain_final"] for r in
        csv.DictReader(open(os.path.join(AN, "06_strain_assignment.csv"),
                            encoding="utf-8-sig"))}

# potential parents: one representative per strain, plus each candidate's nearest
# neighbours (so that a true parent outside the representative panel is still reachable)
PW = np.load(os.path.join(AN, "06_pid_pw.npy")).astype(float)
msa_accs = json.load(open(os.path.join(AN, "06_accession_order.json")))
mi = {a: i for i, a in enumerate(msa_accs)}
if np.nanmax(PW) <= 1.5:
    PW = PW * 100.0

by_strain = defaultdict(list)
for a in KEEP:
    by_strain[sass.get(a, "unassigned")].append(a)
reps = []
for s, v in sorted(by_strain.items()):
    if s == "unassigned":
        continue
    reps.append(sorted(v)[0])
print(f"strain representatives used as candidate parents: {len(reps)}")


def maxchi(x):
    """max 2x2 chi-square over all breakpoints of a 0/1 site series."""
    n = x.size
    if n < 2 * MIN_SIDE:
        return 0.0, -1
    c = np.cumsum(x)
    tot = c[-1]
    k = np.arange(MIN_SIDE, n - MIN_SIDE + 1)
    a = c[k - 1]                      # P1-matches left of breakpoint
    b = tot - a                       # P1-matches right
    left, right = k, n - k
    num = n * (a * (right - b) - (left - a) * b) ** 2
    den = left * right * tot * (n - tot)
    with np.errstate(divide="ignore", invalid="ignore"):
        chi = np.where(den > 0, num / den, 0.0)
    j = int(np.nanargmax(chi))
    return float(chi[j]), int(k[j])


def maxchi_perm_dist(x, nperm):
    """null distribution of max chi-square under random site order."""
    n = x.size
    # vectorised permutation: argsort of a random matrix gives nperm independent orders
    order = rng.random((nperm, n)).argsort(axis=1)
    P = x[order]
    C = np.cumsum(P, axis=1)
    tot = C[:, -1]
    k = np.arange(MIN_SIDE, n - MIN_SIDE + 1)
    a = C[:, k - 1]
    b = tot[:, None] - a
    left, right = k[None, :], (n - k)[None, :]
    num = n * (a * (right - b) - (left - a) * b) ** 2
    den = left * right * tot[:, None] * (n - tot)[:, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        chi = np.where(den > 0, num / den, 0.0)
    out = np.nanmax(chi, axis=1)
    return out


# ------------------------------------------------------------------ stage 1
t0 = time.time()
cand_rows = []
ntrip = 0
for q in KEEP:
    iq = A[idx[q]]
    nn = sorted(((PW[mi[q], mi[b]], b) for b in KEEP if b != q and b in mi),
                reverse=True)[:10]
    parents = sorted({b for _, b in nn} | set(reps) - {q})
    best = None
    for i1 in range(len(parents)):
        for i2 in range(i1 + 1, len(parents)):
            p1, p2 = parents[i1], parents[i2]
            a1, a2 = A[idx[p1]], A[idx[p2]]
            inf = (a1 != a2) & (a1 != GAP) & (a2 != GAP) & (iq != GAP) \
                  & ((iq == a1) | (iq == a2))
            n = int(inf.sum())
            if n < 2 * MIN_SIDE:
                continue
            ntrip += 1
            x = (iq[inf] == a1[inf]).astype(np.int64)
            if x.sum() == 0 or x.sum() == n:
                continue
            chi, kk = maxchi(x)
            if best is None or chi > best[0]:
                cols = np.flatnonzero(inf)
                bp = int(cols[min(kk, n - 1)])
                best = (chi, p1, p2, n, bp, kk)
    if best:
        chi, p1, p2, n, bp, kk = best
        cand_rows.append(dict(accession=q, strain=sass.get(q, "unassigned"),
                              country=md.get(q, {}).get("country", ""),
                              parent1=p1, parent2=p2,
                              parent1_strain=sass.get(p1, "?"),
                              parent2_strain=sass.get(p2, "?"),
                              informative_sites=n, maxchi=round(chi, 3),
                              breakpoint_aln_col=bp))
print(f"stage 1: {ntrip} triplets scored for {len(cand_rows)} candidates "
      f"in {time.time()-t0:.0f}s")

# ------------------------------------------------------------------ stage 2
cand_rows.sort(key=lambda r: -r["maxchi"])
take = cand_rows[:TOP_FOR_PERM]
print(f"stage 2: permutation test ({N_PERM} permutations) on top {len(take)} candidates")
for r in take:
    q, p1, p2 = r["accession"], r["parent1"], r["parent2"]
    iq, a1, a2 = A[idx[q]], A[idx[p1]], A[idx[p2]]
    inf = (a1 != a2) & (a1 != GAP) & (a2 != GAP) & (iq != GAP) \
          & ((iq == a1) | (iq == a2))
    x = (iq[inf] == a1[inf]).astype(np.int64)
    null = maxchi_perm_dist(x, N_PERM)
    obs = r["maxchi"]
    n_ge = int(np.sum(null >= obs))
    if n_ge > 0:
        p = (1.0 + n_ge) / (N_PERM + 1.0)
        r["p_method"] = "empirical"
    else:
        # observed exceeds every permutation: fit a GPD to the null tail
        from scipy import stats
        tail = np.sort(null)[-GPD_EXC:]
        u = tail[0]
        exc = tail - u
        try:
            shp, loc, scl = stats.genpareto.fit(exc[exc > 0], floc=0.0)
            p = (GPD_EXC / N_PERM) * float(
                stats.genpareto.sf(obs - u, shp, loc=loc, scale=scl))
            r["p_method"] = "GPD tail fit"
            if not np.isfinite(p) or p <= 0:
                p = 1.0 / (N_PERM + 1.0); r["p_method"] = "floor (GPD failed)"
        except Exception:
            p = 1.0 / (N_PERM + 1.0); r["p_method"] = "floor (GPD failed)"
    r["p_perm"] = p
    r["p_bonferroni"] = min(1.0, p * ntrip)      # corrected over ALL triplets searched
    r["p_bonf_candidates"] = min(1.0, p * len(KEEP))
    r["significant"] = bool(r["p_bonferroni"] < 0.05)

sig = [r for r in take if r.get("significant")]
print(f"\ncandidates with Bonferroni-corrected p < 0.05: {len(sig)} of {len(KEEP)} "
      f"({100*len(sig)/len(KEEP):.1f}%)")

disc = {r["accession"] for r in
        csv.DictReader(open(os.path.join(AN, "11_strain_assignment_discordance.csv"),
                            encoding="utf-8-sig"))}
hit = {r["accession"] for r in sig}
print(f"previously reported capsid/Rep discordant genomes: {len(disc)}")
print(f"  of those, MaxChi-significant: {len(disc & hit)}/{len(disc)} "
      f"{sorted(disc & hit)}")
print(f"  discordant but NOT MaxChi-significant: {sorted(disc - hit)}")
print(f"  MaxChi-significant but not previously flagged: {len(hit - disc)}")

print("\ntop 20 supported recombinants:")
for r in sig[:20]:
    print(f"  {r['accession']} ({r['strain']:>10s}, {r['country'][:12]:12s}) "
          f"parents {r['parent1']}[{r['parent1_strain']}] x "
          f"{r['parent2']}[{r['parent2_strain']}]  "
          f"maxchi={r['maxchi']:7.1f}  bp~col {r['breakpoint_aln_col']:4d}  "
          f"p_bonf={r['p_bonferroni']:.2e} [{r['p_method']}]")

out = os.path.join(AN, "29_recombination_maxchi.csv")
flds = ["accession", "strain", "country", "parent1", "parent1_strain", "parent2",
        "parent2_strain", "informative_sites", "maxchi", "breakpoint_aln_col",
        "p_perm", "p_method", "p_bonferroni", "p_bonf_candidates", "significant"]
with open(out, "w", newline="", encoding="utf-8-sig") as fh:
    w = csv.DictWriter(fh, fieldnames=flds, extrasaction="ignore")
    w.writeheader()
    for r in cand_rows:
        w.writerow({k: r.get(k, "") for k in flds})
print(f"\nwrote {out}")
json.dump({"n_triplets": ntrip, "n_significant": len(sig),
           "significant_accessions": sorted(hit), "n_perm": N_PERM,
           "min_side": MIN_SIDE},
          open(os.path.join(AN, "29_recombination_summary.json"), "w"), indent=1)
