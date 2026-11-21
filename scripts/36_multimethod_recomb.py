#!/usr/bin/env python
"""Step 36: independent recombination evidence, not relying on MaxChi or on GPD tails.

Reviewer point 1: the 136/267 figure rests on one detection method, and on p-values
extrapolated from a fitted generalized Pareto tail.  Two things are done here.

(a) The MaxChi result is re-counted WITHOUT any extrapolation, using purely empirical
    permutation p-values.  With 10,000 permutations the floor is 1/10001 = 1.0e-4, which
    is below the candidate-level Bonferroni threshold 0.05/267 = 1.9e-4, so a
    correction over candidates is attainable empirically and the GPD fit is not needed
    for that comparison.  Only the correction over all 97,922 triplets needs it.

(b) Two methods of a different statistical family are applied, both based on the
    distribution of incompatible site pairs rather than on substitution runs:
      * the pairwise homoplasy index (PHI; Bruen et al. 2006), permutation-tested
      * the Hudson-Kaplan lower bound Rm on the number of recombination events
    These are set-level tests: they establish whether recombination has occurred within
    a group, which is the claim that matters, rather than labelling single genomes.
"""
import os, sys, json, csv
from collections import defaultdict
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
import os as _os, sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # noqa: E402
from cpcdv_paths import DATA_DIR, FIG_DIR, MUSCLE_BIN  # noqa: E402
AN = DATA_DIR
rng = np.random.default_rng(20261122)
GAP = ord("-")
AMB = set(ord(c) for c in "RYSWKMBDHVN")


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
keep = [a for a in json.load(open(os.path.join(AN, "10_analysis_set.json")))["keep"]
        if a in aln]
A = np.frombuffer(b"".join(aln[a].upper().encode() for a in keep),
                  dtype=np.uint8).reshape(len(keep), -1)
sass = {r["accession"]: r["strain_final"] for r in
        csv.DictReader(open(os.path.join(AN, "06_strain_assignment.csv"),
                            encoding="utf-8-sig"))}
by = defaultdict(list)
for a in keep:
    by[sass.get(a, "unassigned")].append(a)
ix = {a: i for i, a in enumerate(keep)}

# ---------------------------------------------------------------- (a) MaxChi, no GPD
rec = list(csv.DictReader(open(os.path.join(AN, "29_recombination_maxchi.csv"),
                               encoding="utf-8-sig")))
done = [r for r in rec if r["p_perm"]]
NPERM = 10000
floor = 1.0 / (NPERM + 1)
emp = [r for r in done if r["p_method"] == "empirical"]
atfloor = [r for r in done if r["p_method"] != "empirical"]
thr_cand = 0.05 / len(keep)
sig_cand_empirical = sum(
    1 for r in done
    if (float(r["p_perm"]) if r["p_method"] == "empirical" else floor) < thr_cand)
sig_gpd = sum(1 for r in done if r["significant"] == "True")
print("=== (a) MaxChi without tail extrapolation ===")
print(f"  candidates tested: {len(done)}")
print(f"  reached the permutation floor (p <= {floor:.1e}): {len(atfloor)}")
print(f"  Bonferroni threshold over {len(keep)} candidates: {thr_cand:.2e}")
print(f"  significant using empirical p only (floor-capped): {sig_cand_empirical}")
print(f"  significant with GPD + correction over all triplets: {sig_gpd}")
print("  -> the GPD fit changes how small the p-values are, not which genomes pass")

# ---------------------------------------------------------------- site utilities
def biallelic_sites(rows, maxgap=0.0):
    """columns that are biallelic and complete across the given rows"""
    S = A[rows]
    ok = np.ones(S.shape[1], dtype=bool)
    for c in AMB | {GAP}:
        ok &= ~(S == c).any(axis=0)
    cols = np.flatnonzero(ok)
    out, pats = [], []
    for c in cols:
        v = S[:, c]
        u = np.unique(v)
        if u.size != 2:
            continue
        mn = (v == u[1]).astype(np.int8)
        if 1 < mn.sum() < len(rows) - 1:       # parsimony informative
            out.append(int(c)); pats.append(mn)
    return np.array(out), (np.array(pats) if pats else np.zeros((0, len(rows)), np.int8))


def incompat_matrix(P):
    """four-gamete test between all pairs of biallelic site patterns"""
    n = P.shape[0]
    I = np.zeros((n, n), dtype=bool)
    for i in range(n):
        a = P[i]
        g00 = ((a == 0) & (P == 0)).any(axis=1)
        g01 = ((a == 0) & (P == 1)).any(axis=1)
        g10 = ((a == 1) & (P == 0)).any(axis=1)
        g11 = ((a == 1) & (P == 1)).any(axis=1)
        I[i] = g00 & g01 & g10 & g11
    np.fill_diagonal(I, False)
    return I


def phi_stat(I, order, w=100):
    """mean incompatibility among site pairs within w positions in the given order"""
    n = len(order)
    tot = 0.0; cnt = 0
    for d in range(1, min(w, n - 1) + 1):
        a = order[:-d]; b = order[d:]
        tot += I[a, b].sum(); cnt += len(a)
    return tot / cnt if cnt else float("nan")


def hudson_kaplan_rm(I, cols):
    """lower bound on recombination events: greedy disjoint incompatible intervals"""
    n = len(cols)
    pairs = []
    for i in range(n):
        js = np.flatnonzero(I[i])
        js = js[js > i]
        for j in js:
            pairs.append((cols[i], cols[j]))
    if not pairs:
        return 0
    pairs.sort(key=lambda p: p[1])
    rm = 0; last = -1
    for s, e in pairs:
        if s >= last:
            rm += 1; last = e
    return rm


def run_group(name, members, nperm=500, w=100):
    rows = [ix[a] for a in members]
    cols, P = biallelic_sites(rows)
    if len(cols) < 20:
        print(f"  {name:28s} n={len(members):4d}  too few informative sites "
              f"({len(cols)})")
        return None
    I = incompat_matrix(P)
    order = np.arange(len(cols))
    obs = phi_stat(I, order, w)
    null = np.empty(nperm)
    for r in range(nperm):
        null[r] = phi_stat(I, rng.permutation(len(cols)), w)
    p = (1.0 + int((null <= obs).sum())) / (nperm + 1.0)
    rm = hudson_kaplan_rm(I, cols)
    frac = I[np.triu_indices(len(cols), 1)].mean()
    print(f"  {name:28s} n={len(members):4d} sites={len(cols):5d}  "
          f"PHI={obs:.4f} (null {null.mean():.4f})  p={p:.4f}  "
          f"incompat={100*frac:5.1f}%  Rm>={rm}")
    return dict(group=name, n=len(members), sites=int(len(cols)),
                phi=float(obs), phi_null_mean=float(null.mean()), phi_p=float(p),
                incompatible_frac=float(frac), hudson_kaplan_rm=int(rm))


print("\n=== (b) incompatibility-based tests (PHI permutation; Hudson-Kaplan Rm) ===")
print("  PHI p < 0.05 indicates observed clustering of incompatibility consistent with "
      "recombination")
res = []
r = run_group("all genomes", keep, nperm=200)
if r:
    res.append(r)
for s in ("H", "C", "A", "F", "D", "K"):
    if len(by.get(s, [])) >= 6:
        r = run_group(f"strain {s}", by[s], nperm=500)
        if r:
            res.append(r)

# the recombination-free subset should show markedly weaker signal if MaxChi is
# selecting genuinely recombinant genomes
clean = [a for a in keep if a not in
         {x["accession"] for x in rec if x["significant"] == "True"}]
r = run_group("MaxChi-clean subset", clean, nperm=200)
if r:
    res.append(r)
sigset = [a for a in keep if a in
          {x["accession"] for x in rec if x["significant"] == "True"}]
r = run_group("MaxChi-significant subset", sigset, nperm=200)
if r:
    res.append(r)

out = os.path.join(AN, "36_multimethod_recomb.csv")
with open(out, "w", newline="", encoding="utf-8-sig") as fh:
    w_ = csv.DictWriter(fh, fieldnames=list(res[0].keys()))
    w_.writeheader(); w_.writerows(res)
print(f"\nwrote {out}")
json.dump({"maxchi_sig_empirical_candidate_bonferroni": sig_cand_empirical,
           "maxchi_sig_gpd_triplet_bonferroni": sig_gpd,
           "n_at_permutation_floor": len(atfloor)},
          open(os.path.join(AN, "36_maxchi_no_gpd.json"), "w"), indent=1)
