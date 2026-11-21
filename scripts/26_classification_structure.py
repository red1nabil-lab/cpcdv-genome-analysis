#!/usr/bin/env python
"""Step 26: is the strain system inconsistent, or merely non-transitive?

Reviewer criticism 2 is the strongest conceptual objection: pairwise identity
thresholds are not transitive (A-B 96%, B-C 96%, A-C 92%), so NO partition can satisfy
every pairwise constraint, and counting violating pairs may therefore measure a property
of thresholds rather than a defect in the published labels.  The objection is correct in
principle, so counting pairs cannot settle the question.

This script answers it with measures that do not depend on the threshold at all:

  1. Violation RATES against the correct denominators, and how many distinct genomes and
     strains are actually involved (a single misplaced genome generates many bad pairs).
  2. Violation MAGNITUDE.  Non-transitivity produces violations hugging the threshold.
     It cannot produce a same-strain pair at 86.7% or a different-strain pair at 99.7%.
  3. SEPARABILITY (threshold-free): for each strain, is the lowest within-strain identity
     greater than the highest identity to any genome outside it?  If yes the strain is a
     well-separated group under ANY threshold; if no, no threshold can recover it.
  4. NEAREST-NEIGHBOUR consistency (threshold-free): is each genome's single most similar
     genome a member of its own published strain?
  5. Clustering: complete linkage at 94% (the strict reading, every pair must comply) and
     single linkage at 94% (the lenient reading), compared against the published labels.
"""
import os, sys, json, csv
import numpy as np
from collections import Counter, defaultdict

sys.stdout.reconfigure(encoding="utf-8")
import os as _os, sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # noqa: E402
from cpcdv_paths import DATA_DIR, FIG_DIR, MUSCLE_BIN  # noqa: E402
AN = DATA_DIR
SRC = sys.argv[1] if len(sys.argv) > 1 else "msa"   # "msa" or "truepw"

order = json.load(open(os.path.join(AN, "06_accession_order.json")))
accs = order["accessions"] if isinstance(order, dict) else order
if SRC == "msa":
    PW = np.load(os.path.join(AN, "06_pid_pw.npy"))
    GP = np.load(os.path.join(AN, "06_pid_gap.npy"))
    label = "MSA-derived"
else:
    meta = json.load(open(os.path.join(AN, "25_truepw_meta.json")))
    accs = meta["accessions"]
    PW = np.load(os.path.join(AN, "25_truepw_pid_pw_A.npy"))
    GP = np.load(os.path.join(AN, "25_truepw_pid_gap_A.npy"))
    label = "true pairwise (scheme A)"
# 06_pid_*.npy store identities as FRACTIONS, 25_truepw_* as percentages: normalise
if np.nanmax(PW) <= 1.5:
    PW = PW * 100.0; GP = GP * 100.0
assert 50.0 < np.nanmax(PW) <= 100.0001, f"unexpected identity scale: {np.nanmax(PW)}"
idx = {a: i for i, a in enumerate(accs)}

md = {r["accession"]: r for r in
      csv.DictReader(open(os.path.join(AN, "03_curated_metadata.csv"), encoding="utf-8-sig"))}
# restrict to the paper's analysis set: the MSA matrix covers all 270 complete genomes,
# including the two BeYDV sequences and the RefSeq duplicate that the paper excludes
KEEP = set(json.load(open(os.path.join(AN, "10_analysis_set.json")))["keep"])
pub = {a: md[a]["strain"] for a in accs
       if a in KEEP and a in md and len(md[a]["strain"]) == 1 and md[a]["strain"].isalpha()}
L = sorted(pub)                      # genomes carrying a published single-letter strain
by = defaultdict(list)
for a in L:
    by[pub[a]].append(a)

print(f"=== classification structure, identities from {label} ===")
print(f"genomes in matrix: {len(accs)}   with a published strain letter: {len(L)}")
print("strain sizes:", dict(sorted(((k, len(v)) for k, v in by.items()),
                                   key=lambda x: -x[1])))

# ---------------------------------------------------------------- 1. rates
win = [(a, b) for k, v in by.items() for i, a in enumerate(v) for b in v[i+1:]]
btw = [(a, b) for i, a in enumerate(L) for b in L[i+1:] if pub[a] != pub[b]]
wv = [(a, b, PW[idx[a], idx[b]]) for a, b in win if PW[idx[a], idx[b]] < 94.0]
bv = [(a, b, PW[idx[a], idx[b]]) for a, b in btw if PW[idx[a], idx[b]] > 94.0]
print(f"\n1. VIOLATION RATES")
print(f"   within-strain pairs : {len(win):6d}   below 94%: {len(wv):4d}  "
      f"({100*len(wv)/len(win):.2f}%)")
print(f"   between-strain pairs: {len(btw):6d}   above 94%: {len(bv):4d}  "
      f"({100*len(bv)/len(btw):.2f}%)")
gw = sorted({a for a, b, _ in wv} | {b for a, b, _ in wv})
gb = sorted({a for a, b, _ in bv} | {b for a, b, _ in bv})
print(f"   distinct genomes in within-strain violations : {len(gw)} "
      f"({100*len(gw)/len(L):.1f}% of labelled genomes)")
print(f"   distinct genomes in between-strain violations: {len(gb)} "
      f"({100*len(gb)/len(L):.1f}%)")
print(f"   strains implicated (within) : {sorted({pub[a] for a in gw})}")
print(f"   strain PAIRS implicated (between): "
      f"{sorted({tuple(sorted((pub[a], pub[b]))) for a, b, _ in bv})}")

# ---------------------------------------------------------------- 2. magnitude
print(f"\n2. VIOLATION MAGNITUDE  (non-transitivity predicts violations hug 94%)")
for name, v, sgn in (("within-strain, below 94%", wv, -1),
                     ("between-strain, above 94%", bv, +1)):
    if not v:
        continue
    dev = sorted(abs(x[2] - 94.0) for x in v)
    big = [x for x in v if abs(x[2] - 94.0) > 2.0]
    print(f"   {name}: n={len(v)}  median deviation {np.median(dev):.2f} pp, "
          f"max {max(dev):.2f} pp")
    print(f"      pairs deviating >2 pp from the threshold: {len(big)} "
          f"({100*len(big)/len(v):.0f}%)  <- cannot be explained by non-transitivity")
    ext = sorted(v, key=lambda x: sgn * -x[2])[:5]
    for a, b, p in ext:
        print(f"      {a} vs {b}  strain {pub[a]}/{pub[b]}  {p:.2f}%")

# ---------------------------------------------------------------- 3. separability
print(f"\n3. SEPARABILITY  (threshold-free: is min(within) > max(outside)?)")
rows = []
for k in sorted(by, key=lambda x: -len(by[x])):
    v = by[k]
    if len(v) < 2:
        continue
    ii = [idx[a] for a in v]
    oo = [idx[a] for a in L if pub[a] != k]
    sub = PW[np.ix_(ii, ii)]
    mn = sub[np.triu_indices(len(ii), 1)].min()
    mx = PW[np.ix_(ii, oo)].max()
    sep = mn > mx
    am = np.unravel_index(PW[np.ix_(ii, oo)].argmax(), (len(ii), len(oo)))
    rows.append((k, len(v), mn, mx, sep, v[am[0]], [a for a in L if pub[a] != k][am[1]]))
    print(f"   strain {k:2s} n={len(v):3d}  min within {mn:6.2f}%  "
          f"max outside {mx:6.2f}%  separable: {'YES' if sep else 'NO '}"
          f"   (closest outsider {[a for a in L if pub[a]!=k][am[1]]}"
          f" strain {pub[[a for a in L if pub[a]!=k][am[1]]]})")
nsep = sum(1 for r in rows if not r[4])
print(f"   -> {nsep} of {len(rows)} multi-genome strains are NOT separable "
      f"under any identity threshold")

# ---------------------------------------------------------------- 4. 1-NN
print(f"\n4. NEAREST-NEIGHBOUR CONSISTENCY (threshold-free)")
ok = bad = 0
badlist = []
for a in L:
    i = idx[a]
    cand = [(PW[i, idx[b]], b) for b in L if b != a]
    p, nb = max(cand)
    if pub[nb] == pub[a]:
        ok += 1
    else:
        bad += 1; badlist.append((a, pub[a], nb, pub[nb], p))
print(f"   genomes whose most similar genome shares their strain letter: "
      f"{ok}/{len(L)} ({100*ok/len(L):.1f}%)")
print(f"   genomes whose nearest neighbour is in a DIFFERENT strain: {bad} "
      f"({100*bad/len(L):.1f}%)")
for a, pa, nb, pb, p in sorted(badlist, key=lambda x: -x[4])[:12]:
    print(f"      {a} (strain {pa}) -> {nb} (strain {pb}) at {p:.2f}%")

# ---------------------------------------------------------------- 5. clustering
def linkage_clusters(members, thr, mode):
    """complete linkage = every pair within a cluster must meet thr;
       single linkage  = a chain of >=thr links suffices."""
    n = len(members)
    ii = [idx[a] for a in members]
    M = PW[np.ix_(ii, ii)]
    if mode == "single":
        parent = list(range(n))
        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]; x = parent[x]
            return x
        for i in range(n):
            for j in range(i+1, n):
                if M[i, j] >= thr:
                    a, b = find(i), find(j)
                    if a != b:
                        parent[a] = b
        g = defaultdict(list)
        for i in range(n):
            g[find(i)].append(members[i])
        return list(g.values())
    # complete linkage, agglomerative on the minimum-identity criterion
    cl = [[i] for i in range(n)]
    while True:
        best = None
        for x in range(len(cl)):
            for y in range(x+1, len(cl)):
                m = M[np.ix_(cl[x], cl[y])].min()
                if m >= thr and (best is None or m > best[0]):
                    best = (m, x, y)
        if best is None:
            break
        _, x, y = best
        cl[x] = cl[x] + cl[y]; cl.pop(y)
    return [[members[i] for i in c] for c in cl]


def ari(lab1, lab2):
    from math import comb
    keys = list(lab1)
    c = Counter((lab1[k], lab2[k]) for k in keys)
    a = Counter(lab1[k] for k in keys); b = Counter(lab2[k] for k in keys)
    sij = sum(comb(v, 2) for v in c.values())
    sa = sum(comb(v, 2) for v in a.values()); sb = sum(comb(v, 2) for v in b.values())
    n2 = comb(len(keys), 2)
    exp = sa * sb / n2
    mx = 0.5 * (sa + sb)
    return (sij - exp) / (mx - exp) if mx != exp else 1.0


print(f"\n5. CLUSTERING AT 94% vs THE PUBLISHED LABELS")
for mode in ("complete", "single"):
    cl = linkage_clusters(L, 94.0, mode)
    cl.sort(key=len, reverse=True)
    assign = {a: f"c{i}" for i, c in enumerate(cl) for a in c}
    print(f"   {mode} linkage: {len(cl)} clusters "
          f"(published strains among these genomes: {len(by)})  "
          f"adjusted Rand index vs published = {ari(pub, assign):.3f}")
    # how do published strains map onto clusters?
    split = merged = 0
    for k, v in by.items():
        if len({assign[a] for a in v}) > 1:
            split += 1
    for i, c in enumerate(cl):
        if len({pub[a] for a in c}) > 1:
            merged += 1
    print(f"      published strains split across >1 cluster: {split}/{len(by)}")
    print(f"      clusters containing >1 published strain   : {merged}/{len(cl)}")
    if mode == "complete":
        for i, c in enumerate(cl[:8]):
            print(f"      cluster {i} (n={len(c)}): "
                  f"{dict(Counter(pub[a] for a in c))}")

out = os.path.join(AN, f"26_classification_structure_{SRC}.json")
json.dump({"source": label, "n_labelled": len(L),
           "within_pairs": len(win), "within_violations": len(wv),
           "between_pairs": len(btw), "between_violations": len(bv),
           "genomes_in_within_violations": gw,
           "genomes_in_between_violations": gb,
           "nn_consistent": ok, "nn_inconsistent": bad,
           "nn_bad": [[a, pa, nb, pb, round(p, 2)] for a, pa, nb, pb, p in badlist],
           "separability": [[k, n, round(mn, 2), round(mx, 2), bool(s)]
                            for k, n, mn, mx, s, _, _ in rows]},
          open(out, "w"), indent=1)
print(f"\nwrote {out}")
