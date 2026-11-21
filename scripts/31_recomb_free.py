#!/usr/bin/env python
"""Step 31: do the strain-threshold failures survive removal of the recombinants?

Reviewer criticism 4, final bullet.  Step 29 found statistically supported recombination
signal in 136 of 267 genomes.  If the threshold violations and the non-separability are
simply consequences of recombination, they should largely disappear in the
recombination-free subset.  If they persist, they are properties of the designations.

Also maps the inferred breakpoints onto the genome to see whether they cluster.
"""
import os, sys, json, csv
from collections import Counter, defaultdict
from math import comb
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
import os as _os, sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # noqa: E402
from cpcdv_paths import DATA_DIR, FIG_DIR, MUSCLE_BIN  # noqa: E402
AN = DATA_DIR

md = {r["accession"]: r for r in
      csv.DictReader(open(os.path.join(AN, "03_curated_metadata.csv"),
                          encoding="utf-8-sig"))}
rec = list(csv.DictReader(open(os.path.join(AN, "29_recombination_maxchi.csv"),
                               encoding="utf-8-sig")))
SIG = {r["accession"] for r in rec if r["significant"] == "True"}
KEEP = list(json.load(open(os.path.join(AN, "10_analysis_set.json")))["keep"])
CLEAN = [a for a in KEEP if a not in SIG]
print(f"analysis set {len(KEEP)}; recombination-signal genomes {len(SIG)}; "
      f"clean subset {len(CLEAN)}")

accs = KEEP
PW = np.load(os.path.join(AN, "25_truepw_pid_pw_A.npy")).astype(float)
if np.nanmax(PW) <= 1.5:
    PW *= 100.0
idx = {a: i for i, a in enumerate(accs)}


def ari(l1, l2):
    ks = list(l1)
    c = Counter((l1[k], l2[k]) for k in ks)
    a = Counter(l1[k] for k in ks); b = Counter(l2[k] for k in ks)
    sij = sum(comb(v, 2) for v in c.values())
    sa = sum(comb(v, 2) for v in a.values()); sb = sum(comb(v, 2) for v in b.values())
    n2 = comb(len(ks), 2); e = sa * sb / n2; m = .5 * (sa + sb)
    return (sij - e) / (m - e) if m != e else 1.0


def complete_linkage(members, thr=94.0):
    n = len(members); ii = [idx[a] for a in members]
    S = PW[np.ix_(ii, ii)]
    cl = [[i] for i in range(n)]
    while True:
        best = None
        for x in range(len(cl)):
            for y in range(x + 1, len(cl)):
                m = S[np.ix_(cl[x], cl[y])].min()
                if m >= thr and (best is None or m > best[0]):
                    best = (m, x, y)
        if best is None:
            break
        _, x, y = best
        cl[x] += cl[y]; cl.pop(y)
    return [[members[i] for i in c] for c in cl]


def report(name, subset):
    pub = {a: md[a]["strain"] for a in subset
           if len(md[a]["strain"]) == 1 and md[a]["strain"].isalpha()}
    L = sorted(pub)
    by = defaultdict(list)
    for a in L:
        by[pub[a]].append(a)
    multi = {k: v for k, v in by.items() if len(v) > 1}
    win = [(a, b) for v in by.values() for i, a in enumerate(v) for b in v[i+1:]]
    btw = [(a, b) for i, a in enumerate(L) for b in L[i+1:] if pub[a] != pub[b]]
    wv = [p for p in win if PW[idx[p[0]], idx[p[1]]] < 94.0]
    bv = [p for p in btw if PW[idx[p[0]], idx[p[1]]] > 94.0]
    nonsep = []
    for k, v in multi.items():
        ii = [idx[a] for a in v]; oo = [idx[a] for a in L if pub[a] != k]
        if not oo:
            continue
        mn = PW[np.ix_(ii, ii)][np.triu_indices(len(ii), 1)].min()
        mx = PW[np.ix_(ii, oo)].max()
        if mn <= mx:
            nonsep.append(k)
    nnbad = 0
    for a in L:
        i = idx[a]
        _, nb = max((PW[i, idx[b]], b) for b in L if b != a)
        if pub[nb] != pub[a]:
            nnbad += 1
    cl = complete_linkage(L)
    assign = {a: f"c{i}" for i, c in enumerate(cl) for a in c}
    print(f"\n=== {name} ===")
    print(f"  labelled genomes {len(L)}, multi-genome strains {len(multi)} "
          f"{dict(sorted(((k, len(v)) for k, v in multi.items()), key=lambda x:-x[1]))}")
    print(f"  within-strain pairs below 94% : {len(wv):4d}/{len(win)} "
          f"({100*len(wv)/max(len(win),1):.2f}%)")
    print(f"  between-strain pairs above 94%: {len(bv):4d}/{len(btw)} "
          f"({100*len(bv)/max(len(btw),1):.2f}%)")
    print(f"  strains NOT separable (threshold-free): {len(nonsep)}/{len(multi)} "
          f"{sorted(nonsep)}")
    print(f"  nearest neighbour in another strain: {nnbad}/{len(L)} "
          f"({100*nnbad/max(len(L),1):.1f}%)")
    print(f"  complete-linkage 94% clusters {len(cl)}; ARI vs published "
          f"{ari(pub, assign):.3f}")
    return dict(subset=name, n=len(L), within=len(win), within_bad=len(wv),
                between=len(btw), between_bad=len(bv), nonsep=len(nonsep),
                nonsep_strains=";".join(sorted(nonsep)), nn_bad=nnbad,
                clusters=len(cl), ari=round(ari(pub, assign), 4))


rows = [report("all 267 genomes", KEEP),
        report(f"recombination-free subset (n={len(CLEAN)})", CLEAN)]

# ---------------------------------------------------------------- breakpoints
print("\n=== inferred breakpoint positions (significant candidates only) ===")
REGIONS = [("LIR", 0, 181), ("V1 / MP", 181, 507), ("(V1-V2 spacer)", 507, 519),
           ("V2 / CP", 519, 1264), ("SIR", 1264, 1437), ("Rep exon 1", 1437, 1859),
           ("(intron)", 1859, 1945), ("Rep exon 2", 1945, 2632), ("LIR (3')", 2632, 2829)]
bps = [int(r["breakpoint_aln_col"]) for r in rec if r["significant"] == "True"]
cnt = Counter()
for b in bps:
    for nm, s, e in REGIONS:
        if s <= b < e:
            cnt[nm] += 1
            break
tot = sum(cnt.values())
for nm, s, e in REGIONS:
    n = cnt.get(nm, 0)
    width = e - s
    print(f"  {nm:16s} cols {s:4d}-{e:4d} ({width:4d} nt)  breakpoints {n:4d} "
          f"({100*n/tot:5.1f}%)  density {n/width*1000:6.2f} per kb")
print(f"\n  most frequent individual breakpoint columns: "
      f"{Counter(bps).most_common(8)}")

out = os.path.join(AN, "31_recomb_free_comparison.csv")
with open(out, "w", newline="", encoding="utf-8-sig") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
print(f"\nwrote {out}")
json.dump({"significant_recombinants": sorted(SIG), "clean_subset": CLEAN,
           "breakpoint_region_counts": dict(cnt)},
          open(os.path.join(AN, "31_recomb_free_sets.json"), "w"), indent=1)
