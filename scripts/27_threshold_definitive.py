#!/usr/bin/env python
"""Step 27: the definitive threshold evaluation, under every method x gap convention.

This settles three reviewer criticisms at once and corrects an error of our own.

Criticism 1 (ICTV procedure): identity should come from true pairwise alignment with
pairwise deletion of gaps, not from a multiple sequence alignment.
Criticism 2 (non-transitivity): counting violating pairs may measure a property of
thresholds rather than a defect in the labels.
Criticism 3 (species boundary): the L-R value sits 0.13 pp below 78% and so is
hostage to alignment method.

Our own error, found while checking: the previously reported counts (41 within-strain
pairs below 94%, 65 between-strain above) were computed from the gaps-as-differences
matrix, while the manuscript described them as pairwise deletion.  Both conventions are
computed here and labelled correctly.

Four cells:  {MSA, true pairwise} x {pairwise deletion, gaps as differences}
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

KEEP = set(json.load(open(os.path.join(AN, "10_analysis_set.json")))["keep"])
md = {r["accession"]: r for r in
      csv.DictReader(open(os.path.join(AN, "03_curated_metadata.csv"),
                          encoding="utf-8-sig"))}

msa_accs = json.load(open(os.path.join(AN, "06_accession_order.json")))
# 25_truepw_meta.json is only written once both scoring schemes finish; the row order of
# the scheme-A matrix is the analysis-set order that step 25 iterated over
_mp = os.path.join(AN, "25_truepw_meta.json")
if os.path.exists(_mp):
    tpw_accs = json.load(open(_mp))["accessions"]
else:
    tpw_accs = list(json.load(open(os.path.join(AN, "10_analysis_set.json")))["keep"])
    print(f"(25_truepw_meta.json not yet written; using analysis-set order, "
          f"n={len(tpw_accs)})")

SOURCES = {}
for nm, accs, f_pw, f_gp in (
        ("MSA", msa_accs, "06_pid_pw.npy", "06_pid_gap.npy"),
        ("pairwise", tpw_accs, "25_truepw_pid_pw_A.npy", "25_truepw_pid_gap_A.npy")):
    pw = np.load(os.path.join(AN, f_pw)).astype(float)
    gp = np.load(os.path.join(AN, f_gp)).astype(float)
    if np.nanmax(pw) <= 1.5:
        pw = pw * 100.0; gp = gp * 100.0
    SOURCES[nm] = (accs, {"pairwise deletion": pw, "gaps as differences": gp})

STRAIN_REPS = {}   # strain -> representative accession (as used in step 8)
for r in csv.DictReader(open(os.path.join(AN, "08_species_identity.csv"),
                             encoding="utf-8-sig")):
    for c in (r["seq_a"], r["seq_b"]):
        if c.startswith("CpCDV_strain"):
            s, acc = c[len("CpCDV_strain"):].split("_", 1)
            STRAIN_REPS[s] = acc


def ari(lab1, lab2):
    keys = list(lab1)
    c = Counter((lab1[k], lab2[k]) for k in keys)
    a = Counter(lab1[k] for k in keys); b = Counter(lab2[k] for k in keys)
    sij = sum(comb(v, 2) for v in c.values())
    sa = sum(comb(v, 2) for v in a.values()); sb = sum(comb(v, 2) for v in b.values())
    n2 = comb(len(keys), 2); exp = sa * sb / n2; mx = 0.5 * (sa + sb)
    return (sij - exp) / (mx - exp) if mx != exp else 1.0


def complete_linkage(members, M, idx, thr):
    n = len(members)
    ii = [idx[a] for a in members]
    S = M[np.ix_(ii, ii)]
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
        cl[x] = cl[x] + cl[y]; cl.pop(y)
    return [[members[i] for i in c] for c in cl]


summary = []
print("=" * 86)
for src, (accs, mats) in SOURCES.items():
    idx = {a: i for i, a in enumerate(accs)}
    pub = {a: md[a]["strain"] for a in accs
           if a in KEEP and a in md and len(md[a]["strain"]) == 1
           and md[a]["strain"].isalpha()}
    L = sorted(pub)
    by = defaultdict(list)
    for a in L:
        by[pub[a]].append(a)
    multi = {k: v for k, v in by.items() if len(v) > 1}
    for conv, M in mats.items():
        win = [(a, b) for v in by.values() for i, a in enumerate(v) for b in v[i+1:]]
        btw = [(a, b) for i, a in enumerate(L) for b in L[i+1:] if pub[a] != pub[b]]
        wv = [(a, b, M[idx[a], idx[b]]) for a, b in win if M[idx[a], idx[b]] < 94.0]
        bv = [(a, b, M[idx[a], idx[b]]) for a, b in btw if M[idx[a], idx[b]] > 94.0]
        # separability
        nonsep = []
        for k, v in multi.items():
            ii = [idx[a] for a in v]; oo = [idx[a] for a in L if pub[a] != k]
            mn = M[np.ix_(ii, ii)][np.triu_indices(len(ii), 1)].min()
            mx = M[np.ix_(ii, oo)].max()
            if mn <= mx:
                nonsep.append(k)
        # 1-NN
        nnbad = 0
        for a in L:
            i = idx[a]
            p, nb = max((M[i, idx[b]], b) for b in L if b != a)
            if pub[nb] != pub[a]:
                nnbad += 1
        # clustering
        cl = complete_linkage(L, M, idx, 94.0)
        assign = {a: f"c{i}" for i, c in enumerate(cl) for a in c}
        a_ri = ari(pub, assign)
        # key pairs
        lr = (M[idx[STRAIN_REPS['L']], idx[STRAIN_REPS['R']]]
              if STRAIN_REPS.get('L') in idx and STRAIN_REPS.get('R') in idx else float("nan"))
        reps = [s for s in STRAIN_REPS if STRAIN_REPS[s] in idx]
        rp = [(M[idx[STRAIN_REPS[x]], idx[STRAIN_REPS[y]]], x, y)
              for i, x in enumerate(sorted(reps)) for y in sorted(reps)[i+1:]]
        below78 = sum(1 for v, _, _ in rp if v < 78.0)
        cl_max = max((M[idx[a], idx[b]] for a, b in btw
                      if {pub[a], pub[b]} == {"C", "L"}), default=float("nan"))

        print(f"\n--- {src} alignment, {conv} ---")
        print(f"  within-strain pairs below 94% : {len(wv):4d} / {len(win)} "
              f"({100*len(wv)/len(win):.2f}%)")
        print(f"  between-strain pairs above 94%: {len(bv):4d} / {len(btw)} "
              f"({100*len(bv)/len(btw):.2f}%)")
        if bv:
            dev = [abs(v - 94.0) for _, _, v in bv]
            print(f"     between-strain: median deviation {np.median(dev):.2f} pp, "
                  f">2 pp: {sum(1 for d in dev if d > 2)}")
        if wv:
            dev = [abs(v - 94.0) for _, _, v in wv]
            print(f"     within-strain : median deviation {np.median(dev):.2f} pp, "
                  f">2 pp: {sum(1 for d in dev if d > 2)}")
        print(f"  strains not separable (threshold-free): {len(nonsep)}/{len(multi)} "
              f"{sorted(nonsep)}")
        print(f"  genomes whose nearest neighbour is another strain: {nnbad}/{len(L)} "
              f"({100*nnbad/len(L):.1f}%)")
        print(f"  complete-linkage 94% clusters: {len(cl)}  ARI vs published = {a_ri:.3f}")
        print(f"  strain L vs strain R reference identity: {lr:.2f}%   "
              f"(rep pairs below 78%: {below78}/{len(rp)})")
        print(f"  max identity between a strain-C and a strain-L genome: {cl_max:.2f}%")
        summary.append(dict(alignment=src, gap_convention=conv,
                            within_pairs=len(win), within_below94=len(wv),
                            between_pairs=len(btw), between_above94=len(bv),
                            nonseparable=len(nonsep), nonseparable_strains=";".join(sorted(nonsep)),
                            nn_other_strain=nnbad, n_labelled=len(L),
                            cl_clusters=len(cl), ari_vs_published=round(a_ri, 4),
                            L_vs_R=round(float(lr), 2),
                            reps_below78=below78, rep_pairs=len(rp),
                            maxC_L=round(float(cl_max), 2)))

out = os.path.join(AN, "27_threshold_definitive.csv")
with open(out, "w", newline="", encoding="utf-8-sig") as fh:
    w = csv.DictWriter(fh, fieldnames=list(summary[0].keys()))
    w.writeheader(); w.writerows(summary)
print(f"\nwrote {out}")

# ---- how far apart are the MSA and pairwise matrices?
ma, mp = SOURCES["MSA"], SOURCES["pairwise"]
common = [a for a in mp[0] if a in set(ma[0])]
ia = {a: i for i, a in enumerate(ma[0])}; ip = {a: i for i, a in enumerate(mp[0])}
for conv in ("pairwise deletion", "gaps as differences"):
    A = ma[1][conv]; B = mp[1][conv]
    v1, v2 = [], []
    for i, x in enumerate(common):
        for y in common[i+1:]:
            v1.append(A[ia[x], ia[y]]); v2.append(B[ip[x], ip[y]])
    v1 = np.array(v1); v2 = np.array(v2); d = v2 - v1
    print(f"\nMSA vs true pairwise ({conv}): n={len(d)}  "
          f"mean diff {d.mean():+.3f} pp  mean|diff| {np.abs(d).mean():.3f} pp  "
          f"max|diff| {np.abs(d).max():.3f} pp  r={np.corrcoef(v1, v2)[0,1]:.5f}")
