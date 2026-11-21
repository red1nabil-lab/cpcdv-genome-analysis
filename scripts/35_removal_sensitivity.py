#!/usr/bin/env python
"""Step 35: is removing the recombinants unusually effective, or just removing genomes?

Reviewer point 2.  Dropping the 136 genomes with recombination signal eliminates all 28
within-strain threshold violations, but dropping any 136 genomes removes comparisons
too, so the association could be an artefact of sample size.  The claim needs a null.

Two nulls, 10,000 replicates each:
  (1) unstratified - remove the same NUMBER of labelled genomes at random
  (2) strain-stratified - remove the same number FROM EACH STRAIN as the recombinant
      set removes, so strain composition is held exactly constant; this is the stronger
      control, because the violations are concentrated in one strain and an unstratified
      draw could remove strain C by chance

Also reports, for each of the 28 violating pairs, whether one or both members carry
recombination signal - which is what licenses the wording used in the manuscript.
"""
import os, sys, json, csv
from collections import defaultdict, Counter
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
import os as _os, sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # noqa: E402
from cpcdv_paths import DATA_DIR, FIG_DIR, MUSCLE_BIN  # noqa: E402
AN = DATA_DIR
rng = np.random.default_rng(20261121)
NPERM = 10000

cfg = json.load(open(os.path.join(AN, "10_analysis_set.json")))
keep = list(cfg["keep"])
M = np.load(os.path.join(AN, "25_truepw_pid_pw_A.npy")).astype(float)
if np.nanmax(M) <= 1.5:
    M *= 100.0
ix = {a: i for i, a in enumerate(keep)}
meta = {r["accession"]: r for r in
        csv.DictReader(open(os.path.join(AN, "03_curated_metadata.csv"),
                            encoding="utf-8-sig"))}
rec = list(csv.DictReader(open(os.path.join(AN, "29_recombination_maxchi.csv"),
                               encoding="utf-8-sig")))
SIG = {r["accession"] for r in rec if r["significant"] == "True"}

pub = {a: meta[a]["strain"] for a in keep
       if a in meta and len(meta[a]["strain"]) == 1 and meta[a]["strain"].isalpha()}
L = sorted(pub)
by = defaultdict(list)
for a in L:
    by[pub[a]].append(a)

# observed
def violations(subset):
    s = set(subset)
    n = 0
    pairs = []
    for st, mem in by.items():
        m = [a for a in mem if a in s]
        for i, a in enumerate(m):
            for b in m[i + 1:]:
                if M[ix[a], ix[b]] < 94.0:
                    n += 1; pairs.append((a, b, st, M[ix[a], ix[b]]))
    return n, pairs


obs_all, obs_pairs = violations(L)
clean_lab = [a for a in L if a not in SIG]
obs_clean, _ = violations(clean_lab)
n_removed = len(L) - len(clean_lab)
print(f"labelled genomes {len(L)}; with recombination signal {n_removed}; "
      f"clean {len(clean_lab)}")
print(f"within-strain violations: all {obs_all}, recombination-free {obs_clean}")

# ---------------------------------------------------------------- the 28 pairs
print(f"\n=== composition of the {obs_all} violating pairs ===")
both = one = neither = 0
for a, b, st, v in obs_pairs:
    k = (a in SIG) + (b in SIG)
    both += k == 2; one += k == 1; neither += k == 0
print(f"  both members carry recombination signal : {both}")
print(f"  exactly one member                      : {one}")
print(f"  neither member                          : {neither}")
print(f"  -> {both+one} of {obs_all} involve at least one genome with supported signal")
print(f"  strains involved: {dict(Counter(st for _,_,st,_ in obs_pairs))}")

# ---------------------------------------------------------------- null 1
cnt_removed_per_strain = Counter(pub[a] for a in L if a in SIG)
null1 = np.empty(NPERM, dtype=int)
Larr = np.array(L)
for r in range(NPERM):
    keepmask = rng.permutation(len(L))[: len(L) - n_removed]
    null1[r], _ = violations(Larr[keepmask])
p1 = (1.0 + int((null1 <= obs_clean).sum())) / (NPERM + 1.0)
print(f"\n=== null 1: remove {n_removed} labelled genomes at random ===")
print(f"  violations remaining: mean {null1.mean():.2f}, median {np.median(null1):.0f}, "
      f"min {null1.min()}, max {null1.max()}")
print(f"  replicates reaching 0 violations: {int((null1 <= 0).sum())} / {NPERM} "
      f"({100*int((null1<=0).sum())/NPERM:.2f}%)")
print(f"  P(random removal does as well as removing recombinants) = {p1:.4f}")

# ---------------------------------------------------------------- null 2
null2 = np.empty(NPERM, dtype=int)
for r in range(NPERM):
    sel = []
    for st, mem in by.items():
        k = cnt_removed_per_strain.get(st, 0)
        idxs = rng.permutation(len(mem))[: len(mem) - k]
        sel.extend(mem[i] for i in idxs)
    null2[r], _ = violations(sel)
p2 = (1.0 + int((null2 <= obs_clean).sum())) / (NPERM + 1.0)
print(f"\n=== null 2: strain-stratified removal (same number per strain) ===")
print(f"  removals per strain: {dict(cnt_removed_per_strain)}")
print(f"  violations remaining: mean {null2.mean():.2f}, median {np.median(null2):.0f}, "
      f"min {null2.min()}, max {null2.max()}")
print(f"  replicates reaching 0 violations: {int((null2 <= 0).sum())} / {NPERM} "
      f"({100*int((null2<=0).sum())/NPERM:.2f}%)")
print(f"  P(stratified random removal does as well) = {p2:.4f}")

json.dump({"observed_all": obs_all, "observed_clean": obs_clean,
           "n_removed": n_removed, "n_labelled": len(L),
           "pairs_both_recombinant": both, "pairs_one_recombinant": one,
           "pairs_neither": neither,
           "null1_mean": float(null1.mean()), "null1_zero_frac": float((null1 <= 0).mean()),
           "null1_p": p1,
           "null2_mean": float(null2.mean()), "null2_zero_frac": float((null2 <= 0).mean()),
           "null2_p": p2, "n_perm": NPERM},
          open(os.path.join(AN, "35_removal_sensitivity.json"), "w"), indent=1)
print("\nwrote 35_removal_sensitivity.json")
