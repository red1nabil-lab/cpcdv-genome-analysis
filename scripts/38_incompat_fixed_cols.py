#!/usr/bin/env python
"""Step 38: redo PHI and Hudson-Kaplan on IDENTICAL alignment columns for every group.

Reviewer round 3, point 1.  The previous run reported Rm >= 197 for all 267 genomes but
Rm >= 211 for a 131-genome subset of them.  Under the four-gamete method, removing
sequences cannot create incompatibilities, so a subset bound cannot exceed the full-set
bound on the same columns -- the result signalled a bug, and it was one.

Cause: the helper recomputed "biallelic and complete-case" columns WITHIN each group, so
a column gapped in one of the excluded sequences became usable once that sequence was
dropped.  The subsets were therefore scored on 452 and 450 columns against 441 for the
full set, and were not comparable.

Fix: the column set is defined ONCE, from the full analysis set, and every group is
scored on exactly those homologous columns.  A column that becomes monomorphic within a
subset simply contributes no incompatibility.  Rm for any subset must now be <= Rm for
the full set, and that is asserted.

Assumptions documented: the four-gamete test assumes an infinite-sites model, so Rm is a
lower bound on recombination events and recurrent mutation inflates it.
"""
import os, sys, json, csv
from collections import defaultdict
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
import os as _os, sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # noqa: E402
from cpcdv_paths import DATA_DIR, FIG_DIR, MUSCLE_BIN  # noqa: E402
AN = DATA_DIR
rng = np.random.default_rng(20261123)
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
ix = {a: i for i, a in enumerate(keep)}
sass = {r["accession"]: r["strain_final"] for r in
        csv.DictReader(open(os.path.join(AN, "06_strain_assignment.csv"),
                            encoding="utf-8-sig"))}
by = defaultdict(list)
for a in keep:
    by[sass.get(a, "unassigned")].append(a)
rec = list(csv.DictReader(open(os.path.join(AN, "29_recombination_maxchi.csv"),
                               encoding="utf-8-sig")))
SIG = {r["accession"] for r in rec if r["significant"] == "True"}

# ---------------------------------------------------------------- fixed columns
usable = np.ones(A.shape[1], dtype=bool)
for c in AMB | {GAP}:
    usable &= ~(A == c).any(axis=0)          # complete-case across ALL 267
COLS = []
for c in np.flatnonzero(usable):
    u = np.unique(A[:, c])
    if u.size == 2:
        mn = (A[:, c] == u[1]).sum()
        if 1 < mn < len(keep) - 1:           # parsimony informative in the full set
            COLS.append(int(c))
COLS = np.array(COLS)
print(f"alignment columns: {A.shape[1]}")
print(f"complete-case across all {len(keep)} genomes: {int(usable.sum())}")
print(f"biallelic + parsimony-informative (THE FIXED SET): {len(COLS)}")
PAT_ALL = A[:, COLS]


def patterns(rows):
    """0/1 coding of the fixed columns for the given rows, using the full-set
    major/minor allele definition so coding is identical across groups"""
    sub = PAT_ALL[rows]
    out = np.zeros(sub.shape, dtype=np.int8)
    for j in range(sub.shape[1]):
        col_all = PAT_ALL[:, j]
        u = np.unique(col_all)
        out[:, j] = (sub[:, j] == u[1]).astype(np.int8)
    return out


def incompat(P):
    """four-gamete test between all column pairs; columns monomorphic in this group
    cannot be incompatible with anything"""
    n = P.shape[1]
    I = np.zeros((n, n), dtype=bool)
    poly = np.array([len(np.unique(P[:, j])) > 1 for j in range(n)])
    for i in range(n):
        if not poly[i]:
            continue
        a = P[:, i]
        g00 = ((a[:, None] == 0) & (P == 0)).any(axis=0)
        g01 = ((a[:, None] == 0) & (P == 1)).any(axis=0)
        g10 = ((a[:, None] == 1) & (P == 0)).any(axis=0)
        g11 = ((a[:, None] == 1) & (P == 1)).any(axis=0)
        I[i] = g00 & g01 & g10 & g11 & poly
    np.fill_diagonal(I, False)
    return I


def phi_stat(I, order, w=100):
    n = len(order)
    tot = 0.0; cnt = 0
    for d in range(1, min(w, n - 1) + 1):
        a = order[:-d]; b = order[d:]
        tot += I[a, b].sum(); cnt += len(a)
    return tot / cnt if cnt else float("nan")


def rm_bound(I, cols):
    """Hudson-Kaplan lower bound: max set of disjoint intervals each spanning an
    incompatible pair (greedy by right endpoint)."""
    n = len(cols)
    pairs = []
    for i in range(n):
        js = np.flatnonzero(I[i])
        for j in js[js > i]:
            pairs.append((cols[i], cols[j]))
    if not pairs:
        return 0
    pairs.sort(key=lambda p: p[1])
    rm = 0; last = -1
    for s, e in pairs:
        if s >= last:
            rm += 1; last = e
    return rm


def run(name, members, nperm=500, w=100):
    rows = [ix[a] for a in members]
    P = patterns(rows)
    I = incompat(P)
    obs = phi_stat(I, np.arange(len(COLS)), w)
    null = np.array([phi_stat(I, rng.permutation(len(COLS)), w) for _ in range(nperm)])
    p = (1.0 + int((null <= obs).sum())) / (nperm + 1.0)
    rm = rm_bound(I, COLS)
    npoly = int(sum(len(np.unique(P[:, j])) > 1 for j in range(P.shape[1])))
    print(f"  {name:28s} n={len(members):4d}  polymorphic {npoly:4d}/{len(COLS)}  "
          f"PHI={obs:.4f} (null {null.mean():.4f})  p={p:.4f}  Rm>={rm}")
    return dict(group=name, n=len(members), fixed_columns=int(len(COLS)),
                polymorphic_in_group=npoly, phi=float(obs),
                phi_null_mean=float(null.mean()), phi_p=float(p),
                hudson_kaplan_rm=int(rm))


print("\n=== PHI and Hudson-Kaplan on one fixed homologous column set ===")
res = [run("all genomes", keep, nperm=200)]
RM_FULL = res[0]["hudson_kaplan_rm"]
for s in ("H", "C", "A", "F", "D", "K"):
    if len(by.get(s, [])) >= 6:
        res.append(run(f"strain {s}", by[s], nperm=500))
clean = [a for a in keep if a not in SIG]
sig = [a for a in keep if a in SIG]
res.append(run("no event localised", clean, nperm=200))
res.append(run("event localised", sig, nperm=200))

print("\n--- monotonicity check (subset Rm must not exceed full-set Rm) ---")
bad = [r for r in res if r["group"] != "all genomes"
       and r["hudson_kaplan_rm"] > RM_FULL]
print(f"  full set Rm >= {RM_FULL}")
for r in res[1:]:
    flag = "OK " if r["hudson_kaplan_rm"] <= RM_FULL else "VIOLATION"
    print(f"  {flag} {r['group']:28s} Rm >= {r['hudson_kaplan_rm']}")
assert not bad, "subset exceeds full-set bound: implementation still wrong"
print("  -> all subsets respect the bound; the earlier anomaly was the column-set bug")

out = os.path.join(AN, "38_incompat_fixed_cols.csv")
with open(out, "w", newline="", encoding="utf-8-sig") as fh:
    w_ = csv.DictWriter(fh, fieldnames=list(res[0].keys()))
    w_.writeheader(); w_.writerows(res)
print(f"\nwrote {out}")
