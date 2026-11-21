#!/usr/bin/env python
"""Step 30: how much do the diversity statistics depend on the gap filter?

Reviewer criticism 7: columns gapped in more than 10% of sequences were excluded,
reducing 2829 alignment columns to 2201, and the effect of that choice on the diversity
estimates was never evaluated.  Tajima's D in particular is sensitive to how missing
data are handled.

Recomputes pi, S, haplotype diversity and Tajima's D for the whole set and for the
largest strains across gap-filter thresholds from 0% (only fully ungapped columns) to
100% (no filtering at all).  If the headline values move materially the manuscript must
say so.
"""
import os, sys, json, csv
import numpy as np
from collections import defaultdict

sys.stdout.reconfigure(encoding="utf-8")
import os as _os, sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # noqa: E402
from cpcdv_paths import DATA_DIR, FIG_DIR, MUSCLE_BIN  # noqa: E402
AN = DATA_DIR
GAP = ord("-")
AMB = {ord(c) for c in "RYSWKMBDHVN"}


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
KEEP = [a for a in json.load(open(os.path.join(AN, "10_analysis_set.json")))["keep"]
        if a in aln]
A = np.frombuffer(b"".join(aln[a].upper().encode() for a in KEEP),
                  dtype=np.uint8).reshape(len(KEEP), -1)
sass = {r["accession"]: r["strain_final"] for r in
        csv.DictReader(open(os.path.join(AN, "06_strain_assignment.csv"),
                            encoding="utf-8-sig"))}
by = defaultdict(list)
for a in KEEP:
    by[sass.get(a, "unassigned")].append(a)
idx = {a: i for i, a in enumerate(KEEP)}

isbase = (A != GAP)
for c in AMB:
    isbase &= (A != c)
gapfrac = 1.0 - isbase.mean(axis=0)


def stats(rows, cols):
    """pi, S, haplotypes, Hd, Tajima's D over the given rows/columns."""
    S_ = A[np.ix_(rows, cols)]
    base = isbase[np.ix_(rows, cols)]
    n = len(rows)
    if n < 4:
        return None
    # pairwise differences with pairwise deletion
    tot = 0.0; cnt = 0; npair = 0
    for i in range(n):
        for j in range(i + 1, n):
            both = base[i] & base[j]
            nb = int(both.sum())
            if nb == 0:
                continue
            d = int(np.count_nonzero(S_[i][both] != S_[j][both]))
            tot += d / nb; cnt += d; npair += 1
    pi = tot / npair if npair else float("nan")
    k = cnt / npair if npair else float("nan")
    # segregating sites: >1 distinct base among sequences that have a base
    seg = 0
    for c in range(S_.shape[1]):
        col = S_[:, c][base[:, c]]
        if col.size > 1 and np.unique(col).size > 1:
            seg += 1
    # haplotypes
    haps = {S_[i].tobytes() for i in range(n)}
    hc = defaultdict(int)
    for i in range(n):
        hc[S_[i].tobytes()] += 1
    hd = (n / (n - 1.0)) * (1.0 - sum((v / n) ** 2 for v in hc.values()))
    # Tajima's D
    a1 = sum(1.0 / i for i in range(1, n))
    a2 = sum(1.0 / i ** 2 for i in range(1, n))
    b1 = (n + 1.0) / (3.0 * (n - 1)); b2 = 2.0 * (n * n + n + 3) / (9.0 * n * (n - 1))
    c1 = b1 - 1.0 / a1; c2 = b2 - (n + 2.0) / (a1 * n) + a2 / (a1 ** 2)
    e1 = c1 / a1; e2 = c2 / (a1 ** 2 + a2)
    var = e1 * seg + e2 * seg * (seg - 1)
    D = (k - seg / a1) / np.sqrt(var) if var > 0 else float("nan")
    return dict(n=n, sites=len(cols), S=seg, pi=pi, k=k,
                haplotypes=len(haps), hd=hd, tajimaD=D)


THRESH = [0.0, 0.02, 0.05, 0.10, 0.20, 0.50, 1.01]
groups = [("all genomes", KEEP)] + [(f"strain {s}", by[s]) for s in ("H", "C", "A", "F")]

rows_out = []
print(f"alignment columns: {A.shape[1]}")
for gname, members in groups:
    rr = [idx[a] for a in members]
    print(f"\n=== {gname} (n={len(members)}) ===")
    print(f"  {'gap filter':>12s} {'cols':>6s} {'S':>5s} {'pi':>8s} {'k':>8s} "
          f"{'Hd':>7s} {'Tajima D':>9s}")
    for t in THRESH:
        cols = np.flatnonzero(gapfrac <= t)
        if cols.size == 0:
            continue
        st = stats(rr, cols)
        if st is None:
            continue
        lab = "none" if t > 1 else f"<={t:.0%}"
        print(f"  {lab:>12s} {st['sites']:6d} {st['S']:5d} {st['pi']:8.4f} "
              f"{st['k']:8.2f} {st['hd']:7.4f} {st['tajimaD']:+9.3f}")
        rows_out.append(dict(group=gname, gap_filter=lab, **{k: (round(v, 5)
                        if isinstance(v, float) else v) for k, v in st.items()}))

out = os.path.join(AN, "30_gap_sensitivity.csv")
with open(out, "w", newline="", encoding="utf-8-sig") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows_out[0].keys()))
    w.writeheader(); w.writerows(rows_out)
print(f"\nwrote {out}")

# summarise the spread for the headline numbers
allr = [r for r in rows_out if r["group"] == "all genomes"]
pis = [r["pi"] for r in allr]; ds = [r["tajimaD"] for r in allr]
print(f"\nacross all gap filters, all genomes: pi {min(pis):.4f}-{max(pis):.4f} "
      f"(spread {max(pis)-min(pis):.4f}), Tajima D {min(ds):+.3f} to {max(ds):+.3f}")
for s in ("H", "C", "A", "F"):
    g = [r for r in rows_out if r["group"] == f"strain {s}"]
    if g:
        d = [r["tajimaD"] for r in g]
        print(f"  strain {s}: Tajima D range {min(d):+.3f} to {max(d):+.3f} "
              f"-> {'always negative' if max(d) < 0 else 'SIGN CHANGES'}")
