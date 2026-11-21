#!/usr/bin/env python
"""Step 6: pairwise identity, species/strain threshold checks, strain assignment.

Species and strain limits in Mastrevirus are 78% and 94% genome-wide pairwise
nucleotide identity (Muhire et al. 2013).  Identity is reported two ways:
  pid_pw  - pairwise deletion: matches / columns where both sequences have a base
  pid_gap - gaps counted as differences: matches / columns where either has a base
pid_gap is the more conservative measure and is used for the threshold tests.
"""
import csv, os, json
import numpy as np

import os as _os, sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # noqa: E402
from cpcdv_paths import DATA_DIR, FIG_DIR, MUSCLE_BIN  # noqa: E402
OUT = DATA_DIR

def read_fasta(p):
    d, k, buf = {}, None, []
    for line in open(p, encoding="utf8"):
        line = line.strip()
        if line.startswith(">"):
            if k: d[k] = "".join(buf)
            k, buf = line[1:].split()[0], []
        elif line:
            buf.append(line)
    if k: d[k] = "".join(buf)
    return d

aln = read_fasta(os.path.join(OUT, "05_full_aln.fasta"))
meta = {r["accession"]: r for r in csv.DictReader(open(os.path.join(OUT, "03_curated_metadata.csv"), encoding="utf-8-sig"))}
accs = sorted(aln)
N, L = len(accs), len(aln[accs[0]])
print(f"alignment: {N} sequences x {L} columns")

A = np.frombuffer(b"".join(aln[a].upper().encode() for a in accs), dtype=np.uint8).reshape(N, L)
GAP = ord("-")
AMB = {ord(c) for c in "RYSWKMBDHVN"}
isbase = (A != GAP)
for c in AMB:
    isbase &= (A != c)                      # ambiguity treated as missing

pid_pw = np.zeros((N, N)); pid_gap = np.zeros((N, N))
cmp_pw = np.zeros((N, N), dtype=int)
for i in range(N):
    same = (A == A[i])
    both = isbase & isbase[i]
    either = isbase | isbase[i]
    m = (same & both).sum(1)
    cmp_pw[i] = both.sum(1)
    pid_pw[i] = m / np.maximum(both.sum(1), 1)
    pid_gap[i] = m / np.maximum(either.sum(1), 1)
np.fill_diagonal(pid_pw, 1.0); np.fill_diagonal(pid_gap, 1.0)

off = ~np.eye(N, dtype=bool)
print(f"\npid_pw  : min={pid_pw[off].min()*100:.2f}%  mean={pid_pw[off].mean()*100:.2f}%  max={pid_pw[off].max()*100:.2f}%")
print(f"pid_gap : min={pid_gap[off].min()*100:.2f}%  mean={pid_gap[off].mean()*100:.2f}%  max={pid_gap[off].max()*100:.2f}%")
print(f"mean columns compared per pair: {cmp_pw[off].mean():.0f}")

# ---------------------------------------------------- species threshold (78%)
idx = {a: i for i, a in enumerate(accs)}
print("\n=== genomes falling below the 78% species threshold against the bulk ===")
below = []
for a in accs:
    i = idx[a]
    row = pid_gap[i][off[i]]
    frac = (row < 0.78).mean()
    if frac > 0.5:
        below.append((a, row.max(), frac))
if below:
    for a, mx, frac in sorted(below, key=lambda t: t[1]):
        m = meta[a]
        print(f"  {a}: max identity to any other genome {mx*100:.2f}%; "
              f"below 78% against {frac*100:.0f}% of genomes")
        print(f"      iso={m['isolate']} country={m['country']} org={m['gb_organism']}")
        print(f"      ref: {m['reference_title'][:100]}")
else:
    print("  none")

# per-genome nearest neighbour, to spot any other outliers
nn = [(a, max(pid_gap[idx[a]][j] for j in range(N) if j != idx[a])) for a in accs]
nn.sort(key=lambda t: t[1])
print("\nten genomes most isolated from the rest (highest identity to any other genome):")
for a, v in nn[:10]:
    m = meta[a]
    print(f"  {a} {v*100:6.2f}%  {m['country']:<12s} {m['host_common']:<22s} "
          f"strain={m['strain'] or '-':<4s} {m['isolate'][:26]}")

# ---------------------------------------------------- strain consistency
ann = {a: meta[a]["strain"] for a in accs
       if len(meta[a]["strain"]) == 1 and meta[a]["strain"].isalpha()}
from collections import defaultdict
bys = defaultdict(list)
for a, s in ann.items():
    bys[s].append(a)
print(f"\n=== annotated strains among the {N} full genomes ===")
print({k: len(v) for k, v in sorted(bys.items())})

print("\nwithin-strain identity range (pid_gap) and worst within-strain pair:")
rows_ws = []
for s, members in sorted(bys.items()):
    if len(members) < 2:
        print(f"  {s}: n=1, not assessable")
        rows_ws.append(dict(strain=s, n=1, min_within="", mean_within="", breaches=""))
        continue
    ii = [idx[a] for a in members]
    sub = pid_gap[np.ix_(ii, ii)]
    m = ~np.eye(len(ii), dtype=bool)
    mn, mean = sub[m].min(), sub[m].mean()
    nb = int((sub[m] < 0.94).sum() / 2)
    print(f"  {s}: n={len(members):>3} min={mn*100:6.2f}% mean={mean*100:6.2f}% "
          f"pairs below 94%: {nb}")
    rows_ws.append(dict(strain=s, n=len(members), min_within=round(mn*100, 2),
                        mean_within=round(mean*100, 2), breaches=nb))

print("\nbetween-strain maxima that exceed 94% (annotated strains that merge under the rule):")
merge = []
ks = sorted(bys)
for x in range(len(ks)):
    for y in range(x+1, len(ks)):
        ii = [idx[a] for a in bys[ks[x]]]; jj = [idx[a] for a in bys[ks[y]]]
        mx = pid_gap[np.ix_(ii, jj)].max()
        if mx > 0.94:
            merge.append((ks[x], ks[y], mx))
for x, y, mx in sorted(merge, key=lambda t: -t[2]):
    print(f"  {x} vs {y}: max {mx*100:.2f}%")
if not merge:
    print("  none")

# ---------------------------------------------------- assign unannotated
print("\n=== strain assignment of genomes without an annotated strain ===")
assign = []
for a in accs:
    if a in ann:
        assign.append(dict(accession=a, strain_final=ann[a], basis="annotation",
                           nearest="", nearest_identity=""))
        continue
    i = idx[a]
    best = max(((pid_gap[i][idx[b]], b) for b in ann), key=lambda t: t[0])
    v, b = best
    if v >= 0.94:
        assign.append(dict(accession=a, strain_final=ann[b],
                           basis=f">=94% to a strain-{ann[b]} genome",
                           nearest=b, nearest_identity=round(v*100, 2)))
    else:
        assign.append(dict(accession=a, strain_final="unassigned",
                           basis="<94% to every annotated strain",
                           nearest=b, nearest_identity=round(v*100, 2)))
from collections import Counter
print(Counter(r["strain_final"] for r in assign).most_common())
print("\ngenomes that remain unassigned (candidate unclassified lineages):")
for r in assign:
    if r["strain_final"] == "unassigned":
        m = meta[r["accession"]]
        print(f"  {r['accession']} nearest {r['nearest']} ({m and ann[r['nearest']]}) "
              f"{r['nearest_identity']}%  {m['country']:<12s} {m['host_common']:<20s} {m['isolate'][:28]}")

# ---------------------------------------------------- outputs
np.save(os.path.join(OUT, "06_pid_gap.npy"), pid_gap)
np.save(os.path.join(OUT, "06_pid_pw.npy"), pid_pw)
with open(os.path.join(OUT, "06_accession_order.json"), "w") as fh:
    json.dump(accs, fh)
with open(os.path.join(OUT, "06_strain_assignment.csv"), "w", newline="", encoding="utf-8-sig") as fh:
    w = csv.DictWriter(fh, fieldnames=list(assign[0].keys())); w.writeheader(); w.writerows(assign)
with open(os.path.join(OUT, "06_within_strain_identity.csv"), "w", newline="", encoding="utf-8-sig") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows_ws[0].keys())); w.writeheader(); w.writerows(rows_ws)
# full matrix as csv (percent, 2 dp)
with open(os.path.join(OUT, "06_pairwise_identity_matrix.csv"), "w", newline="", encoding="utf-8-sig") as fh:
    w = csv.writer(fh); w.writerow([""] + accs)
    for i, a in enumerate(accs):
        w.writerow([a] + [f"{pid_gap[i, j]*100:.2f}" for j in range(N)])
print("\nwrote 06_* outputs")
