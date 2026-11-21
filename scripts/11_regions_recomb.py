#!/usr/bin/env python
"""Step 11: region-wise diversity (LIR/SIR handled correctly), sliding-window
diversity, and a direct test of whether coat-protein-based strain assignment
agrees with whole-genome and Rep-based assignment.
"""
import csv, os, json, itertools
import numpy as np
from collections import Counter, defaultdict

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
strainf = {r["accession"]: r["strain_final"] for r in csv.DictReader(open(os.path.join(OUT, "06_strain_assignment.csv"), encoding="utf-8-sig"))}
cfg = json.load(open(os.path.join(OUT, "10_analysis_set.json")))
keep = cfg["keep"]

A = np.frombuffer(b"".join(aln[a].upper().encode() for a in keep), dtype=np.uint8).reshape(len(keep), len(aln[keep[0]]))
isb = A != ord("-")
for c in "RYSWKMBDHVN":
    isb &= (A != ord(c))
NCOL = A.shape[1]

# region boundaries in alignment columns, from step 10's mapping of NC_011058
# (rotated genome begins at the Rep nick inside the large intergenic region)
REGIONS = [
    ("LIR (large intergenic region)", [(0, 181), (2632, NCOL)]),
    ("V1 / movement protein",          [(181, 507)]),
    ("V2 / capsid protein",            [(519, 1264)]),
    ("SIR (small intergenic region)",  [(1264, 1437)]),
    ("C1 / RepA",                      [(1650, 2632)]),
    ("C1:C2 / Rep (spliced)",          [(1437, 2632)]),
    ("Rep exon 1",                     [(1437, 1859)]),
    ("Rep exon 2",                     [(1945, 2632)]),
]

def cols_of(blocks):
    c = []
    for lo, hi in blocks:
        c.extend(range(lo, hi))
    return np.array(c, dtype=int)

def stats(rows, cols, mincov=0.90):
    sub, sb = A[np.ix_(rows, cols)], isb[np.ix_(rows, cols)]
    n = len(rows)
    ok = sb.mean(0) >= mincov            # column present in >=90% of sequences
    if ok.sum() < 30 or n < 2:
        return None
    S_, B_ = sub[:, ok], sb[:, ok]
    L = int(ok.sum())
    seg = 0
    for c in range(L):
        v = set(S_[B_[:, c], c].tolist())
        if len(v) > 1:
            seg += 1
    tot = cnt = 0
    for i, j in itertools.combinations(range(n), 2):
        m = B_[i] & B_[j]
        if m.sum() < 30:
            continue
        tot += (S_[i][m] != S_[j][m]).sum() / m.sum(); cnt += 1
    pi = tot / cnt if cnt else 0.0
    haps = len({S_[i].tobytes() for i in range(n)})
    freq = Counter(S_[i].tobytes() for i in range(n))
    hd = (n / (n - 1)) * (1 - sum((v / n) ** 2 for v in freq.values()))
    a1 = sum(1 / i for i in range(1, n)); a2 = sum(1 / i**2 for i in range(1, n))
    b1 = (n + 1) / (3 * (n - 1)); b2 = 2 * (n**2 + n + 3) / (9 * n * (n - 1))
    c1 = b1 - 1 / a1; c2 = b2 - (n + 2) / (a1 * n) + a2 / a1**2
    e1 = c1 / a1; e2 = c2 / (a1**2 + a2)
    k = pi * L
    var = e1 * seg + e2 * seg * (seg - 1)
    D = (k - seg / a1) / np.sqrt(var) if var > 0 and n >= 4 else float("nan")
    return dict(n=n, sites=L, S=seg, pi=pi, k=k, haplotypes=haps,
                hap_diversity=hd, tajimaD=D)

rowsall = list(range(len(keep)))
print("=== diversity by genome region (columns present in >=90% of genomes) ===")
print(f"{'region':<34s} {'sites':>6} {'S':>6} {'S/site':>7} {'pi':>8} {'Hd':>7} {'D':>7}")
rows_out = []
for name, blocks in REGIONS:
    st = stats(rowsall, cols_of(blocks))
    if st is None:
        print(f"{name:<34s} (insufficient)"); continue
    print(f"{name:<34s} {st['sites']:>6} {st['S']:>6} {st['S']/st['sites']:>7.3f} "
          f"{st['pi']:>8.4f} {st['hap_diversity']:>7.4f} {st['tajimaD']:>7.3f}")
    rows_out.append(dict(region=name, **st))
with open(os.path.join(OUT, "11_region_diversity.csv"), "w", newline="", encoding="utf-8-sig") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows_out[0].keys())); w.writeheader(); w.writerows(rows_out)

# ---------------------------------------------------- sliding window
print("\ncomputing sliding-window diversity (200 nt window, 25 nt step)...")
win, step = 200, 25
wr = []
for lo in range(0, NCOL - win + 1, step):
    st = stats(rowsall, np.arange(lo, lo + win), mincov=0.80)
    if st:
        wr.append(dict(window_start=lo, window_mid=lo + win // 2,
                       sites=st["sites"], S=st["S"], pi=round(st["pi"], 5)))
with open(os.path.join(OUT, "11_sliding_window_pi.csv"), "w", newline="", encoding="utf-8-sig") as fh:
    w = csv.DictWriter(fh, fieldnames=list(wr[0].keys())); w.writeheader(); w.writerows(wr)
top = sorted(wr, key=lambda r: -r["pi"])[:6]
print("  most variable windows (alignment coords):",
      [(r["window_mid"], r["pi"]) for r in top])
print(f"  wrote 11_sliding_window_pi.csv ({len(wr)} windows)")

# ------------------------------- strain assignment: CP vs Rep vs whole genome
print("\n=== does coat-protein-based strain assignment agree with the genome? ===")
ann = {a: meta[a]["strain"] for a in keep
       if len(meta[a]["strain"]) == 1 and meta[a]["strain"].isalpha()}
PARTS = {
    "whole genome": cols_of([(0, NCOL)]),
    "V2 / capsid protein": cols_of([(519, 1264)]),
    "C1:C2 / Rep": cols_of([(1437, 2632)]),
}
ri = {a: i for i, a in enumerate(keep)}

def assign(cols):
    sub, sb = A[:, cols], isb[:, cols]
    out = {}
    for a in keep:
        i = ri[a]
        best = (-1, None)
        for b in ann:
            j = ri[b]
            if b == a:
                continue
            m = sb[i] & sb[j]
            if m.sum() < 100:
                continue
            v = (sub[i][m] == sub[j][m]).mean()
            if v > best[0]:
                best = (v, b)
        out[a] = (ann[best[1]] if best[1] else "?", best[0], best[1])
    return out

res = {k: assign(v) for k, v in PARTS.items()}
disc = []
for a in keep:
    wg, cp, rp = res["whole genome"][a], res["V2 / capsid protein"][a], res["C1:C2 / Rep"][a]
    if cp[0] != rp[0] or cp[0] != wg[0]:
        disc.append(dict(accession=a, country=meta[a]["country"],
                         host=meta[a]["host_common"],
                         annotated_strain=meta[a]["strain"],
                         nearest_by_whole_genome=f"{wg[0]} ({wg[2]}, {wg[1]*100:.1f}%)",
                         nearest_by_capsid=f"{cp[0]} ({cp[2]}, {cp[1]*100:.1f}%)",
                         nearest_by_Rep=f"{rp[0]} ({rp[2]}, {rp[1]*100:.1f}%)",
                         cp_vs_rep_conflict=cp[0] != rp[0],
                         cp_vs_genome_conflict=cp[0] != wg[0]))
n_cr = sum(1 for d in disc if d["cp_vs_rep_conflict"])
n_cg = sum(1 for d in disc if d["cp_vs_genome_conflict"])
print(f"  genomes analysed: {len(keep)}")
print(f"  nearest-strain call differs between capsid and Rep: {n_cr} "
      f"({100*n_cr/len(keep):.1f}%)")
print(f"  capsid-based call differs from whole-genome call:    {n_cg} "
      f"({100*n_cg/len(keep):.1f}%)")
print("\n  examples:")
for d in disc[:14]:
    print(f"    {d['accession']} {d['country']:<12s} ann={d['annotated_strain'] or '-':<3s} "
          f"genome->{d['nearest_by_whole_genome']:<22s} CP->{d['nearest_by_capsid']:<22s} "
          f"Rep->{d['nearest_by_Rep']}")
if disc:
    with open(os.path.join(OUT, "11_strain_assignment_discordance.csv"), "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(disc[0].keys())); w.writeheader(); w.writerows(disc)
    print(f"\n  wrote 11_strain_assignment_discordance.csv ({len(disc)} rows)")

# summary of how often CP and Rep nearest-neighbour strains disagree, by strain
cm = Counter((res["V2 / capsid protein"][a][0], res["C1:C2 / Rep"][a][0]) for a in keep)
print("\n  capsid-strain x Rep-strain combinations seen more than once "
      "(off-diagonal = mosaic signal):")
for (cp, rp), n in sorted(cm.items(), key=lambda kv: -kv[1]):
    if cp != rp and n >= 1:
        print(f"    CP={cp}  Rep={rp}  n={n}")
