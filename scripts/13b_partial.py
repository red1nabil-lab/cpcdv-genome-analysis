#!/usr/bin/env python
"""Step 13b: the 77 partial records, assessed by local pairwise alignment.

A short fragment placed in a 96-sequence progressive alignment can be scattered
across the whole genome, which produces meaningless identities (one 388 nt
fragment appeared to span 1750 alignment columns at 56% identity).  Each
fragment is therefore aligned locally against each strain reference genome, and
both its genome position and its identity are read from that alignment.
"""
import csv, os, json, itertools
import numpy as np
from collections import Counter
from Bio import Align

import os as _os, sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # noqa: E402
from cpcdv_paths import DATA_DIR, FIG_DIR, MUSCLE_BIN  # noqa: E402
OUT = DATA_DIR
RC = str.maketrans("ACGTRYSWKMBDHVNacgtryswkmbdhvn", "TGCAYRSWMKVHDBNtgcayrswmkvhdbn")
def rc(s): return s.translate(RC)[::-1]

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

allseq = read_fasta(os.path.join(OUT, "01_all_dedup.fasta"))
ori = read_fasta(os.path.join(OUT, "04_full_oriented.fasta"))
meta = {r["accession"]: r for r in csv.DictReader(open(os.path.join(OUT, "03_curated_metadata.csv"), encoding="utf-8-sig"))}

part = [a for a in allseq if meta[a]["molecule_class"] != "full genome"]
reps = {}
for a in sorted(ori):
    s = meta[a]["strain"]
    if len(s) == 1 and s.isalpha() and s not in reps:
        reps[s] = a
print(f"{len(part)} fragments vs {len(reps)} strain references: {sorted(reps)}")

al = Align.PairwiseAligner()
al.mode = "local"
al.match_score = 2; al.mismatch_score = -3
al.open_gap_score = -5; al.extend_gap_score = -2

# feature spans in the rotated coordinate frame (from NC_011058, step 10)
GENES = [("V1 / movement protein", 131, 410), ("V2 / capsid protein", 422, 1160),
         ("SIR", 1160, 1317), ("C1:C2 / Rep + RepA", 1317, 2408)]
REFLEN = len(ori["NC_011058"])

def describe(lo, hi):
    """which annotated features does reference span [lo,hi) overlap?"""
    out = []
    for name, gs, ge in GENES:
        ov = max(0, min(hi, ge) - max(lo, gs))
        if ov >= 50:
            out.append(f"{name} ({ov} nt)")
    return "; ".join(out) or "intergenic / unassigned"

rows = []
for a in part:
    s0 = allseq[a].upper()
    best = None
    for st, s in (("+", s0), ("-", rc(s0))):
        for lab, racc in reps.items():
            ref = ori[racc]
            aln = al.align(ref + ref, s)[0]
            # coordinates and identity of the aligned block
            tgt, qry = str(aln[0]), str(aln[1])
            m = sum(1 for x, y in zip(tgt, qry) if x == y and x != "-")
            cols = sum(1 for x, y in zip(tgt, qry) if x != "-" or y != "-")
            ident = m / cols if cols else 0
            rs = int(aln.coordinates[0][0]) % REFLEN
            re_ = rs + (int(aln.coordinates[0][-1]) - int(aln.coordinates[0][0]))
            cand = (ident, aln.score, st, lab, racc, rs, re_, cols,
                    int(aln.coordinates[1][0]), int(aln.coordinates[1][-1]))
            if best is None or cand[1] > best[1]:
                best = cand
    ident, score, st, lab, racc, rs, re_, cols, qs, qe = best
    # identity to every reference, in the chosen orientation
    sq = s0 if st == "+" else rc(s0)
    allid = {}
    for lab2, racc2 in reps.items():
        aln = al.align(ori[racc2] + ori[racc2], sq)[0]
        t, q = str(aln[0]), str(aln[1])
        mm = sum(1 for x, y in zip(t, q) if x == y and x != "-")
        cc = sum(1 for x, y in zip(t, q) if x != "-" or y != "-")
        allid[lab2] = mm / cc if cc else 0
    srt = sorted(allid.items(), key=lambda kv: -kv[1])
    rows.append(dict(accession=a, length=len(s0), country=meta[a]["country"],
                     host=meta[a]["host_common"], host_family=meta[a]["host_family"],
                     year=meta[a]["collection_year"],
                     orientation_vs_reference=st,
                     ref_start=rs, ref_end=re_, aligned_nt=cols,
                     fragment_covers=describe(rs, re_),
                     nearest_strain=srt[0][0],
                     identity_to_nearest=round(srt[0][1]*100, 2),
                     second_strain=srt[1][0],
                     identity_to_second=round(srt[1][1]*100, 2),
                     margin=round((srt[0][1]-srt[1][1])*100, 2),
                     assignment_confident="yes" if (srt[0][1] >= 0.94 and (srt[0][1]-srt[1][1]) >= 0.01) else "no"))

print(f"\nreverse-complemented relative to the reference frame: "
      f"{sum(1 for r in rows if r['orientation_vs_reference']=='-')}")
print("\ngenome region covered:")
for k, v in Counter(r["fragment_covers"] for r in rows).most_common():
    print(f"  {k}: {v}")
ids = [r["identity_to_nearest"] for r in rows]
print(f"\nidentity to nearest strain reference: min={min(ids):.1f}% "
      f"median={sorted(ids)[len(ids)//2]:.1f}% max={max(ids):.1f}%")
print("nearest strain (fragment evidence only):")
for k, v in Counter(r["nearest_strain"] for r in rows).most_common():
    print(f"  strain {k}: {v}")
print(f"\nfragments below the 94% strain threshold against every reference: "
      f"{sum(1 for v in ids if v < 94)}")
print(f"fragments whose top two strain matches differ by <1%: "
      f"{sum(1 for r in rows if r['margin'] < 1)}")
print(f"fragments with a confident strain call: "
      f"{sum(1 for r in rows if r['assignment_confident']=='yes')}/{len(rows)}")

print("\nlowest-identity fragments:")
for r in sorted(rows, key=lambda r: r["identity_to_nearest"])[:6]:
    print(f"  {r['accession']} len={r['length']:>5} {r['identity_to_nearest']:6.2f}% "
          f"ref {r['ref_start']}-{r['ref_end']} {r['country']:<9s} {r['fragment_covers'][:44]}")

with open(os.path.join(OUT, "13b_partial_assignment.csv"), "w", newline="", encoding="utf-8-sig") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("\nwrote 13b_partial_assignment.csv")
