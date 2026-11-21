#!/usr/bin/env python
"""Step 7: verify the divergent genomes by dedicated pairwise global alignment,
and place them against reference genomes of other dicot-infecting mastreviruses.

The matrix in step 6 comes from one 270-sequence progressive alignment, which can
understate identity for a divergent or short sequence.  Here each outlier is
realigned pairwise (Needleman-Wunsch, affine gaps) against every CpCDV genome and
against outgroup references, which is the defensible basis for a threshold claim.
"""
import csv, os, json, time
import numpy as np
from Bio import Entrez, SeqIO, Align

import os as _os, sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # noqa: E402
from cpcdv_paths import DATA_DIR, FIG_DIR, MUSCLE_BIN  # noqa: E402
OUT = DATA_DIR
Entrez.email = "plant-virology-client@example.org"
Entrez.tool = "cpcdv-diversity-analysis"

REFS = {           # other mastreviruses, for species-level context
    "NC_003493": "bean yellow dwarf virus (BeYDV)",
    "NC_014740": "chickpea chlorosis virus (CpCV)",
    "NC_023038": "chickpea redleaf virus (CpRLV)",
    "NC_003822": "tobacco yellow dwarf virus (TYDV)",
    "NC_001346": "maize streak virus (MSV, monocot outgroup)",
}
OUTLIERS = ["MN178605", "DQ458791", "AM849096", "HE610413", "KY817246",
            "KM229911", "KM229910"]

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

seqs = read_fasta(os.path.join(OUT, "04_full_oriented.fasta"))
meta = {r["accession"]: r for r in csv.DictReader(open(os.path.join(OUT, "03_curated_metadata.csv"), encoding="utf-8-sig"))}

# ---- fetch reference genomes -------------------------------------------
cache = os.path.join(OUT, "07_refs.gb")
if not (os.path.exists(cache) and os.path.getsize(cache) > 5000):
    h = Entrez.efetch(db="nuccore", id=",".join(REFS), rettype="gb", retmode="text")
    open(cache, "w", encoding="utf8").write(h.read()); h.close()
refseqs = {}
for rec in SeqIO.parse(cache, "genbank"):
    refseqs[rec.id.split(".")[0]] = str(rec.seq).upper()
print("reference genomes:", {k: len(v) for k, v in refseqs.items()})

# ---- aligner ------------------------------------------------------------
al = Align.PairwiseAligner()
al.mode = "global"
al.match_score = 1; al.mismatch_score = -1
al.open_gap_score = -5; al.extend_gap_score = -1
al.target_end_gap_score = 0.0; al.query_end_gap_score = 0.0   # free terminal gaps

def identity(a, b):
    """return (pid_pairwise_deletion, pid_gaps_as_difference, aligned_len)"""
    aln = al.align(a, b)[0]
    s1, s2 = str(aln[0]), str(aln[1])
    core = [(x, y) for x, y in zip(s1, s2)]
    # strip leading/trailing columns where either sequence is absent
    while core and ("-" in core[0]):
        core.pop(0)
    while core and ("-" in core[-1]):
        core.pop()
    both = [(x, y) for x, y in core if x != "-" and y != "-"]
    m = sum(1 for x, y in both if x == y)
    pid_pw = m / len(both) if both else 0.0
    pid_gap = m / len(core) if core else 0.0
    return pid_pw * 100, pid_gap * 100, len(core)

rows = []
for q in OUTLIERS:
    qs = seqs[q]
    best = []
    for a, s in seqs.items():
        if a == q:
            continue
        pw, pg, L = identity(qs, s)
        best.append((pw, pg, L, a))
    best.sort(reverse=True)
    m = meta[q]
    print(f"\n=== {q}  ({m['isolate']}, {m['country']}, {m['host_common']}, "
          f"strain {m['strain'] or '-'}, {len(qs)} nt) ===")
    print("  closest CpCDV genomes (pairwise realignment):")
    for pw, pg, L, a in best[:5]:
        am = meta[a]
        print(f"    {a} pid_pw={pw:6.2f}%  pid_gap={pg:6.2f}%  alnlen={L}  "
              f"strain={am['strain'] or '-':<3s} {am['country']:<12s} {am['isolate'][:24]}")
    for rk, name in REFS.items():
        pw, pg, L = identity(qs, refseqs[rk])
        print(f"    [ref] {rk} pid_pw={pw:6.2f}% pid_gap={pg:6.2f}%  {name}")
        rows.append(dict(query=q, subject=rk, subject_label=name,
                         pid_pw=round(pw, 2), pid_gap=round(pg, 2), aln_len=L,
                         comparison="outgroup reference"))
    pw, pg, L, a = best[0]
    rows.append(dict(query=q, subject=a, subject_label=f"closest CpCDV genome "
                     f"(strain {meta[a]['strain'] or '-'}, {meta[a]['country']})",
                     pid_pw=round(pw, 2), pid_gap=round(pg, 2), aln_len=L,
                     comparison="closest CpCDV genome"))

with open(os.path.join(OUT, "07_outlier_verification.csv"), "w", newline="", encoding="utf-8-sig") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("\nwrote 07_outlier_verification.csv")
