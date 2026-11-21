#!/usr/bin/env python
"""Step 8: place the divergent genomes in a species context using a MUSCLE
alignment of CpCDV strain representatives plus reference dicot-infecting
mastreviruses.  MUSCLE + pairwise identity is the basis of the SDT procedure
used to set the 78%/94% Mastrevirus thresholds, so identities are computed the
same way here rather than from an ad hoc pairwise aligner.
"""
import csv, os, subprocess, json
import numpy as np
from Bio import Entrez, SeqIO

import os as _os, sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # noqa: E402
from cpcdv_paths import DATA_DIR, FIG_DIR, MUSCLE_BIN  # noqa: E402
OUT = DATA_DIR
MUSCLE = MUSCLE_BIN
Entrez.email = "plant-virology-client@example.org"
Entrez.tool = "cpcdv-diversity-analysis"

# reference genomes, verified by their GenBank ORGANISM line after download
WANT = ["NC_003493",   # Bean yellow dwarf virus
        "NC_014740",   # Chickpea chlorosis virus-A
        "NC_003822",   # Tobacco yellow dwarf virus
        "NC_001346",   # Maize streak virus (monocot outgroup)
        "NC_038901",   # candidate: chickpea redleaf virus
        "NC_038902"]   # candidate: further dicot mastrevirus

cache = os.path.join(OUT, "08_refs.gb")
if not (os.path.exists(cache) and os.path.getsize(cache) > 5000):
    h = Entrez.efetch(db="nuccore", id=",".join(WANT), rettype="gb", retmode="text")
    open(cache, "w", encoding="utf8").write(h.read()); h.close()
refs = {}
for rec in SeqIO.parse(cache, "genbank"):
    acc = rec.id.split(".")[0]
    org = rec.annotations.get("organism", "?")
    refs[acc] = (org, str(rec.seq).upper())
    print(f"  {acc}  {len(rec.seq):>5} bp  {org}")

# keep only mastrevirus genomes of plausible size
refs = {k: v for k, v in refs.items() if 2400 <= len(v[1]) <= 2900}
print("\nreferences retained:", {k: v[0] for k, v in refs.items()})

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
strainf = {r["accession"]: r["strain_final"] for r in csv.DictReader(open(os.path.join(OUT, "06_strain_assignment.csv"), encoding="utf-8-sig"))}

# one representative per annotated strain (the longest genome), plus the outliers
reps = {}
for a in seqs:
    s = meta[a]["strain"]
    if len(s) == 1 and s.isalpha():
        if s not in reps or len(seqs[a]) > len(seqs[reps[s]]):
            reps[s] = a
OUTLIERS = ["MN178605", "DQ458791", "AM849096", "HE610413", "KY817246", "NC_011058"]
sel = {}
for s, a in sorted(reps.items()):
    sel[f"CpCDV_strain{s}_{a}"] = seqs[a]
for a in OUTLIERS:
    lab = f"QUERY_{a}_{(meta[a]['strain'] or 'nostrain')}_{meta[a]['country'].replace(' ','')}"
    sel[lab] = seqs[a]
for k, (org, s) in refs.items():
    sel[f"REF_{k}_{org.replace(' ', '')[:28]}"] = s

fin = os.path.join(OUT, "08_species_set.fasta")
with open(fin, "w", encoding="utf8") as fh:
    for k, s in sel.items():
        fh.write(f">{k}\n")
        for i in range(0, len(s), 60):
            fh.write(s[i:i+60] + "\n")
print(f"\n{len(sel)} sequences in the species-context set")

fout = os.path.join(OUT, "08_species_aln.fasta")
subprocess.run([MUSCLE, "-in", fin, "-out", fout, "-maxiters", "16", "-quiet"], check=True)
aln = read_fasta(fout)
labs = list(sel)
L = len(next(iter(aln.values())))
print(f"alignment length {L}")

A = np.frombuffer(b"".join(aln[k].upper().encode() for k in labs), dtype=np.uint8).reshape(len(labs), L)
GAP = ord("-")
isb = A != GAP
for c in "RYSWKMBDHVN":
    isb &= (A != ord(c))

n = len(labs)
pw = np.zeros((n, n)); pg = np.zeros((n, n))
for i in range(n):
    same = A == A[i]
    both = isb & isb[i]; either = isb | isb[i]
    m = (same & both).sum(1)
    pw[i] = m / np.maximum(both.sum(1), 1)
    pg[i] = m / np.maximum(either.sum(1), 1)

rows = []
for i, a in enumerate(labs):
    for j, b in enumerate(labs):
        if i < j:
            rows.append(dict(seq_a=a, seq_b=b, pid_pairwise_deletion=round(pw[i, j]*100, 2),
                             pid_gaps_as_difference=round(pg[i, j]*100, 2)))
with open(os.path.join(OUT, "08_species_identity.csv"), "w", newline="", encoding="utf-8-sig") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

qi = [i for i, k in enumerate(labs) if k.startswith("QUERY_")]
ri = [i for i, k in enumerate(labs) if k.startswith("REF_")]
ci = [i for i, k in enumerate(labs) if k.startswith("CpCDV_")]
print("\n=== each query against CpCDV strain representatives and references "
      "(pid_pairwise_deletion %) ===")
for i in qi:
    print(f"\n{labs[i]}  ({len(sel[labs[i]])} nt)")
    best = sorted(((pw[i, j], labs[j]) for j in ci), reverse=True)[:3]
    for v, lab in best:
        print(f"    best CpCDV: {lab:<34s} {v*100:6.2f}%")
    for j in ri:
        print(f"    reference : {labs[j]:<34s} {pw[i, j]*100:6.2f}%")

print("\n=== CpCDV RefSeq vs references, for threshold context ===")
k = labs.index([l for l in labs if l.startswith("QUERY_NC_011058")][0])
for j in ri:
    print(f"  NC_011058 (CpCDV) vs {labs[j]:<34s} {pw[k, j]*100:6.2f}% "
          f"(gaps-as-diff {pg[k, j]*100:.2f}%)")
print("\nCpCDV strain representatives: lowest pairwise identity within CpCDV")
sub = pw[np.ix_(ci, ci)]
m = ~np.eye(len(ci), dtype=bool)
print(f"  min={sub[m].min()*100:.2f}%  max={sub[m].max()*100:.2f}%")
print("wrote 08_species_identity.csv")
