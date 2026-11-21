#!/usr/bin/env python
"""Step 4: check strand and rotation of the circular genomes, normalise the origin.

Mastrevirus genomes carry the conserved nonanucleotide TAATATTAC in the hairpin
loop of the large intergenic region; the Rep nick falls between T7 and A8.  To
compare circular molecules sensibly every genome is reverse-complemented if
necessary and rotated to begin at the A of ...TAATATT|AC.  Variant
nonanucleotides are accepted at up to one mismatch; any genome in which no
variant can be located is rotated instead by the best ungapped match to a
40-nt probe taken from the oriented RefSeq.
"""
import csv, os, re
from collections import Counter

import os as _os, sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # noqa: E402
from cpcdv_paths import DATA_DIR, FIG_DIR, MUSCLE_BIN  # noqa: E402
OUT = DATA_DIR
MOTIF = "TAATATTAC"
REFACC = "NC_011058"
RC = str.maketrans("ACGTRYSWKMBDHVNacgtryswkmbdhvn", "TGCAYRSWMKVHDBNtgcayrswmkvhdbn")

def rc(s):
    return s.translate(RC)[::-1]

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

def find_motif(s, maxmm):
    """positions of MOTIF in circular s with <= maxmm mismatches, best first"""
    dbl = s + s
    hits = []
    for i in range(len(s)):
        mm = sum(1 for a, b in zip(dbl[i:i+9], MOTIF) if a != b)
        if mm <= maxmm:
            hits.append((mm, i))
    return sorted(hits)

seqs = read_fasta(os.path.join(OUT, "01_all_dedup.fasta"))
meta = {r["accession"]: r for r in csv.DictReader(open(os.path.join(OUT, "03_curated_metadata.csv"), encoding="utf-8-sig"))}
full = {a: s.upper() for a, s in seqs.items() if meta[a]["molecule_class"] == "full genome"}
print(f"{len(full)} full genomes to orient")

report, oriented, unresolved = [], {}, []
for a, s in full.items():
    f, r = find_motif(s, 1), find_motif(rc(s), 1)
    if f and (not r or f[0][0] <= r[0][0]):
        strand, (mm, pos) = "+", f[0]
    elif r:
        strand, (mm, pos) = "-", r[0]
    else:
        strand, mm, pos = "?", -1, -1
        unresolved.append(a)
    motif_seen = ""
    if strand != "?":
        base = s if strand == "+" else rc(s)
        motif_seen = (base + base)[pos:pos+9]
        off = (pos + 7) % len(base)              # nick between T7 and A8
        oriented[a] = base[off:] + base[:off]
    report.append(dict(accession=a, length=len(s), strand=strand, motif_pos=pos,
                       motif_mismatches=mm, motif_observed=motif_seen,
                       method="nonanucleotide"))

# ---- fallback for genomes with no locatable nonanucleotide ---------------
if unresolved:
    probe = oriented[REFACC][:40]
    for a in unresolved:
        best = None
        for cand, lab in ((full[a], "+"), (rc(full[a]), "-")):
            dbl = cand + cand
            for i in range(len(cand)):
                sc = sum(1 for x, y in zip(dbl[i:i+40], probe) if x == y)
                if best is None or sc > best[0]:
                    best = (sc, i, lab, cand)
        sc, i, lab, cand = best
        oriented[a] = cand[i:] + cand[:i]
        for row in report:
            if row["accession"] == a:
                row.update(strand=lab, motif_pos=i, motif_mismatches="",
                           motif_observed=f"{sc}/40 identity to RefSeq probe",
                           method="RefSeq probe alignment")
        print(f"  {a}: rotated by RefSeq probe, strand {lab}, offset {i}, {sc}/40 match")

print("\nstrand as submitted:", dict(Counter(r["strand"] for r in report)))
print("nonanucleotide variants observed:",
      dict(Counter(r["motif_observed"] for r in report if len(str(r["motif_observed"])) == 9)))
pp = [r["motif_pos"] for r in report if r["motif_pos"] >= 0]
print(f"origin position before rotation: min={min(pp)} max={max(pp)} distinct={len(set(pp))}")
print("first 4 bases after rotation:", Counter(s[:4] for s in oriented.values()).most_common(6))

with open(os.path.join(OUT, "04_orientation_report.csv"), "w", newline="", encoding="utf-8-sig") as fh:
    w = csv.DictWriter(fh, fieldnames=list(report[0].keys())); w.writeheader(); w.writerows(report)
with open(os.path.join(OUT, "04_full_oriented.fasta"), "w", encoding="utf8") as fh:
    for a in sorted(oriented):
        fh.write(f">{a}\n")
        s = oriented[a]
        for i in range(0, len(s), 60):
            fh.write(s[i:i+60] + "\n")
print(f"\nwrote 04_full_oriented.fasta ({len(oriented)} genomes), 04_orientation_report.csv")
