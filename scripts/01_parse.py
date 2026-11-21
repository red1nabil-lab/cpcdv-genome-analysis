#!/usr/bin/env python
"""Step 1: parse the CpCDV FASTA, dedupe, classify, write a clean inventory."""
import re, csv, json, sys, os

FASTA = os.path.join(DATA_DIR, "CpCDV_all_seq.fasta")
import os as _os, sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # noqa: E402
from cpcdv_paths import DATA_DIR, FIG_DIR, MUSCLE_BIN  # noqa: E402
OUT = DATA_DIR
os.makedirs(OUT, exist_ok=True)

recs = []
acc = None
hdr = None
seq = []
with open(FASTA, encoding="utf8", errors="replace") as fh:
    for line in fh:
        line = line.rstrip("\n\r")
        if line.startswith(">"):
            if acc is not None:
                recs.append((acc, hdr, "".join(seq)))
            hdr = line[1:].strip()
            acc = hdr.split()[0]
            seq = []
        else:
            seq.append(line.strip())
if acc is not None:
    recs.append((acc, hdr, "".join(seq)))

print(f"records parsed: {len(recs)}")

# --- duplicate accessions -------------------------------------------------
from collections import Counter, defaultdict
c = Counter(a for a, _, _ in recs)
dups = {a: n for a, n in c.items() if n > 1}
print("duplicate accessions:", dups)

by_acc = defaultdict(list)
for a, h, s in recs:
    by_acc[a].append((h, s))
for a, lst in by_acc.items():
    if len(lst) > 1:
        seqs = {s for _, s in lst}
        print(f"  {a}: {len(lst)} entries, {len(seqs)} distinct sequence(s), "
              f"lengths={[len(s) for _, s in lst]}")

# keep the longest entry per accession
uniq = {}
for a, lst in by_acc.items():
    h, s = max(lst, key=lambda t: len(t[1]))
    uniq[a] = (h, s)
print(f"unique accessions: {len(uniq)}")

# --- parse header fields --------------------------------------------------
def parse_header(acc, h):
    rest = h[len(acc):].strip()
    # strip the "-CpCDV-" marker
    rest = re.sub(r"^-?\s*CpCDV\s*-?", "", rest).strip()
    strain = None
    isolate = None
    clone = None
    m = re.search(r"\bstrain\s+([A-Za-z0-9_\-/\.\[\]:]+)", rest)
    if m:
        strain = m.group(1).rstrip(",")
    m = re.search(r"\bisolate\s+(.+?)(?:\s+(?:capsid|coat|movement|replication|RepA|Rep,|L RepA)\b|,\s*complete|,\s*$|$)", rest)
    if m:
        isolate = m.group(1).strip().rstrip(",")
    m = re.search(r"\bclone\s+([A-Za-z0-9_\-/\.]+)", rest)
    if m:
        clone = m.group(1).rstrip(",")
    # molecule completeness
    low = rest.lower()
    if "complete genome" in low or "complete sequence" in low:
        completeness = "complete genome"
    elif "partial cds" in low:
        completeness = "partial"
    else:
        completeness = "unstated"
    # which genes the record covers
    genes = []
    if re.search(r"capsid protein|coat protein|\bCp\b|\bV2\b|\bV1\b|av1", rest, re.I):
        genes.append("CP")
    if re.search(r"movement protein|\bMP\b", rest, re.I):
        genes.append("MP")
    if re.search(r"RepA|replication-associated|\bRep\b|\bC1\b", rest, re.I):
        genes.append("Rep")
    return strain, isolate, clone, completeness, "+".join(genes) or None, rest

rows = []
for a, (h, s) in uniq.items():
    strain, isolate, clone, completeness, genes, rest = parse_header(a, h)
    L = len(s)
    # length-based classification is the reliable one
    if L >= 2400:
        cls = "full genome"
    elif L >= 1000:
        cls = "large fragment"
    else:
        cls = "partial (gene fragment)"
    amb = sum(1 for ch in s.upper() if ch not in "ACGT")
    rows.append(dict(accession=a, strain=strain or "", isolate=isolate or "",
                     clone=clone or "", completeness=completeness,
                     genes=genes or "", length=L, class_=cls,
                     ambiguous=amb, header_rest=rest))

rows.sort(key=lambda r: (-r["length"], r["accession"]))

with open(os.path.join(OUT, "01_inventory_raw.csv"), "w", newline="", encoding="utf8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

print("\nlength classes:")
for k, v in Counter(r["class_"] for r in rows).most_common():
    print(f"  {k}: {v}")
print("\nlength distribution:")
import statistics
L = sorted(r["length"] for r in rows)
print(f"  min={L[0]} q1={L[len(L)//4]} median={statistics.median(L)} "
      f"q3={L[3*len(L)//4]} max={L[-1]}")
print("\nstrain labels present in headers:")
for k, v in Counter(r["strain"] for r in rows if r["strain"]).most_common():
    print(f"  {k}: {v}")
print("\nrecords with ambiguous bases:",
      sum(1 for r in rows if r["ambiguous"] > 0))

# write the accession list for the GenBank fetch
with open(os.path.join(OUT, "01_accessions.txt"), "w", encoding="utf8") as fh:
    for r in sorted(rows, key=lambda r: r["accession"]):
        fh.write(r["accession"] + "\n")

# write a deduplicated fasta with clean single-token ids
with open(os.path.join(OUT, "01_all_dedup.fasta"), "w", encoding="utf8") as fh:
    for r in rows:
        fh.write(f">{r['accession']}\n")
        s = uniq[r["accession"]][1]
        for i in range(0, len(s), 60):
            fh.write(s[i:i+60] + "\n")

print(f"\nwrote {OUT}\\01_inventory_raw.csv, 01_accessions.txt, 01_all_dedup.fasta")
