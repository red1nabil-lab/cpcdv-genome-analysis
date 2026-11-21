#!/usr/bin/env python
"""Step 10: redundancy screen, gene coordinates, and population-genetic statistics.

Analysis set = 270 full genomes, less
  - the bean yellow dwarf virus lineage (verified in step 8/9), and
  - RefSeq records that duplicate a primary submission.
"""
import csv, os, json, itertools
import numpy as np
from collections import Counter, defaultdict
from Bio import SeqIO

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
ori = read_fasta(os.path.join(OUT, "04_full_oriented.fasta"))
meta = {r["accession"]: r for r in csv.DictReader(open(os.path.join(OUT, "03_curated_metadata.csv"), encoding="utf-8-sig"))}
strainf = {r["accession"]: r["strain_final"] for r in csv.DictReader(open(os.path.join(OUT, "06_strain_assignment.csv"), encoding="utf-8-sig"))}
pid = np.load(os.path.join(OUT, "06_pid_gap.npy"))
accs0 = json.load(open(os.path.join(OUT, "06_accession_order.json")))
idx0 = {a: i for i, a in enumerate(accs0)}

# ---------------------------------------------------- 1. redundancy screen
print("=== identical and near-identical records (>=99.9% identity) ===")
groups = []
seen = set()
for i, a in enumerate(accs0):
    if a in seen:
        continue
    g = [a]
    for j, b in enumerate(accs0):
        if j > i and pid[i, j] >= 0.999:
            g.append(b)
    if len(g) > 1:
        groups.append(g)
        seen.update(g)
for g in groups:
    print("  " + " = ".join(g))
    for a in g:
        m = meta[a]
        print(f"      {a}: {m['isolate'][:28]:<30s} {m['country']:<12s} "
              f"{m['collection_year']:<6s} {m['reference_title'][:52]}")

BEYDV = ["AM849096", "DQ458791"]
refseq_dupes = [a for g in groups for a in g
                if a.startswith("NC_") and any(not x.startswith("NC_") for x in g)]
print(f"\nbean yellow dwarf virus lineage removed: {BEYDV}")
print(f"RefSeq duplicates of a primary submission removed: {refseq_dupes}")

keep = [a for a in accs0 if a not in set(BEYDV) | set(refseq_dupes)]
print(f"analysis set: {len(keep)} genomes")

# ---------------------------------------------------- 2. gene coordinates
rec = next(r for r in SeqIO.parse(os.path.join(OUT, "02_genbank_cache.gb"), "genbank")
           if r.id.split(".")[0] == "NC_011058")
rot = ori["NC_011058"]
raw = str(rec.seq).upper()
print(f"\nNC_011058: raw {len(raw)} nt, rotated {len(rot)} nt")

def locate(sub):
    """position of sub in the rotated circular genome"""
    p = (rot + rot).find(sub)
    return p if p >= 0 else None

feats = []
for f in rec.features:
    if f.type not in ("CDS", "gene"):
        continue
    name = (f.qualifiers.get("gene") or f.qualifiers.get("product") or ["?"])[0]
    for part in (f.location.parts if hasattr(f.location, "parts") else [f.location]):
        s, e = int(part.start), int(part.end)
        seg = raw[s:e]
        p = locate(seg)
        if p is None:
            continue
        feats.append(dict(feature=f.type, name=name, strand=part.strand,
                          rot_start=p, rot_end=p + len(seg), length=len(seg)))
# deduplicate gene/CDS pairs
uniq = {}
for f in feats:
    k = (f["name"], f["rot_start"], f["rot_end"])
    if k not in uniq or f["feature"] == "CDS":
        uniq[k] = f
feats = sorted(uniq.values(), key=lambda f: f["rot_start"])
print("\nannotated features mapped onto the rotated genome:")
for f in feats:
    print(f"  {f['name']:<10s} {f['feature']:<4s} strand={f['strand']:>2}  "
          f"{f['rot_start']:>5}-{f['rot_end']:<5} ({f['length']} nt)")

# map rotated coords -> alignment columns via the NC_011058 row
refrow = aln["NC_011058"]
r2c, pos = {}, 0
for col, ch in enumerate(refrow):
    if ch != "-":
        r2c[pos] = col
        pos += 1
print(f"\nNC_011058 occupies {pos} of {len(refrow)} alignment columns")

REGIONS = []
for f in feats:
    s, e = f["rot_start"], min(f["rot_end"], pos) - 1
    if s in r2c and e in r2c:
        REGIONS.append((f"{f['name']} ({f['length']} nt)", r2c[s], r2c[e] + 1))
# intergenic regions: everything not covered by a feature
cov = np.zeros(len(refrow), dtype=bool)
for _, s, e in REGIONS:
    cov[s:e] = True
gaps_un = []
st = None
for c in range(len(refrow)):
    if not cov[c] and st is None:
        st = c
    elif cov[c] and st is not None:
        if c - st >= 60:
            gaps_un.append((st, c))
        st = None
if st is not None and len(refrow) - st >= 60:
    gaps_un.append((st, len(refrow)))
for s, e in gaps_un:
    REGIONS.append((f"intergenic {s}-{e}", s, e))

# ---------------------------------------------------- 3. diversity statistics
A = np.frombuffer(b"".join(aln[a].upper().encode() for a in keep), dtype=np.uint8).reshape(len(keep), len(refrow))
isb = A != ord("-")
for c in "RYSWKMBDHVN":
    isb &= (A != ord(c))

def stats(rows, lo=0, hi=None):
    """pi, S, haplotypes, Tajima's D for a set of rows over columns [lo,hi)"""
    hi = hi if hi is not None else A.shape[1]
    sub, sb = A[np.ix_(rows, range(lo, hi))], isb[np.ix_(rows, range(lo, hi))]
    n = len(rows)
    # columns usable in every sequence
    ok = sb.all(0)
    if ok.sum() < 50 or n < 2:
        return None
    S_ = sub[:, ok]
    L = ok.sum()
    seg = sum(1 for c in range(L) if len(set(S_[:, c].tolist())) > 1)
    # pairwise differences
    tot = cnt = 0
    for i, j in itertools.combinations(range(n), 2):
        tot += int((S_[i] != S_[j]).sum()); cnt += 1
    pi = tot / cnt / L if cnt else 0.0
    haps = len({S_[i].tobytes() for i in range(n)})
    # haplotype diversity
    freq = Counter(S_[i].tobytes() for i in range(n))
    hd = (n / (n - 1)) * (1 - sum((v / n) ** 2 for v in freq.values())) if n > 1 else 0.0
    # Tajima's D
    a1 = sum(1 / i for i in range(1, n)); a2 = sum(1 / i**2 for i in range(1, n))
    b1 = (n + 1) / (3 * (n - 1)); b2 = 2 * (n**2 + n + 3) / (9 * n * (n - 1))
    c1 = b1 - 1 / a1; c2 = b2 - (n + 2) / (a1 * n) + a2 / a1**2
    e1 = c1 / a1; e2 = c2 / (a1**2 + a2)
    k = tot / cnt if cnt else 0.0
    var = e1 * seg + e2 * seg * (seg - 1)
    D = (k - seg / a1) / np.sqrt(var) if var > 0 else float("nan")
    return dict(n=n, sites=int(L), S=seg, pi=pi, k=k, haplotypes=haps,
                hap_diversity=hd, tajimaD=D)

rowsall = list(range(len(keep)))
print("\n=== whole-genome diversity, full analysis set ===")
s = stats(rowsall)
print(f"  n={s['n']}  sites={s['sites']}  segregating={s['S']}  "
      f"pi={s['pi']:.4f}  mean pairwise diff k={s['k']:.1f}  "
      f"haplotypes={s['haplotypes']}  Hd={s['hap_diversity']:.4f}  "
      f"Tajima's D={s['tajimaD']:.3f}")

out_rows = [dict(group="all genomes", level="whole genome", **s)]

def group_stats(key, label):
    print(f"\n=== diversity by {label} ===")
    g = defaultdict(list)
    for i, a in enumerate(keep):
        v = meta[a][key] if key in meta[a] else strainf[a]
        g[v].append(i)
    for v, rows in sorted(g.items(), key=lambda kv: -len(kv[1])):
        st = stats(rows)
        if st is None:
            print(f"  {v or '(blank)':<28s} n={len(rows):<4} (too few)")
            continue
        print(f"  {v or '(blank)':<28s} n={st['n']:<4} S={st['S']:<5} "
              f"pi={st['pi']:.4f}  Hd={st['hap_diversity']:.3f}  "
              f"haps={st['haplotypes']:<4} D={st['tajimaD']:+.3f}")
        out_rows.append(dict(group=f"{label}: {v}", level="whole genome", **st))

group_stats("country", "country")
group_stats("subregion", "subregion")
group_stats("host_family", "host family")

print("\n=== diversity by assigned strain ===")
g = defaultdict(list)
for i, a in enumerate(keep):
    g[strainf[a]].append(i)
for v, rows in sorted(g.items(), key=lambda kv: -len(kv[1])):
    st = stats(rows)
    if st is None:
        print(f"  strain {v:<12s} n={len(rows):<4} (too few)")
        continue
    print(f"  strain {v:<12s} n={st['n']:<4} S={st['S']:<5} pi={st['pi']:.4f}  "
          f"Hd={st['hap_diversity']:.3f}  D={st['tajimaD']:+.3f}")
    out_rows.append(dict(group=f"strain: {v}", level="whole genome", **st))

print("\n=== diversity by genome region ===")
for name, lo, hi in REGIONS:
    st = stats(rowsall, lo, hi)
    if st is None:
        print(f"  {name:<34s} (too few usable columns)")
        continue
    print(f"  {name:<34s} cols {lo:>5}-{hi:<5} sites={st['sites']:<5} "
          f"S={st['S']:<5} pi={st['pi']:.4f} D={st['tajimaD']:+.3f}")
    out_rows.append(dict(group=f"region: {name}", level="genome region", **st))

with open(os.path.join(OUT, "10_diversity_statistics.csv"), "w", newline="", encoding="utf-8-sig") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys())); w.writeheader(); w.writerows(out_rows)
json.dump(dict(keep=keep, regions=REGIONS, excluded_beydv=BEYDV,
               excluded_refseq_dupes=refseq_dupes),
          open(os.path.join(OUT, "10_analysis_set.json"), "w"), indent=1)
print("\nwrote 10_diversity_statistics.csv, 10_analysis_set.json")
