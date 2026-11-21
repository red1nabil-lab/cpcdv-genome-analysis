#!/usr/bin/env python
"""Step 37: are the breakpoints independent observations, or one event counted often?

Reviewer point 3.  Forty-five genomes share an inferred breakpoint at alignment column
1820.  If those genomes are close relatives, that is one ancestral recombination event
inherited by descendants, not 45 events, and the ~45-fold regional enrichment is
inflated by pseudo-replication.  Three checks:

  1. collapse genomes to distinct EVENTS, defined as a shared breakpoint window plus a
     shared pair of parental strains, and recount the regional enrichment on events
  2. test whether the genomes sharing the dominant breakpoint are phylogenetically
     clustered, which would indicate common descent rather than repeated events
  3. check whether the dominant columns sit in gappy or unusually variable alignment
     regions, which would suggest an alignment artefact
"""
import os, sys, json, csv
from collections import Counter, defaultdict
import numpy as np
from Bio import Phylo

sys.stdout.reconfigure(encoding="utf-8")
import os as _os, sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # noqa: E402
from cpcdv_paths import DATA_DIR, FIG_DIR, MUSCLE_BIN  # noqa: E402
AN = DATA_DIR
GAP = ord("-")
AMB = set(ord(c) for c in "RYSWKMBDHVN")

rec = [r for r in csv.DictReader(open(os.path.join(AN, "29_recombination_maxchi.csv"),
                                      encoding="utf-8-sig"))
       if r["significant"] == "True"]
print(f"genomes with supported recombination signal: {len(rec)}")

sass = {r["accession"]: r["strain_final"] for r in
        csv.DictReader(open(os.path.join(AN, "06_strain_assignment.csv"),
                            encoding="utf-8-sig"))}

# ------------------------------------------------- 1. collapse to distinct events
WIN = 50
events = defaultdict(list)
for r in rec:
    bp = int(r["breakpoint_aln_col"])
    key = (bp // WIN, tuple(sorted((r["parent1_strain"], r["parent2_strain"]))))
    events[key].append(r["accession"])
print(f"\n=== 1. collapsing genomes to distinct events "
      f"(same {WIN} nt window + same parental strain pair) ===")
print(f"  {len(rec)} genomes -> {len(events)} distinct events")
big = sorted(events.items(), key=lambda kv: -len(kv[1]))[:8]
for (wnd, par), accs in big:
    print(f"    window {wnd*WIN}-{wnd*WIN+WIN}  parents {par[0]}x{par[1]}  "
          f"{len(accs)} genomes")

REGIONS = [("LIR", 0, 181), ("V1 / MP", 181, 507), ("spacer", 507, 519),
           ("V2 / CP", 519, 1264), ("SIR", 1264, 1437), ("Rep exon 1", 1437, 1859),
           ("intron", 1859, 1945), ("Rep exon 2", 1945, 2632), ("LIR 3'", 2632, 2829)]


def tally(positions):
    c = Counter()
    for b in positions:
        for nm, s, e in REGIONS:
            if s <= b < e:
                c[nm] += 1
                break
    return c


per_genome = tally([int(r["breakpoint_aln_col"]) for r in rec])
per_event = tally([k[0] * WIN + WIN // 2 for k in events])
print(f"\n  {'region':14s} {'width':>6s} {'genomes':>8s} {'per kb':>8s} "
      f"{'events':>7s} {'per kb':>8s}")
tg, te = sum(per_genome.values()), sum(per_event.values())
for nm, s, e in REGIONS:
    w = e - s
    g, ev = per_genome.get(nm, 0), per_event.get(nm, 0)
    print(f"  {nm:14s} {w:6d} {g:8d} {1000*g/w:8.1f} {ev:7d} {1000*ev/w:8.1f}")
cod = sum(per_event.get(n, 0) for n in ("V1 / MP", "V2 / CP"))
codw = (507 - 181) + (1264 - 519)
hot = sum(per_event.get(n, 0) for n in ("SIR", "Rep exon 1"))
hotw = (1437 - 1264) + (1859 - 1437)
print(f"\n  enrichment on EVENTS (SIR+Rep exon 1 vs MP+CP): "
      f"{(hot/hotw)/(cod/codw) if cod else float('inf'):.1f}-fold")
cod_g = sum(per_genome.get(n, 0) for n in ("V1 / MP", "V2 / CP"))
hot_g = sum(per_genome.get(n, 0) for n in ("SIR", "Rep exon 1"))
print(f"  enrichment on GENOMES (the figure reported earlier):      "
      f"{(hot_g/hotw)/(cod_g/codw) if cod_g else float('inf'):.1f}-fold")

# ------------------------------------------------- 2. phylogenetic clustering
print("\n=== 2. are genomes sharing the dominant breakpoint close relatives? ===")
tr = Phylo.read(os.path.join(AN, "12_ml_tree.newick"), "newick")
tips = [t.name for t in tr.get_terminals()]
dom_key, dom_accs = big[0]
dom = [a for a in dom_accs if a in tips]
mrca = tr.common_ancestor(dom)
clade = [t.name for t in mrca.get_terminals()]
print(f"  dominant event: window {dom_key[0]*WIN}-{dom_key[0]*WIN+WIN}, "
      f"parents {dom_key[1][0]}x{dom_key[1][1]}, {len(dom)} genomes")
print(f"  their MRCA clade contains {len(clade)} tips, of which {len(dom)} carry the "
      f"breakpoint ({100*len(dom)/len(clade):.0f}%)")
print(f"  strains represented: {dict(Counter(sass.get(a,'?') for a in dom))}")
if len(clade) <= 2 * len(dom):
    print("  -> strongly clustered: consistent with ONE ancestral event inherited by "
          "descendants, not independent events")
else:
    print("  -> dispersed across the tree: more consistent with repeated independent "
          "events")

# ------------------------------------------------- 3. alignment artefact check
print("\n=== 3. do the dominant columns sit in gappy or hyper-variable regions? ===")


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
isbase = (A != GAP)
for c in AMB:
    isbase &= (A != c)
gapfrac = 1.0 - isbase.mean(axis=0)
nalleles = np.array([len(np.unique(A[:, c][isbase[:, c]])) if isbase[:, c].any() else 0
                     for c in range(A.shape[1])])
print(f"  genome-wide: mean gap fraction {gapfrac.mean():.3f}, "
      f"polymorphic columns {100*(nalleles>1).mean():.1f}%")
for col in [k[0] * WIN + WIN // 2 for k, _ in big[:3]]:
    lo, hi = max(0, col - 25), min(A.shape[1], col + 25)
    print(f"  around column {col:4d}: gap fraction {gapfrac[lo:hi].mean():.3f}, "
          f"polymorphic {100*(nalleles[lo:hi]>1).mean():.0f}%")

json.dump({"n_genomes": len(rec), "n_distinct_events": len(events),
           "per_genome_region": dict(per_genome), "per_event_region": dict(per_event),
           "dominant_event_genomes": len(dom), "dominant_event_clade_size": len(clade)},
          open(os.path.join(AN, "37_breakpoint_independence.json"), "w"), indent=1)
print("\nwrote 37_breakpoint_independence.json")
