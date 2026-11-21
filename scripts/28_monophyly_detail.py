#!/usr/bin/env python
"""Step 28: which accessions break strain monophyly, and with what node support?

Reviewer criticism 5: the manuscript reports strains C, H and L as non-monophyletic but
does not identify the responsible genomes nor the support of the nodes involved.  A
non-monophyly resting on a node with 64% bootstrap support is not evidence of anything.

Two things also have to be got right, and the first version of this script got the
second one wrong:

  (a) Monophyly is ROOT-DEPENDENT.  An ML tree is unrooted, so the primary test must be
      rooting-independent: a strain forms a clade iff its tip set is one side of a
      bipartition induced by some internal edge.  Reading the unrooted Newick as though
      it were rooted makes any group straddling the arbitrary root look polyphyletic
      (this is what made strain K appear non-monophyletic).
  (b) Support must be reported for the specific nodes that create the conflict.
"""
import os, sys, csv
from collections import Counter, defaultdict
from Bio import Phylo

sys.stdout.reconfigure(encoding="utf-8")
import os as _os, sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # noqa: E402
from cpcdv_paths import DATA_DIR, FIG_DIR, MUSCLE_BIN  # noqa: E402
AN = DATA_DIR

sass = {r["accession"]: r["strain_final"] for r in
        csv.DictReader(open(os.path.join(AN, "06_strain_assignment.csv"),
                            encoding="utf-8-sig"))}
md = {r["accession"]: r for r in
      csv.DictReader(open(os.path.join(AN, "03_curated_metadata.csv"),
                          encoding="utf-8-sig"))}


def support(clade):
    v = clade.confidence
    if v is None:
        try:
            v = float(clade.name)
        except (TypeError, ValueError):
            return None
    return v / 100.0 if v > 1.5 else v


def analyse(path, tag):
    tr = Phylo.read(path, "newick")
    tips = [t.name for t in tr.get_terminals()]
    allt = frozenset(tips)
    strain = {t: sass.get(t, "unassigned") for t in tips}
    by = defaultdict(list)
    for t, s in strain.items():
        by[s].append(t)

    # --- rooting-independent bipartitions
    splits = {}
    for nd in tr.get_nonterminals():
        half = frozenset(t.name for t in nd.get_terminals())
        if 1 < len(half) < len(tips):
            splits[half] = support(nd)

    print("\n" + "=" * 78)
    print(f"{tag}: {len(tips)} tips, {len(splits)} informative internal edges")
    rows = []
    for s in sorted(by, key=lambda x: -len(by[x])):
        mem = by[s]
        if s == "unassigned" or len(mem) < 2:
            continue
        S = frozenset(mem)
        comp = allt - S
        clade_unrooted = (S in splits) or (comp in splits)
        sup_split = splits.get(S, splits.get(comp))

        mrca = tr.common_ancestor(mem)
        ctips = [t.name for t in mrca.get_terminals()]
        intr = [t for t in ctips if strain.get(t) != s]
        straddles_root = len(ctips) == len(tips) and not clade_unrooted

        verdict = ("clade (rooting-independent)" if clade_unrooted else
                   "NOT a clade in any rooting")
        print(f"\nstrain {s:3s} n={len(mem):3d}  {verdict}")
        if clade_unrooted:
            print(f"   bipartition support = "
                  f"{('%.2f' % sup_split) if sup_split is not None else 'n/a'}")
        else:
            print(f"   MRCA clade holds {len(ctips)} tips; "
                  f"{len(intr)} non-members nested inside")
            print(f"   intruders by assigned strain: "
                  f"{dict(Counter(strain[t] for t in intr))}")
            for t in sorted(intr, key=lambda x: strain[x]):
                path_ = tr.get_path(t)
                ns = None
                for nd in reversed(path_[:-1]):
                    ns = support(nd)
                    if ns is not None:
                        break
                print(f"      {t:10s} strain={strain[t]:10s} "
                      f"published={md.get(t,{}).get('strain','-') or '-':3s} "
                      f"{md.get(t,{}).get('country','?'):13s} "
                      f"node support={('%.2f' % ns) if ns is not None else 'n/a'}")
        rows.append(dict(strain=s, n=len(mem),
                         clade_rooting_independent=clade_unrooted,
                         bipartition_support=("" if sup_split is None
                                              else round(sup_split, 3)),
                         mrca_clade_size=len(ctips),
                         n_intruders=(0 if clade_unrooted else len(intr)),
                         intruders=("" if clade_unrooted else ";".join(sorted(intr))),
                         straddles_root=straddles_root))
    sups = [v for v in splits.values() if v is not None]
    print(f"\n   support distribution over {len(sups)} internal edges:")
    for thr in (0.95, 0.90, 0.70, 0.50):
        print(f"     >= {thr:.2f}: {sum(1 for v in sups if v >= thr):4d} "
              f"({100*sum(1 for v in sups if v>=thr)/len(sups):5.1f}%)")
    return rows


rows = analyse(os.path.join(AN, "12_ml_tree.newick"), "ML tree (as inferred, unrooted)")
out = os.path.join(AN, "28_monophyly_detail.csv")
with open(out, "w", newline="", encoding="utf-8-sig") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
print(f"\nwrote {out}")
