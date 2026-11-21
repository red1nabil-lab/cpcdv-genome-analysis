#!/usr/bin/env python
"""Step 34: regenerate Figure 2 from the PRIMARY identity matrix.

The original Figure 2 was drawn from the multiple-alignment, gaps-as-differences matrix
and so annotated 41 and 65 violating pairs, while the revised text reports 28 and 97
from true pairwise alignment with pairwise deletion. A figure disagreeing with its own
text is indefensible, so the figure is rebuilt from the matrix the text uses.

Panel a: distribution of all pairwise identities among the analysis-set genomes.
Panel b: within- vs between-strain identity for the genomes carrying a published
         single-letter strain designation.

Colour carries one meaning: the two distributions in panel b are the only categorical
encoding, and the threshold lines are a single neutral accent used for both panels.
"""
import os, sys, csv, json
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.stdout.reconfigure(encoding="utf-8")
import os as _os, sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # noqa: E402
from cpcdv_paths import DATA_DIR, FIG_DIR, MUSCLE_BIN  # noqa: E402
AN = DATA_DIR
OUTD = FIG_DIR

INK, INK2, INK3 = "#0b0b0b", "#3d3d3a", "#6b6a66"
SURFACE, GRID = "#fcfcfb", "#e8e7e2"
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6"]
ACCENT = "#6250d6"

cfg = json.load(open(os.path.join(AN, "10_analysis_set.json")))
keep = list(cfg["keep"])
M = np.load(os.path.join(AN, "25_truepw_pid_pw_A.npy")).astype(float)
if np.nanmax(M) <= 1.5:
    M *= 100.0
ix = {a: i for i, a in enumerate(keep)}
meta = {r["accession"]: r for r in
        csv.DictReader(open(os.path.join(AN, "03_curated_metadata.csv"),
                            encoding="utf-8-sig"))}
pub = {a: meta[a]["strain"] for a in keep
       if a in meta and len(meta[a]["strain"]) == 1 and meta[a]["strain"].isalpha()}
L = sorted(pub)
print(f"genomes {len(keep)}, labelled {len(L)}")

off = ~np.eye(len(keep), dtype=bool)
allvals = M[off]

within, between = [], []
for i, a in enumerate(L):
    for b in L[i + 1:]:
        (within if pub[a] == pub[b] else between).append(M[ix[a], ix[b]])
within = np.array(within); between = np.array(between)
n_breach = int((within < 94).sum())
n_merge = int((between > 94).sum())
print(f"within {len(within)} (below 94%: {n_breach}); "
      f"between {len(between)} (above 94%: {n_merge})")

fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.9),
                         gridspec_kw=dict(wspace=0.32, width_ratios=[1.15, 1]))
fig.patch.set_facecolor(SURFACE)

# ------------------------------------------------- panel a
ax = axes[0]; ax.set_facecolor(SURFACE)
ax.hist(allvals, bins=np.arange(74, 100.4, 0.4), color=CAT[0], zorder=3, rwidth=0.9)
for x, lab in ((94, "94%  strain"), (78, "78%  species")):
    ax.axvline(x, color=ACCENT, lw=1.3, ls=(0, (4, 2)), zorder=4)
    ax.text(x, ax.get_ylim()[1] * 0.96, lab, fontsize=6.8, color=ACCENT,
            rotation=90, va="top", ha="right")
ax.set_xlabel("pairwise identity (%)", fontsize=8, color=INK2)
ax.set_ylabel(f"genome pairs (n = {len(allvals)//2:,})", fontsize=8, color=INK2)
ax.set_title("a   all pairwise comparisons", loc="left", color=INK, fontsize=8.6)
ax.set_axisbelow(True); ax.yaxis.grid(True, color=GRID, lw=0.5)
ax.tick_params(labelsize=7, colors=INK3, length=2)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
for sp in ("left", "bottom"):
    ax.spines[sp].set_color(GRID)

# ------------------------------------------------- panel b
ax = axes[1]; ax.set_facecolor(SURFACE)
bp = ax.boxplot([within, between], widths=0.5, patch_artist=True, showfliers=True,
                medianprops=dict(color=SURFACE, lw=1.4),
                whiskerprops=dict(color=INK3, lw=0.7),
                capprops=dict(color=INK3, lw=0.7),
                flierprops=dict(marker="o", markersize=2.2, markerfacecolor=INK3,
                                markeredgecolor="none", alpha=0.45),
                boxprops=dict(lw=0), zorder=3)
for patch, c in zip(bp["boxes"], (CAT[2], CAT[1])):
    patch.set_facecolor(c)
ax.axhline(94, color=ACCENT, lw=1.3, ls=(0, (4, 2)), zorder=4)
ax.text(2.56, 94.4, "94%", fontsize=7, color=ACCENT, va="bottom")
ax.annotate(f"{n_breach} pairs\nbelow 94%", xy=(1.17, 89.5), fontsize=6.8,
            color=INK, va="center", ha="left")
ax.annotate(f"{n_merge} pairs\nabove 94%", xy=(2.18, 97.6), fontsize=6.8,
            color=INK, va="center", ha="left")
ax.set_xticks([1, 2])
ax.set_xticklabels([f"within strain\n(n = {len(within):,})",
                    f"between strains\n(n = {len(between):,})"], fontsize=7.4)
ax.set_ylabel("pairwise identity (%)", fontsize=8, color=INK2)
ax.set_title(f"b   strains as published (n = {len(L)} genomes)",
             loc="left", color=INK, fontsize=8.6)
ax.set_ylim(70, 101); ax.set_xlim(0.5, 3.15)
ax.set_axisbelow(True); ax.yaxis.grid(True, color=GRID, lw=0.5)
ax.tick_params(labelsize=7, colors=INK3, length=2)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
for sp in ("left", "bottom"):
    ax.spines[sp].set_color(GRID)

for ext in ("pdf", "png"):
    p = os.path.join(OUTD, f"Figure_2_pairwise_identity.{ext}")
    fig.savefig(p, dpi=300, bbox_inches="tight", facecolor=SURFACE)
    print("wrote", p)
