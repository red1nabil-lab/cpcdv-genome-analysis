#!/usr/bin/env python
"""Step 33: Figure 6 - recombination breakpoints along the CpCDV genome.

Two panels, sharing one x axis in alignment coordinates:
  (a) genome map, so the reader can see which feature each peak sits in
  (b) breakpoint count per 50 nt window, with the hotspot regions direct-labelled

Colour carries one meaning only: a single sequential hue for breakpoint count in
panel b, and the categorical palette in panel a for genome features. No hue means two
things across the figure, and the regions that matter are labelled directly rather than
relegated to a legend.
"""
import os, sys, csv, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

sys.stdout.reconfigure(encoding="utf-8")
import os as _os, sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # noqa: E402
from cpcdv_paths import DATA_DIR, FIG_DIR, MUSCLE_BIN  # noqa: E402
AN = DATA_DIR
OUTD = FIG_DIR
os.makedirs(OUTD, exist_ok=True)

INK, INK2, INK3 = "#0b0b0b", "#3d3d3a", "#6b6a66"
SURFACE, GRIDC = "#fcfcfb", "#e8e7e2"
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6"]
SEQ = "#2a78d6"

rec = list(csv.DictReader(open(os.path.join(AN, "29_recombination_maxchi.csv"),
                               encoding="utf-8-sig")))
bps = [int(r["breakpoint_aln_col"]) for r in rec if r["significant"] == "True"]
print(f"breakpoints plotted: {len(bps)}")

NCOL = 2829
GENES = [("V1\nMP", 181, 507, CAT[0]),
         ("V2  capsid protein", 519, 1264, CAT[1]),
         ("Rep exon 1", 1437, 1859, CAT[2]),
         ("Rep exon 2", 1945, 2632, CAT[2])]
NONCODE = [("LIR", 0, 181), ("SIR", 1264, 1437), ("LIR", 2632, NCOL)]

fig = plt.figure(figsize=(7.2, 3.5))
gs = fig.add_gridspec(2, 1, height_ratios=[0.30, 1.0], hspace=0.42)
fig.patch.set_facecolor(SURFACE)

# the two panels must share the x axis exactly, or the reader cannot map a peak in
# panel b onto the feature above it in panel a
ax_b = fig.add_subplot(gs[1])

# ---------------------------------------------------- panel a: genome map
ax = fig.add_subplot(gs[0], sharex=ax_b)
ax.set_facecolor(SURFACE)
for nm, s, e in NONCODE:
    ax.add_patch(FancyBboxPatch((s, 0.26), e - s, 0.48, boxstyle="round,pad=0,rounding_size=4",
                                linewidth=0, facecolor="#d8d7d2"))
    ax.text((s + e) / 2, 0.5, nm, ha="center", va="center", fontsize=6.4, color=INK2)
for nm, s, e, c in GENES:
    ax.add_patch(FancyBboxPatch((s, 0.22), e - s, 0.56, boxstyle="round,pad=0,rounding_size=4",
                                linewidth=0, facecolor=c))
    ax.text((s + e) / 2, 0.5, nm, ha="center", va="center", fontsize=6.4,
            color="white", linespacing=0.9)
ax.set_ylim(0, 1)
ax.set_yticks([])
ax.tick_params(axis="x", labelbottom=False, length=0)
for sp in ax.spines.values():
    sp.set_visible(False)
ax.set_title("a   genome organisation (origin-normalised alignment)",
             loc="left", fontsize=8.2, color=INK, pad=3)

# ---------------------------------------------------- panel b: breakpoints
ax = ax_b
ax.set_facecolor(SURFACE)
W = 50
edges = np.arange(0, NCOL + W, W)
cnt, _ = np.histogram(bps, bins=edges)
centres = edges[:-1] + W / 2
norm = cnt / cnt.max() if cnt.max() else cnt
for x, h, f in zip(centres, cnt, norm):
    if h == 0:
        continue
    # one hue, light to dark with magnitude
    col = plt.matplotlib.colors.to_rgb(SEQ)
    shade = tuple(1 - (1 - c) * (0.35 + 0.65 * f) for c in col)
    ax.add_patch(FancyBboxPatch((x - W * 0.42, 0), W * 0.84, h,
                                boxstyle="round,pad=0,rounding_size=1.5",
                                linewidth=0, facecolor=shade, zorder=3))
for y in ax.get_yticks():
    ax.axhline(y, color=GRIDC, lw=0.6, zorder=0)
ax.set_xlim(0, NCOL); ax.set_ylim(0, max(cnt) * 1.28)
ax.set_xlabel("position in alignment (nt)", fontsize=7.6, color=INK2)
ax.set_ylabel(f"breakpoints per {W} nt", fontsize=7.6, color=INK2)
ax.tick_params(labelsize=7, colors=INK3, length=2)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
for sp in ("left", "bottom"):
    ax.spines[sp].set_color(GRIDC)

# Direct labels on the two hotspots. Each names the number of GENOMES and the number of
# distinct EVENTS behind it: the tall peak is 38 strain-H genomes that share one
# ancestral event, so counting genomes would overstate the enrichment threefold.
top = list(np.argsort(cnt)[::-1][:2])
EVENT_NOTE = {0: "Rep exon 1\n38 genomes = 1 event (strain H)",
              1: "SIR\n19 genomes = 1 event"}
for rank, t in enumerate(top):
    left = rank == 1
    ax.annotate(EVENT_NOTE[rank],
                xy=(centres[t], cnt[t]),
                xytext=(centres[t] + (-150 if left else 140), cnt[t] * 1.00),
                fontsize=6.6, color=INK, va="center", linespacing=1.25,
                ha="right" if left else "left",
                arrowprops=dict(arrowstyle="-", color=INK3, lw=0.6))
ax.set_title("b   inferred recombination breakpoints "
             f"({len(bps)} genomes, 42 distinct events)",
             loc="left", fontsize=8.2, color=INK, pad=3)

for ext in ("pdf", "png"):
    p = os.path.join(OUTD, f"Figure_6_recombination_breakpoints.{ext}")
    fig.savefig(p, dpi=300, bbox_inches="tight", facecolor=SURFACE)
    print("wrote", p)
