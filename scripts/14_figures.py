#!/usr/bin/env python
"""Step 14: publication figures for the CpCDV diversity analysis.

Print figures: light surface only, 300 dpi (PDF vector + PNG raster).
Categorical hues are taken in fixed slot order from the validated palette and
are never cycled; a ninth category folds into "other".  Single-series panels
use one hue with direct value labels, which also supplies the relief the
palette's contrast warning requires.
"""
import csv, os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from collections import Counter, defaultdict, OrderedDict

import os as _os, sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # noqa: E402
from cpcdv_paths import DATA_DIR, FIG_DIR, MUSCLE_BIN  # noqa: E402
OUT = DATA_DIR
FIG = os.path.join(OUT, "figures")
os.makedirs(FIG, exist_ok=True)

# ---- design tokens (validated palette, light mode) ----------------------
SURFACE = "#fcfcfb"
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8984"
GRID = "#e4e3df"
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4",
       "#008300", "#6250d6", "#e34948"]
OTHER = "#b8b7b1"
SEQ = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#2a78d6", "#256abf",
       "#1c5cab", "#184f95", "#104281", "#0d366b"]

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8, "axes.labelsize": 8.5, "axes.titlesize": 9.5,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5,
    "axes.edgecolor": INK3, "axes.linewidth": 0.6,
    "xtick.color": INK2, "ytick.color": INK2,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "axes.grid": False, "legend.frameon": False,
    "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), dpi=300,
                    bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    print("  wrote", name + ".pdf/.png")

def grid_x(ax):
    ax.set_axisbelow(True)
    ax.xaxis.grid(True, color=GRID, lw=0.5)

def grid_y(ax):
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color=GRID, lw=0.5)

# ---- data ---------------------------------------------------------------
meta = {r["accession"]: r for r in csv.DictReader(open(os.path.join(OUT, "03_curated_metadata.csv"), encoding="utf-8-sig"))}
sassign = {r["accession"]: r["strain_final"] for r in csv.DictReader(open(os.path.join(OUT, "06_strain_assignment.csv"), encoding="utf-8-sig"))}
cfg = json.load(open(os.path.join(OUT, "10_analysis_set.json")))
keep = cfg["keep"]
pid = np.load(os.path.join(OUT, "06_pid_gap.npy"))
accs0 = json.load(open(os.path.join(OUT, "06_accession_order.json")))
i0 = {a: i for i, a in enumerate(accs0)}
region_div = list(csv.DictReader(open(os.path.join(OUT, "11_region_diversity.csv"), encoding="utf-8-sig")))
win = list(csv.DictReader(open(os.path.join(OUT, "11_sliding_window_pi.csv"), encoding="utf-8-sig")))

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

# alignment column -> NC_011058 genome coordinate
col2pos, p = {}, 0
for c, ch in enumerate(aln["NC_011058"]):
    col2pos[c] = p
    if ch != "-":
        p += 1
GENOME_LEN = p

# =========================================================== FIGURE 1
print("Figure 1: sampling by geography and time")
fig = plt.figure(figsize=(7.2, 5.6))
gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.05], height_ratios=[1, 0.92],
                      hspace=0.52, wspace=0.58)

# 1a country counts, split complete vs partial
ax = fig.add_subplot(gs[0, :])
cc = Counter(meta[a]["country"] for a in meta)
order = [c for c, _ in cc.most_common()]
comp = [sum(1 for a in meta if meta[a]["country"] == c and meta[a]["molecule_class"] == "full genome") for c in order]
part = [cc[c] - n for c, n in zip(order, comp)]
x = np.arange(len(order))
ax.bar(x, comp, 0.74, color=CAT[0], label="complete genome", zorder=3)
ax.bar(x, part, 0.74, bottom=comp, color=CAT[1], label="partial sequence", zorder=3)
for xi, (c, p_) in enumerate(zip(comp, part)):
    ax.text(xi, c + p_ + 3, str(c + p_), ha="center", va="bottom",
            fontsize=7, color=INK2)
ax.set_xticks(x); ax.set_xticklabels(order, rotation=38, ha="right")
ax.set_ylabel("sequences in GenBank set")
ax.set_title("a   CpCDV sequences by country of origin", loc="left", color=INK)
ax.legend(loc="upper right", ncol=1)
ax.set_ylim(0, max(cc.values()) * 1.16)
grid_y(ax)

# 1b cumulative accumulation by subregion
ax = fig.add_subplot(gs[1, 0])
sub_years = defaultdict(list)
for a in meta:
    y = meta[a]["collection_year"]
    if y.isdigit():
        sub_years[meta[a]["subregion"]].append(int(y))
subs = [s for s, _ in Counter(meta[a]["subregion"] for a in meta).most_common()]
yrs = np.arange(1993, 2025)
for i, s in enumerate(subs[:6]):
    v = np.array([sum(1 for y in sub_years[s] if y <= Y) for Y in yrs])
    ax.plot(yrs, v, lw=2, color=CAT[i], label=s, solid_capstyle="round", zorder=3)
ax.set_xlabel("collection year"); ax.set_ylabel("cumulative sequences")
ax.set_title("b   accumulation of records", loc="left", color=INK)
ax.set_xlim(1993, 2024)
ax.legend(loc="upper left", fontsize=6.8, handlelength=1.2,
          labelspacing=0.32, handletextpad=0.5, borderpad=0.2)
grid_y(ax)

# 1c host family
ax = fig.add_subplot(gs[1, 1])
hf = Counter(meta[a]["host_family"] for a in meta if meta[a]["host_family"])
hfo = [h for h, _ in hf.most_common()]
vals = [hf[h] for h in hfo]
ax.barh(np.arange(len(hfo))[::-1], vals, 0.72, color=CAT[2], zorder=3)
for i, v in enumerate(vals):
    ax.text(v + 3, len(hfo) - 1 - i, str(v), va="center", fontsize=7, color=INK2)
ax.set_yticks(np.arange(len(hfo))[::-1]); ax.set_yticklabels(hfo)
ax.set_xlabel("sequences"); ax.set_xlim(0, max(vals) * 1.17)
ax.set_title("c   host family", loc="left", color=INK)
grid_x(ax)
save(fig, "Fig_S1_sampling_geography_time_host")

# =========================================================== FIGURE 2
print("Figure 2: pairwise identity and demarcation thresholds")
fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.9),
                         gridspec_kw=dict(wspace=0.3, width_ratios=[1.15, 1]))
ki = [i0[a] for a in keep]
sub = pid[np.ix_(ki, ki)] * 100
m = ~np.eye(len(ki), dtype=bool)
vals = sub[m]

ax = axes[0]
ax.hist(vals, bins=np.arange(70, 100.5, 0.5), color=CAT[0], zorder=3)
for thr, lab in ((78, "78%  species"), (94, "94%  strain")):
    ax.axvline(thr, color=CAT[7], lw=1.3, ls=(0, (4, 2)), zorder=4)
    ax.text(thr, ax.get_ylim()[1] * 0.97, " " + lab, rotation=90,
            va="top", ha="left", fontsize=7, color=CAT[7])
ax.set_xlabel("pairwise nucleotide identity (%)")
ax.set_ylabel("number of genome pairs")
ax.set_title(f"a   all {len(ki)} complete genomes ({len(vals)//2:,} pairs)",
             loc="left", color=INK)
grid_y(ax)

ax = axes[1]
# Strain labels as published (GenBank strain qualifier / isolate designation).
# The identity-based reassignment is NOT used here: it is defined by the 94%
# rule, so testing it against that rule would be circular.
within, between = [], []
for x in range(len(ki)):
    for y in range(x + 1, len(ki)):
        sx, sy = meta[keep[x]]["strain"], meta[keep[y]]["strain"]
        if not (len(sx) == 1 and sx.isalpha() and len(sy) == 1 and sy.isalpha()):
            continue
        (within if sx == sy else between).append(sub[x, y])
n_breach = sum(1 for v in within if v < 94)
n_merge = sum(1 for v in between if v > 94)
print(f"   annotated strains: {n_breach} within-strain pairs <94%, "
      f"{n_merge} between-strain pairs >94%")
bp = ax.boxplot([within, between], widths=0.5, patch_artist=True,
                showfliers=True, medianprops=dict(color=SURFACE, lw=1.4),
                whiskerprops=dict(color=INK3, lw=0.7),
                capprops=dict(color=INK3, lw=0.7),
                flierprops=dict(marker="o", markersize=2.2, markerfacecolor=INK3,
                                markeredgecolor="none", alpha=0.45),
                boxprops=dict(lw=0))
for patch, c in zip(bp["boxes"], (CAT[2], CAT[1])):
    patch.set_facecolor(c)
ax.axhline(94, color=CAT[7], lw=1.3, ls=(0, (4, 2)), zorder=4)
ax.text(2.46, 94, "94%", fontsize=7, color=CAT[7], va="center")
ax.annotate(f"{n_breach} pairs\nbelow 94%", xy=(1.17, 90), fontsize=6.8,
            color=CAT[7], va="center", ha="left")
ax.annotate(f"{n_merge} pairs\nabove 94%", xy=(2.2, 97.5), fontsize=6.8,
            color=CAT[7], va="center", ha="left")
ax.set_xticks([1, 2])
ax.set_xticklabels([f"within strain\n(n={len(within):,})",
                    f"between strains\n(n={len(between):,})"])
ax.set_ylabel("pairwise identity (%)")
ax.set_title("b   strains as published in GenBank", loc="left", color=INK)
ax.set_ylim(70, 101); ax.set_xlim(0.5, 3.1)
grid_y(ax)
save(fig, "Fig_S2_pairwise_identity_thresholds")

# =========================================================== FIGURE 3
print("Figure 3: diversity along the genome")
fig = plt.figure(figsize=(7.2, 4.9))
gs = fig.add_gridspec(3, 1, height_ratios=[0.32, 1.0, 1.05], hspace=0.75)

GENES = [("V1\nMP", 131, 410, CAT[0]), ("V2  capsid protein (CP)", 422, 1160, CAT[1]),
         ("C1:C2  Rep / RepA", 1317, 2408, CAT[2])]
ax = fig.add_subplot(gs[0])
ax.add_patch(Rectangle((0, 0.32), GENOME_LEN, 0.12, color="#d9d8d3", lw=0))
for name, s, e, c in GENES:
    ax.add_patch(Rectangle((s, 0.2), e - s, 0.36, color=c, lw=0))
    ax.text((s + e) / 2, 0.38, name, ha="center", va="center",
            fontsize=6.8, color="white", linespacing=0.9)
ax.text(0, 0.78, "LIR", ha="left", fontsize=6.8, color=INK2)
ax.text(1238, 0.78, "SIR", ha="center", fontsize=6.8, color=INK2)
ax.text(GENOME_LEN, 0.78, "LIR", ha="right", fontsize=6.8, color=INK2)
ax.set_xlim(0, GENOME_LEN); ax.set_ylim(0, 1)
ax.axis("off")
ax.set_title("a   genome organization (coordinates of NC_011058, "
             "origin at the Rep nick)", loc="left", color=INK, pad=2)

ax = fig.add_subplot(gs[1])
wx = [col2pos[int(r["window_mid"])] for r in win]
wy = [float(r["pi"]) for r in win]
ax.plot(wx, wy, lw=2, color=CAT[0], solid_capstyle="round", zorder=3)
ax.set_xlim(0, GENOME_LEN)
ax.set_ylabel("nucleotide diversity $\\pi$")
ax.set_xlabel("genome position (nt)")
ax.set_title("b   sliding-window nucleotide diversity (200 nt window, 25 nt step)",
             loc="left", color=INK)
for name, s, e, c in GENES:
    ax.axvspan(s, e, color=c, alpha=0.07, lw=0, zorder=1)
grid_y(ax)

ax = fig.add_subplot(gs[2])
want = ["LIR (large intergenic region)", "SIR (small intergenic region)",
        "V1 / movement protein", "V2 / capsid protein",
        "C1:C2 / Rep (spliced)", "C1 / RepA"]
lab = {"LIR (large intergenic region)": "LIR", "SIR (small intergenic region)": "SIR",
       "V1 / movement protein": "V1 (MP)", "V2 / capsid protein": "V2 (CP)",
       "C1:C2 / Rep (spliced)": "Rep (C1:C2)", "C1 / RepA": "RepA (C1)"}
rd = {r["region"]: r for r in region_div}
sel = [r for r in want if r in rd]
pis = [float(rd[r]["pi"]) for r in sel]
ordr = np.argsort(pis)[::-1]
xs = np.arange(len(sel))
cols = [CAT[3] if "LIR" in sel[i] or "SIR" in sel[i] else CAT[0] for i in ordr]
ax.bar(xs, [pis[i] for i in ordr], 0.68, color=cols, zorder=3)
for xi, i in enumerate(ordr):
    ax.text(xi, pis[i] + 0.004, f"{pis[i]:.3f}", ha="center", fontsize=7, color=INK2)
    ax.text(xi, 0.004, f"S={rd[sel[i]]['S']}", ha="center", fontsize=6.3,
            color="white", va="bottom")
ax.set_xticks(xs); ax.set_xticklabels([lab[sel[i]] for i in ordr])
ax.set_ylabel("nucleotide diversity $\\pi$")
ax.set_ylim(0, max(pis) * 1.2)
ax.set_title("c   diversity by genome region (yellow = non-coding)",
             loc="left", color=INK)
grid_y(ax)
save(fig, "Fig_S3_diversity_along_genome")

# =========================================================== FIGURE 4
print("Figure 4: strain distribution")
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.3),
                         gridspec_kw=dict(wspace=0.45, width_ratios=[1.25, 1]))
sc = Counter(sassign[a] for a in keep)
topS = [s for s, _ in sc.most_common() if s != "unassigned"][:7]
def sbin(s):
    return s if s in topS else ("unassigned" if s == "unassigned" else "other strain")
cmap = {s: CAT[i] for i, s in enumerate(topS)}
cmap["other strain"] = OTHER; cmap["unassigned"] = "#6b6a66"

ax = axes[0]
subs2 = [s for s, _ in Counter(meta[a]["subregion"] for a in keep).most_common()]
cats = topS + ["other strain", "unassigned"]
bottom = np.zeros(len(subs2))
for c in cats:
    v = np.array([sum(1 for a in keep if meta[a]["subregion"] == s and sbin(sassign[a]) == c)
                  for s in subs2], dtype=float)
    if v.sum() == 0:
        continue
    ax.bar(np.arange(len(subs2)), v, 0.74, bottom=bottom,
           color=cmap[c], label=c if c in ("other strain", "unassigned") else f"strain {c}",
           zorder=3, edgecolor=SURFACE, linewidth=1.1)
    bottom += v
for i, t in enumerate(bottom):
    ax.text(i, t + 2, f"{int(t)}", ha="center", fontsize=7, color=INK2)
ax.set_xticks(np.arange(len(subs2)))
ax.set_xticklabels([s.replace(" ", "\n") for s in subs2], fontsize=7)
ax.set_ylabel("complete genomes")
ax.set_ylim(0, bottom.max() * 1.16)
ax.set_title("a   strain composition by subregion", loc="left", color=INK)
ax.legend(ncol=2, loc="upper right", fontsize=6.8, handlelength=1.1,
          columnspacing=0.9, handletextpad=0.5)
grid_y(ax)

ax = axes[1]
hosts = [h for h, _ in Counter(meta[a]["host_family"] for a in keep).most_common() if h][:5]
M = np.zeros((len(hosts), len(subs2)))
for a in keep:
    hf, sr = meta[a]["host_family"], meta[a]["subregion"]
    if hf in hosts:
        M[hosts.index(hf), subs2.index(sr)] += 1
from matplotlib.colors import LinearSegmentedColormap, LogNorm
cm = LinearSegmentedColormap.from_list("seqblue", SEQ)
im = ax.imshow(np.where(M > 0, M, np.nan), cmap=cm, aspect="auto",
               norm=LogNorm(vmin=1, vmax=max(M.max(), 2)))
for i in range(len(hosts)):
    for j in range(len(subs2)):
        if M[i, j] > 0:
            ax.text(j, i, int(M[i, j]), ha="center", va="center", fontsize=6.8,
                    color="white" if M[i, j] > M.max() * 0.25 else INK)
ABBR = {"Arabian Peninsula": "Arabian\nPen."}
ax.set_xticks(range(len(subs2)))
ax.set_xticklabels([ABBR.get(s, s.replace(" ", "\n")) for s in subs2], fontsize=6.2)
ax.set_yticks(range(len(hosts))); ax.set_yticklabels(hosts, fontsize=7)
ax.set_title("b   host family by subregion", loc="left", color=INK)
for sp in ax.spines.values():
    sp.set_visible(False)
ax.tick_params(length=0)
cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
cb.set_label("genomes (log)", fontsize=7); cb.ax.tick_params(labelsize=6.5, length=0)
cb.outline.set_visible(False)
save(fig, "Fig_S4_strain_distribution")
print("\nfigures written to", FIG)
