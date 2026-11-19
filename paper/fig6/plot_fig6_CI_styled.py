"""Fig. 6 panels C-I in the style of panels A-B (Liberation Sans; bold panel letters 18 pt,
bold centred titles 10 pt, axis labels 10 pt), width = A4 text width of panels A-B (595.4 pt).
Usage: python3 plot_fig6_CI_styled.py DATA_DIR OUT.pdf"""
import glob, sys
from pathlib import Path
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
data, out = Path(sys.argv[1]), Path(sys.argv[2])
plt.rcParams.update({"font.family": "Liberation Sans", "mathtext.fontset": "custom",
    "mathtext.rm": "Liberation Sans", "mathtext.it": "Liberation Sans:italic", "mathtext.bf": "Liberation Sans:bold",
    "font.size": 8, "axes.labelsize": 10, "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.8, "pdf.fonttype": 42})
GRID = "#e4e3df"
COL = {"L1": "#2a78d6", "L2": "#eb6834", "L3": "#1baf7a", "background": "#52514e"}
DCOL = {20: "#86b6ef", 60: "#2a78d6", 150: "#0d366b"}
W = 595.446 / 72; H = 5.2
fig = plt.figure(figsize=(W, H))
v = pd.read_csv(data / "vaf_levels.csv.gz"); reps = v.rep.nunique()
bins = np.linspace(0, 1, 101)
def letter(ax, L, title):
    pos = ax.get_position()
    fig.text(pos.x0 - 0.045, pos.y1 + 0.075, L, fontsize=18, fontweight="bold", va="bottom", ha="left")
    ax.set_title(title, fontsize=10, fontweight="bold", pad=6)
left, right = 0.085, 0.985
top = [(left + i * 0.183, 0.135) for i in range(5)]
panels = [("1_layer_phased", "C", "Layer-specific,\nphased"), ("2_layer_unphased", "D", "Layer-specific,\nunphased"),
          ("3_bulk_unphased_exact", "E", "Bulk, unphased\n(exact)"), ("4_reads_d60", "F", "+ Read sampling\n(depth 60)"),
          ("5_reads_background_d60", "G", "+ Background\nartefacts")]
for i, (lev, L, title) in enumerate(panels):
    ax = fig.add_axes([left + i * 0.183, 0.625, 0.135, 0.245]); d = v[v.level == lev]
    kinds = ["L2"] if lev.startswith(("1", "2")) else ["L1", "L2", "L3"] + (["background"] if lev.startswith("5") else [])
    for k in kinds:
        x = d.vaf[d.kind == ("somatic" if lev.startswith(("1", "2")) else k)]
        h, _ = np.histogram(x, bins); ax.stairs(h / reps, bins, color=COL[k], lw=1.1, label=k.capitalize() if k == "background" else k)
    ax.set_yscale("log"); ax.set_ylim(0.5, 3e4); ax.set_xlim(0, 1); ax.set_xticks([0, 0.5, 1]); ax.set_xticklabels(["0", "0.5", "1"])
    ax.set_xlabel("VAF")
    if i == 0: ax.set_ylabel("Variants per tree")
    else: ax.set_yticklabels([])
    if lev.startswith("3"):
        for c in (0.05, 0.35, 0.10): ax.axvline(c, color=GRID, lw=0.8, zorder=0)
    ax.grid(axis="y", color=GRID, lw=0.5)
    if i in (2, 4): ax.legend(frameon=False, fontsize=7, handlelength=1.2, borderaxespad=0.1, loc="upper right")
    letter(ax, L, title)
ax = fig.add_axes([left, 0.085, 0.40, 0.30])
for dpt in (20, 60, 150):
    d = v[(v.level == f"4_reads_d{dpt}") & (v.kind == "L2")]
    h, _ = np.histogram(d.vaf, bins); ax.stairs(h / reps, bins, color=DCOL[dpt], lw=1.1, label=f"Mean depth {dpt}")
ax.axvline(0.35, color=GRID, lw=0.8, zorder=0); ax.set_yscale("log"); ax.set_ylim(0.5, 1e4); ax.set_xlim(0, 0.8)
ax.set_xlabel("Observed VAF"); ax.set_ylabel("Variants per tree"); ax.legend(frameon=False, loc="upper right"); ax.grid(axis="y", color=GRID, lw=0.5)
letter(ax, "H", "L2 variants by depth (L2-fixed: 0.35)")
ax = fig.add_axes([0.585, 0.085, 0.40, 0.30])
edges = np.geomspace(0.001, 0.5, 25)
for dpt in (20, 60, 150):
    d = pd.concat([pd.read_csv(f) for f in glob.glob(str(data / f"detection_rep*_d{dpt}_bg0.csv.gz"))])
    d = d[(d.true_vaf >= edges[0]) & (d.true_vaf < edges[-1])]
    idx = np.digitize(d.true_vaf, edges); g = d.groupby(idx).detected.mean()
    mid = np.sqrt(edges[g.index - 1] * edges[g.index])
    ax.plot(mid, g.values, color=DCOL[dpt], lw=1.3, marker="o", ms=2.5, label=f"Mean depth {dpt}")
ax.set_xscale("log"); ax.set_xlim(1e-3, 0.5); ax.set_ylim(0, 1.02)
ax.set_xlabel(r"True assay VAF ($a_k v/2$)"); ax.set_ylabel("Fraction detected")
ax.legend(frameon=False, loc="upper left"); ax.grid(color=GRID, lw=0.5)
letter(ax, "I", "Detection of true variants by depth")
fig.savefig(out); print("wrote", out)
