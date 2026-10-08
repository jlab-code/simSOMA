"""Data panels C-I of the revised Fig. 6 (panels A-B: original schematic).
Usage: python3 paper/fig6/plot_fig6_draft.py DATA_DIR OUT.png"""
import glob, sys
from pathlib import Path
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt

data, out = Path(sys.argv[1]), Path(sys.argv[2])
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
COL = {"L1": "#2a78d6", "L2": "#eb6834", "L3": "#1baf7a", "background": "#52514e"}
DCOL = {20: "#86b6ef", 60: "#2a78d6", 150: "#0d366b"}
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.labelcolor": INK})
v = pd.read_csv(data / "vaf_levels.csv.gz")
reps = v.rep.nunique()
panels = [("1_layer_phased", "C  Layer-specific, phased"), ("2_layer_unphased", "D  Layer-specific, unphased"),
          ("3_bulk_unphased_exact", "E  Bulk, unphased (exact)"), ("4_reads_d60", "F  + read sampling (depth 60)"),
          ("5_reads_background_d60", "G  + background artefacts")]
fig = plt.figure(figsize=(11, 5.6))
gs = fig.add_gridspec(2, 5, height_ratios=[1, 1], hspace=0.55, wspace=0.35)
bins = np.linspace(0, 1, 101)
for i, (lev, title) in enumerate(panels):
    ax = fig.add_subplot(gs[0, i]); d = v[v.level == lev]
    kinds = ["L2"] if lev.startswith(("1", "2")) else ["L1", "L2", "L3"] + (["background"] if lev.startswith("5") else [])
    for k in kinds:
        x = d.vaf[d.kind == ("somatic" if lev.startswith(("1", "2")) else k)]
        h, _ = np.histogram(x, bins); ax.stairs(h / reps, bins, color=COL[k], lw=1.4, label=k)
    ax.set_yscale("log"); ax.set_ylim(0.5, None); ax.set_xlim(0, 1)
    ax.set_title(title, loc="left", fontsize=8.5, color=INK); ax.set_xlabel("VAF")
    if i == 0: ax.set_ylabel("variants per tree (organs pooled)")
    if lev.startswith("3"):
        for c in (0.05, 0.35, 0.10): ax.axvline(c, color=GRID, lw=0.8, zorder=0)
    ax.grid(axis="y", color=GRID, lw=0.5)
    if i in (2, 4): ax.legend(frameon=False, fontsize=7)
# depth series
ax = fig.add_subplot(gs[1, 0:2])
for dpt in (20, 60, 150):
    d = v[(v.level == f"4_reads_d{dpt}") & (v.kind == "L2")]
    h, _ = np.histogram(d.vaf, bins); ax.stairs(h / reps, bins, color=DCOL[dpt], lw=1.4, label=f"mean depth {dpt}")
ax.axvline(0.35, color=GRID, lw=0.8, zorder=0); ax.set_yscale("log"); ax.set_ylim(0.5, None); ax.set_xlim(0, 0.8)
ax.set_title("H  L2 variants by depth (L2-fixed: 0.35)", loc="left", fontsize=8.5)
ax.set_xlabel("observed VAF"); ax.set_ylabel("variants per tree"); ax.legend(frameon=False, fontsize=7); ax.grid(axis="y", color=GRID, lw=0.5)
ax = fig.add_subplot(gs[1, 2:4])
edges = np.geomspace(0.001, 0.5, 25)
for dpt in (20, 60, 150):
    d = pd.concat([pd.read_csv(f) for f in glob.glob(str(data / f"detection_rep*_d{dpt}_bg0.csv.gz"))])
    d = d[(d.true_vaf >= edges[0]) & (d.true_vaf < edges[-1])]
    idx = np.digitize(d.true_vaf, edges); g = d.groupby(idx).detected.mean()
    mid = np.sqrt(edges[g.index - 1] * edges[g.index])
    ax.plot(mid, g.values, color=DCOL[dpt], lw=1.6, marker="o", ms=3, label=f"mean depth {dpt}")
ax.set_xscale("log"); ax.set_xlim(1e-3, 0.5); ax.set_ylim(0, 1.02); ax.set_xlabel("true assay VAF (c_k f / 2)"); ax.set_ylabel("fraction detected")
ax.set_title("I  Detection of true variants by depth", loc="left", fontsize=8.5); ax.legend(frameon=False, fontsize=7); ax.grid(color=GRID, lw=0.5)
fig.savefig(out, dpi=200, bbox_inches="tight"); print("wrote", out)
