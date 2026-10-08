"""Schematic figures for the simSOMA Supplementary Text (revision).

FigureS1_selfrenewal.pdf   : one self-renewal round (division, displacement with independent
                            daughter, victim locality lambda_vic, favored clone w_B)
FigureS2_founding.pdf      : pre-branching amplification into sectors, index notation
                            (x_k, alpha_i, beta_i), branch founding and organ founding blocks,
                            monoclonal vs polyclonal founding
FigureS4_observation.pdf   : levels of the observation model
Usage: python3 si_schematics.py OUTDIR
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, FancyArrowPatch, Rectangle, Circle

OUT = sys.argv[1] if len(sys.argv) > 1 else "."
PAL = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3", "#937860"]
plt.rcParams.update({"font.size": 8, "font.family": "sans-serif"})


def ring(ax, cx, cy, colors, r=1.0, w=0.38, labels=None, hl=None, start=90, ec="white"):
    n = len(colors)
    pos = []
    for i, c in enumerate(colors):
        a0 = start - (i + 1) * 360 / n
        a1 = start - i * 360 / n
        lw = 2.2 if (hl and i in hl) else 0.8
        e = "black" if (hl and i in hl) else ec
        ax.add_patch(Wedge((cx, cy), r, a0, a1, width=w, facecolor=c, edgecolor=e, lw=lw))
        am = np.deg2rad((a0 + a1) / 2)
        pos.append((cx + (r - w / 2) * np.cos(am), cy + (r - w / 2) * np.sin(am)))
        if labels:
            ax.text(*pos[-1], labels[i], ha="center", va="center", fontsize=7, color="white", weight="bold")
    return pos


def arrow(ax, p, q, **kw):
    kw.setdefault("arrowstyle", "-|>")
    kw.setdefault("mutation_scale", 10)
    kw.setdefault("lw", 1.0)
    kw.setdefault("color", "black")
    ax.add_patch(FancyArrowPatch(p, q, **kw))


# ----------------------------------------------------------------------------- self-renewal
def fig_selfrenewal():
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.7))
    for ax in axes:
        ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-1.6, 1.6); ax.set_ylim(-1.75, 1.45)
    m = 6
    base = [PAL[i] for i in range(m)]

    ax = axes[0]
    ax.set_title("A  Niche state X(t-1)", loc="left", fontsize=9, weight="bold")
    ring(ax, 0, 0, base, labels=[str(i + 1) for i in range(m)], hl=[0])
    ax.text(0, 0, "m = 6", ha="center", va="center", fontsize=8)
    ax.text(0, -1.3, "ordered ring of m ASC positions;\nposition 1 = favored seed (q = 1, weight 1);\n"
            "others weight 1 - w_B", ha="center", va="top", fontsize=6.5)

    ax = axes[1]
    ax.set_title("B  Division (all m cells)", loc="left", fontsize=9, weight="bold")
    ring(ax, 0, 0, base, labels=[str(i + 1) for i in range(m)])
    for i in range(m):
        am = np.deg2rad(90 - (i + 0.5) * 360 / m)
        ax.plot(0.62 * np.cos(am), 0.62 * np.sin(am), "*", color="gold", ms=7, mec="k", mew=0.4)
    ax.text(0, 0, "+Pois(μ_div)\nper edge", ha="center", va="center", fontsize=6.5)
    ax.text(0, -1.3, "each cell divides once; the daughter\nkeeping the position receives\n"
            "K ~ Poisson(μ_div) new mutations", ha="center", va="top", fontsize=6.5)

    ax = axes[2]
    ax.set_title("C  Displacement (prob. ρ)", loc="left", fontsize=9, weight="bold")
    new = list(base)
    new[2] = base[1]
    pos = ring(ax, 0, 0, new, labels=["1", "2", "2'", "4", "5", "6"], hl=[2])
    arrow(ax, pos[1], pos[2], connectionstyle="arc3,rad=-0.5", color="black", lw=1.3)
    arrow(ax, pos[1], pos[5], connectionstyle="arc3,rad=0.35", color="grey", lw=0.9, linestyle="--")
    ax.text(0, -0.22, "λ_vic: neighbour\n1-λ_vic: any", ha="center", va="center", fontsize=6.5)
    ax.text(0, -1.3, "displacer 2 divides again; its second\ndaughter 2' (own mutation draw)\n"
            "replaces the victim", ha="center", va="top", fontsize=6.5)
    fig.tight_layout()
    fig.savefig(f"{OUT}/FigureS1_selfrenewal.pdf"); fig.savefig(f"{OUT}/FigureS1_selfrenewal.png", dpi=150)


# ---------------------------------------------------------------------------- founding
def strip(ax, y, sectors, C, x0=0.0, h=0.45, width=10.0, hl=None, hlcolor="black", label=None):
    w = width / C
    k = 0
    for i, s in enumerate(sectors):
        for j in range(s):
            ax.add_patch(Rectangle((x0 + k * w, y), w, h, facecolor=PAL[i % len(PAL)],
                                   edgecolor="white", lw=0.3))
            k += 1
    if hl is not None:
        a, b = hl
        ax.add_patch(Rectangle((x0 + a * w, y - 0.06), (b - a) * w, h + 0.12, fill=False,
                               edgecolor=hlcolor, lw=1.8))
    if label:
        ax.text(x0 - 0.2, y + h / 2, label, ha="right", va="center", fontsize=7)
    return w


def fig_founding():
    fig = plt.figure(figsize=(7.2, 4.6))
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0]); ax.axis("off")
    ax.set_xlim(-2.6, 12.4); ax.set_ylim(-0.1, 9.2)

    m, C = 4, 24
    sec = [C // m] * m
    # A: niche -> sectors
    ax.text(-2.5, 8.9, "A  Pre-branching amplification", fontsize=9, weight="bold", va="top")
    for i in range(m):
        ax.add_patch(Rectangle((3.0 + i * 1.0, 7.6), 0.9, 0.45, facecolor=PAL[i], edgecolor="white"))
        ax.text(3.45 + i * 1.0, 7.82, f"c{i + 1}", ha="center", va="center", color="white", fontsize=7, weight="bold")
    ax.text(2.8, 7.82, "niche X (m = 4)", ha="right", va="center", fontsize=7)
    for i in range(m):
        arrow(ax, (3.45 + i * 1.0, 7.55), (1.25 + 2.5 * i, 6.75), mutation_scale=7, lw=0.7, color="grey")
    w = strip(ax, 6.2, sec, C, label="SAM-boundary\npopulation R (C = 24)")
    ax.text(10.3, 7.15, "balanced binary\ngrowth; mutations\non every division\nedge", fontsize=6.5, va="top")
    # indices
    for i in range(m):
        a = i * sec[0]
        ax.plot([a * w, a * w], [6.05, 6.75], color="k", lw=0.6)
        ax.text(a * w, 5.85, f"α{i + 1}={a}", fontsize=6, ha="center", va="top")
        ax.text((a + sec[0] / 2) * w, 6.42, f"S{i + 1}", ha="center", va="center", fontsize=7, color="white", weight="bold")
    ax.plot([C * w, C * w], [6.05, 6.75], color="k", lw=0.6)
    ax.text(C * w, 5.85, f"β{m}={C}", fontsize=6, ha="center", va="top")
    ax.text(0.0, 5.35, "cells x0, x1, ..., x(C-1); sector S_i occupies indices [α_i, β_i), "
            "σ_i = β_i - α_i ≈ C/m (ring: x(C-1) is adjacent to x0)", fontsize=6.5, va="top")

    # B: branch founding
    ax.text(-2.5, 4.75, "B  Branch founding (P_b,eff = min(P_b, m, C))", fontsize=9, weight="bold", va="top")
    strip(ax, 3.6, sec, C, hl=(1, 4), label="monoclonal\n(P_b,eff = 3)")
    strip(ax, 2.7, sec, C, hl=(11, 14), hlcolor="red", label="polyclonal\n(P_b,eff = 3)")
    # child niche
    ax.text(10.3, 4.05, "child niche\n(m = 4):", fontsize=6.5, va="center")
    for j, col in enumerate([0, 0, 0, 0]):
        ax.add_patch(Rectangle((10.3 + j * 0.2, 3.55), 0.18, 0.3, facecolor=PAL[col], edgecolor="white", lw=0.3))
    for j, col in enumerate([1, 1, 2, 2]):
        ax.add_patch(Rectangle((10.3 + j * 0.2, 2.75), 0.18, 0.3, facecolor=PAL[col], edgecolor="white", lw=0.3))

    # C: organ founding
    ax.text(-2.5, 2.25, "C  Organ founding (P_o precursors, grown to O cells)", fontsize=9, weight="bold", va="top")
    strip(ax, 1.1, sec, C, hl=(13, 21), hlcolor="red", label="P_o = 8\n(polyclonal)")
    ax.text(10.3, 1.33, "organ:\nO cells,\nn_seq sampled", fontsize=6.5, va="center")
    ax.text(0.0, 0.75,
            "A block of P contiguous cells is polyclonal if it spans a sector boundary. For uniformly placed blocks,\n"
            "Pr(polyclonal) = min{1, φ} with φ = (P - 1) m / C (exact if C is divisible by m and P ≤ C/m + 1);\n"
            "the expected number of sectors in the block is 1 + Pr(polyclonal).",
            fontsize=6.5, va="top")
    fig.savefig(f"{OUT}/FigureS2_founding.pdf"); fig.savefig(f"{OUT}/FigureS2_founding.png", dpi=150)


# --------------------------------------------------------------------------- observation
def fig_observation():
    fig = plt.figure(figsize=(7.2, 1.9))
    ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off"); ax.set_xlim(-0.1, 15.1); ax.set_ylim(0.9, 4.8)
    boxes = [
        ("Developmental VAF\n$v^{(k)}_{io}$", "per layer k\n(independent lineage\nhistories, one topology)"),
        ("Assay VAF\n$\\tilde v_{io}=\\eta\\sum_k a_{ko} v^{(k)}_{io}$", "layer contributions $a_{ko}$\nphasing $\\eta\\in\\{1,1/2\\}$"),
        ("Read depth\n$D_{io}$", "fixed | normal | Poisson |\nlog-normal site, site-by-organ\nand organ factors (apricot)"),
        ("Variant reads\n$A_{io}\\sim\\mathrm{Bin}(D_{io},p_{io})$", "$p=\\tilde v(1-e)+(1-\\tilde v)e$\n(or beta-binomial)\n+ background sites"),
        ("Calls / ascertainment", "min depth, min alt reads,\nmin VAF; keep sites called\nin >= 1 sample"),
    ]
    xs = np.linspace(1.5, 13.5, len(boxes))
    for x, (t, s) in zip(xs, boxes):
        ax.add_patch(Rectangle((x - 1.3, 2.4), 2.6, 1.6, facecolor="#E8EEF7", edgecolor="#2F4B7C", lw=1))
        ax.text(x, 3.2, t, ha="center", va="center", fontsize=7)
        ax.text(x, 2.2, s, ha="center", va="top", fontsize=6.3, color="#333333")
    for a, b in zip(xs[:-1], xs[1:]):
        arrow(ax, (a + 1.32, 3.2), (b - 1.32, 3.2), mutation_scale=9)
    ax.text(xs[0], 4.45, "simSOMA (developmental)", ha="center", fontsize=7, weight="bold")
    ax.text((xs[1] + xs[4]) / 2, 4.45, "observation model (plantsoma_obs, versioned)", ha="center", fontsize=7, weight="bold")
    ax.plot([xs[1] - 1.3, xs[4] + 1.3], [4.3, 4.3], color="#2F4B7C", lw=0.8)
    fig.savefig(f"{OUT}/FigureS4_observation.pdf"); fig.savefig(f"{OUT}/FigureS4_observation.png", dpi=150)


if __name__ == "__main__":
    fig_selfrenewal(); fig_founding(); fig_observation()
    print("ok")
