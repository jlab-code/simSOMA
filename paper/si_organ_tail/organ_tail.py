"""SI figure: neutral 1/q tail of organ-private variants and the organ mutation-rate multiplier.

Theory (balanced binary expansion from P_o precursors to O cells, rate mu_div per daughter edge,
multiplier g for organ divisions): number of organ-private variants with cellular frequency >= q,
for 1/O <= q <= 1/P_o,  M(>= q) ~= 2 g mu_div (1/q - P_o).
Simulated with the full pipeline (SI topology 01). Usage: python3 organ_tail.py OUT.png
"""
import sys, collections
from pathlib import Path
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "simSOMA_corefunc"))
import topology_io, pipeline_wrapper, self_renewal  # noqa: E402
out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("organ_tail.png")
topo = topology_io.load_topology_auto(ROOT / "simSOMA_inputs/topologies/topology_01_multibranch_20organs_order1to5.json",
                                      mapping={"unit": "years", "rate": 5.0, "mode": "deterministic"})
O, mu_unit, kappa = 4096, 1.0, 5.0
mu_div = mu_unit / kappa
COL = {(4, 1.0): "#86b6ef", (16, 1.0): "#2a78d6", (4, 3.0): "#f0a07a", (16, 3.0): "#eb6834"}
fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
for (Po, g), col in COL.items():
    curves = []
    for s in range(5):
        sr = self_renewal.SelfRenewalParams(m=4, rho=0.0, mu_unit=mu_unit, kappa_sr=kappa)
        r = pipeline_wrapper.run_pipeline(topo, sr, sam_boundary_cells=64, branch_precursor_number=4,
                                          organ_precursor_number=Po, organ_total_cells=O, sequenced_cells=O,
                                          seed=100 + s, organ_mu_multiplier=g)
        deg = collections.Counter(mid for oe in r["organ_events"] for mid in oe["allele_counts_by_mutation"])
        for oe in r["organ_events"]:
            q = np.array([a / O for mid, a in oe["allele_counts_by_mutation"].items() if deg[mid] == 1])
            curves.append(q)
    qs = 2.0 ** np.arange(0, int(np.log2(O / Po)) + 1) / O   # dyadic frequencies (balanced growth, O/P_o a power of 2)
    M = np.array([[np.sum(q >= x - 1e-12) - np.sum(q >= 1 / Po - 1e-12) for x in qs] for q in curves]).mean(0)
    th = 2 * g * mu_div * (1 / qs - Po)
    lab = f"P_o={Po}, organ multiplier {g:g}"
    ax[0].plot(1 / qs, M, "o", ms=4, color=col, label=lab); ax[0].plot(1 / qs, th, "-", lw=1, color=col)
    ax[1].plot(qs, M, "o", ms=3, color=col); ax[1].plot(qs, th, "-", lw=1, color=col)
ax[0].set_xlabel("1 / q"); ax[0].set_ylabel("organ-private variants with frequency >= q")
ax[0].set_title("A  Points: simulation; lines: 2 g $\\mu_{div}$ (1/q - $P_o$)", loc="left", fontsize=8)
ax[0].legend(frameon=False, fontsize=7)
ax[1].set_xscale("log"); ax[1].set_yscale("log"); ax[1].set_xlabel("cellular frequency q")
ax[1].set_title("B  Same data, log-log", loc="left", fontsize=8)
for a in ax:
    a.spines[["top", "right"]].set_visible(False); a.grid(color="#e4e3df", lw=0.5)
fig.tight_layout(); fig.savefig(out, dpi=200); print("wrote", out)
