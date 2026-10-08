"""Data for the revised Fig. 6: from developmental to observed VAF spectra.

One set of simulated trees (SI topology 01, 20 organs), three SAM layers simulated independently
(layered.py), then successive observation levels:
  level 1  layer-specific, phased      : v = f                      (developmental VAF, L2)
  level 2  layer-specific, unphased    : v = f / 2
  level 3  bulk, unphased (exact)      : v = c_k f / 2              (c = 0.10, 0.70, 0.20)
  level 4  + reads                     : depth ~ log-normal site x organ x organ-factor model (plantsoma_obs 1.1.0, mean 60), binomial reads,
                                         called if >= 2 alt reads in >= 1 organ
  level 5  + background artefacts      : 2000 sites, VAF ~ gamma(shape 2, mean 0.023), same in all organs
Depth series: level 4 at mean depth 20, 60, 150.

Usage: python3 paper/fig6/make_fig6_data.py [--reps 10] [--out paper/fig6/data]
"""
import argparse, json, sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "simSOMA_corefunc"))
import layered, topology_io, plantsoma_obs  # noqa: E402,E401

ap = argparse.ArgumentParser()
ap.add_argument("--reps", type=int, default=10)
ap.add_argument("--out", type=Path, default=ROOT / "paper" / "fig6" / "data")
ap.add_argument("--topology", type=Path, default=ROOT / "simSOMA_inputs" / "topologies" / "topology_01_multibranch_20organs_order1to5.json")
a = ap.parse_args()
a.out.mkdir(parents=True, exist_ok=True)

SETTINGS = dict(
    mapping={"unit": "years", "rate": 5.0, "mode": "deterministic"},
    params=dict(m=4, rho=0.0, kappa_sr=5.0, sam_boundary_cells=64, branch_precursor_number=4,
                organ_precursor_number=8, organ_total_cells=10000, sequenced_cells=1000),
    layers={"L1": {"mu_unit": 1.0}, "L2": {"mu_unit": 1.0}, "L3": {"mu_unit": 1.0}},
    contributions={"L1": 0.10, "L2": 0.70, "L3": 0.20},
    reads={"depth": {"mode": "lognormal_site_sample", "mean": 60},
           "reads": {"type": "binomial", "sequencing_error": 0.0},
           "caller": {"min_depth": 0, "min_alt_reads": 2, "retain_called_any": True}},
    background={"n_sites": 2000, "distribution": "gamma:2", "mean": 0.023},
    depth_series=[20, 60, 150],
    seed=20261008,
)

topo = topology_io.load_topology_auto(a.topology, mapping=SETTINGS["mapping"])
rows = []
for rep in range(a.reps):
    car = layered.simulate_layers(topo, SETTINGS["layers"], seed=SETTINGS["seed"], replicate=rep, **SETTINGS["params"])
    organs = sorted(car.organ_id.unique())
    l2 = car[car.layer == "L2"]
    rows += [dict(rep=rep, level="1_layer_phased", vaf=v, called=True, kind="somatic") for v in l2.carrier_fraction]
    rows += [dict(rep=rep, level="2_layer_unphased", vaf=v / 2, called=True, kind="somatic") for v in l2.carrier_fraction]
    cm = layered.resolve_contributions(SETTINGS["contributions"], organs, list(SETTINGS["layers"]))
    sites, organs, V, F = layered.assay_vaf_matrix(car, cm, "unphased")
    nz = F > 0
    lay = np.repeat(sites.layer.to_numpy()[:, None], len(organs), axis=1)
    rows += [dict(rep=rep, level="3_bulk_unphased_exact", vaf=v, called=True, kind=k) for v, k in zip(V[nz], lay[nz])]
    for depth in SETTINGS["depth_series"]:
        for with_bg in ([False, True] if depth == 60 else [False]):
            cfg = json.loads(json.dumps(SETTINGS["reads"])); cfg["depth"]["mean"] = depth
            if with_bg:
                cfg["background"] = SETTINGS["background"]
            obs = plantsoma_obs.observe(V, cfg, np.random.default_rng(SETTINGS["seed"] + 1000 * rep + depth + 7 * with_bg))
            n = len(V)
            obs_vaf, keep, bg = obs["observed_vaf"], obs["keep"], obs["is_background"]
            present = np.vstack([nz, np.ones((bg.sum(), len(organs)), bool)]) if bg.any() else nz
            # observed VAF per (site, organ) where the site is ascertained and alt > 0 in that organ
            sel = keep[:, None] & (obs["alt"] > 0)
            kind = np.vstack([lay, np.full((bg.sum(), len(organs)), "background")]) if bg.any() else lay
            level = f"5_reads_background_d{depth}" if with_bg else f"4_reads_d{depth}"
            rows += [dict(rep=rep, level=level, vaf=v, called=True, kind=k) for v, k in zip(obs_vaf[sel], kind[sel])]
            # detection of true (site, organ) presences, by true assay VAF
            tv = V[nz]; det = (keep[:n, None] & (obs["alt"][:n] > 0))[nz]
            pd.DataFrame({"rep": rep, "depth": depth, "background": with_bg, "true_vaf": tv, "detected": det,
                          "layer": lay[nz]}).to_csv(a.out / f"detection_rep{rep}_d{depth}_bg{int(with_bg)}.csv.gz", index=False)
    print(f"replicate {rep}: {len(car)} carrier rows, {len(V)} sites")
pd.DataFrame(rows).to_csv(a.out / "vaf_levels.csv.gz", index=False)
(a.out / "settings.json").write_text(json.dumps({**SETTINGS, "reps": a.reps, "topology": str(a.topology),
                                                 "plantsoma_obs_version": plantsoma_obs.OBSERVATION_MODEL_VERSION,
                                                 "reads_config_sha256": plantsoma_obs.config_sha256(SETTINGS["reads"])}, indent=2))
print("wrote", a.out)
