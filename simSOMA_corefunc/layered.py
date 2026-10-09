"""Per-layer simulation and read-level observation of layered organs.

Model
-----
The SAM layers (e.g. L1, L2, L3) are simulated as independent cell-lineage histories on the
same plant topology: one simSOMA run per layer, each with its own random seed and its own
mutation rate (mu_unit), and shared developmental parameters (m, rho, C, P_b, P_o, O, ...).
The topology, including any phyllotactic event angles, is loaded once and shared, so branch
and organ initiation sites occupy the same angular positions in every layer. There is no
cell exchange between layers, so a mutation arises in, and can become fixed in, one layer only.

For a mutation i of layer k with carrier-cell fraction f_iko in organ o (fraction of sampled
cells of that layer carrying it), the assay VAF of a bulk sample of organ o is

    v_io = eta * c_ko * f_iko,     eta = 1/2 (unphased diploid) or 1 (phased),

with layer contributions c_ko (sum_k c_ko = 1; one global vector or one per organ). Layer-
specific sampling of layer k* is the special case c_ko = 1(k = k*). Read counts follow from
the shared plantsoma_obs model.

This replaces the post-hoc approximation of simSOMA <= 0.1.x, which copied one layer-equivalent
spectrum into every layer (identical lineage histories in all layers).
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence

import numpy as np
import pandas as pd

import pipeline_wrapper
import self_renewal
import topology_io

try:
    import plantsoma_obs
except ImportError:  # repository checkout without installation
    import sys as _sys
    _sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    import plantsoma_obs

LAYERED_API_VERSION = "1.0.0"
PHASE_FACTORS = {"phased": 1.0, "unphased": 0.5}


def layer_seed(seed: int, replicate: int, layer: str) -> int:
    h = hashlib.sha256(f"simSOMA-layer|{int(seed)}|{int(replicate)}|{layer}".encode()).hexdigest()
    return int(h[:15], 16)


def simulate_layers(
    topo: Dict[str, Any],
    layers: Mapping[str, Mapping[str, Any]],
    *,
    m: int,
    rho: float,
    kappa_sr: float,
    sam_boundary_cells: int,
    branch_precursor_number: int,
    organ_precursor_number: int,
    organ_total_cells: int,
    sequenced_cells: Optional[int] = None,
    organ_mu_multiplier: float = 1.0,
    victim_locality: float = 0.0,
    seed: int = 0,
    replicate: int = 0,
) -> pd.DataFrame:
    """One independent simSOMA run per layer. Returns one row per (mutation, layer, organ)
    with carrier_fraction = carrier cells / sampled cells of that layer."""
    if not layers:
        raise ValueError("at least one layer is required")
    rows = []
    for name, spec in layers.items():
        if "mu_unit" not in spec:
            raise ValueError(f"layer {name!r}: mu_unit is required")
        sr = self_renewal.SelfRenewalParams(m=int(m), rho=float(rho), mu_unit=float(spec["mu_unit"]),
                                            kappa_sr=float(kappa_sr), victim_locality=float(victim_locality))
        res = pipeline_wrapper.run_pipeline(
            topo, sr, sam_boundary_cells=int(sam_boundary_cells),
            branch_precursor_number=int(branch_precursor_number),
            organ_precursor_number=int(organ_precursor_number), organ_total_cells=int(organ_total_cells),
            sequenced_cells=sequenced_cells, seed=layer_seed(seed, replicate, name),
            organ_mu_multiplier=float(spec.get("organ_mu_multiplier", organ_mu_multiplier)))
        for oe in res["organ_events"]:
            n = int(oe["sequenced_cells"])
            for mid, a in oe["allele_counts_by_mutation"].items():
                rows.append((f"{name}:{mid}", name, str(oe["organ_id"]), int(a), n, int(a) / n))
    return pd.DataFrame(rows, columns=["mutation_id", "layer", "organ_id", "carrier_cells", "sampled_cells",
                                       "carrier_fraction"])


def resolve_contributions(contributions: Any, organs: Sequence[str], layers: Sequence[str]) -> pd.DataFrame:
    """Organ x layer contribution matrix from {layer: c} (all organs) or {organ: {layer: c}}."""
    if not isinstance(contributions, Mapping) or not contributions:
        raise ValueError("layer_contributions must be {layer: c} or {organ: {layer: c}}")
    per_organ = all(isinstance(v, Mapping) for v in contributions.values())
    rows = {}
    for o in organs:
        spec = contributions.get(o, contributions.get("default")) if per_organ else contributions
        if spec is None:
            raise ValueError(f"no layer contributions for organ {o!r} (and no 'default')")
        unknown = set(spec) - set(layers)
        if unknown:
            raise ValueError(f"layer contributions name unknown layers {sorted(unknown)}; simulated: {list(layers)}")
        vec = [float(spec.get(k, 0.0)) for k in layers]
        if min(vec) < 0 or abs(sum(vec) - 1.0) > 1e-6:
            raise ValueError(f"layer contributions for organ {o!r} must be >= 0 and sum to 1, got {dict(zip(layers, vec))}")
        rows[o] = vec
    return pd.DataFrame.from_dict(rows, orient="index", columns=list(layers))


def assay_vaf_matrix(carriers: pd.DataFrame, contributions: pd.DataFrame, phase: str = "unphased"):
    """Site x organ assay-VAF matrix. Returns (sites DataFrame [mutation_id, layer], organs, V, F)."""
    if phase not in PHASE_FACTORS:
        raise ValueError(f"phase must be one of {list(PHASE_FACTORS)}")
    organs = list(contributions.index)
    piv = carriers.pivot_table(index=["mutation_id", "layer"], columns="organ_id", values="carrier_fraction",
                               fill_value=0.0).reindex(columns=organs, fill_value=0.0)
    lay = piv.index.get_level_values("layer")
    c = contributions.to_numpy()[:, [list(contributions.columns).index(k) for k in lay]].T   # sites x organs
    F = piv.to_numpy()
    V = PHASE_FACTORS[phase] * c * F
    sites = piv.index.to_frame(index=False)
    return sites, organs, V, F


def observe_layered(carriers: pd.DataFrame, contributions: Any, observation_config: Mapping[str, Any], *,
                    phase: str = "unphased", seed: int = 0, layers: Optional[Sequence[str]] = None):
    """Read-level observation of layered organs. Returns (evidence long table, observe() result,
    sites, organs)."""
    layers = list(layers) if layers is not None else sorted(carriers.layer.unique())
    organs = sorted(carriers.organ_id.unique())
    cmat = resolve_contributions(contributions, organs, layers)
    sites, organs, V, F = assay_vaf_matrix(carriers, cmat, phase)
    rng = np.random.default_rng(int(seed))
    obs = plantsoma_obs.observe(V, observation_config, rng)
    n_som, n_bg = len(sites), int(obs["is_background"].sum())
    site_id = list(sites.mutation_id) + [f"bg:{i}" for i in range(n_bg)]
    layer = list(sites.layer) + ["background"] * n_bg
    F_all = np.vstack([F, np.zeros((n_bg, len(organs)))]) if n_bg else F
    k = len(organs)
    ev = pd.DataFrame({
        "site_id": np.repeat(site_id, k), "site_layer": np.repeat(layer, k), "organ_id": np.tile(organs, len(site_id)),
        "carrier_fraction": F_all.ravel(), "assay_vaf": obs["assay_vaf"].ravel(), "depth": obs["depth"].ravel(),
        "alt_count": obs["alt"].ravel(), "observed_vaf": obs["observed_vaf"].ravel(),
        "called": obs["called"].ravel(), "ascertained": np.repeat(obs["keep"], k),
    })
    return ev, obs, site_id, organs


# --------------------------------------------------------------------------- config-driven run
LAYERED_CONFIG_TEMPLATE = {
    "_about": ("Template for `simsoma layers`. Layer contributions are the mean of published leaf "
               "compositions; per-layer mutation rates are illustrative relative values. Adapt both to "
               "your tissue. "
               "layer_contributions can also be given per organ: {organ: {layer: c}, 'default': {...}}. "
               "Sources and details: simSOMA_docs/layered_simulation.md"),
    "topology": {"topology_json": "path/to/topology.json", "mapping_unit": "years", "mapping_rate": 5.0,
                 "phyllotaxy": {"mode": "off"}},
    "simulation": {"seed": 0, "n_replicates": 1,
                   "parameters": {"m": 3, "rho": 0.0, "sam_boundary_cells": 40, "branch_precursor_number": 3,
                                  "organ_precursor_number": 10, "organ_total_cells": 1000, "sequenced_cells": 1000,
                                  "organ_mu_multiplier": 1.0},
                   "layers": {"L1": {"mu_unit": 0.95}, "L2": {"mu_unit": 0.55}, "L3": {"mu_unit": 0.55}}},
    "observation": {"phase": "unphased", "layer_contributions": {"L1": 0.13, "L2": 0.84, "L3": 0.03},
                    "model": {"depth": {"mode": "lognormal_site_sample", "mean": 60},
                              "reads": {"type": "binomial", "sequencing_error": 0.0},
                              "background": {"n_sites": 2000, "distribution": "gamma:2", "mean": 0.023},
                              "caller": {"min_depth": 0, "min_alt_reads": 2, "retain_called_any": True}},
                    "vafsoma_tables": {"first_cut_fraction": 0.5, "n_tiers": 4}},
    "output": {"dir": "simSOMA_output/layered_example"},
}

_REQUIRED_PARAMS = ("m", "rho", "sam_boundary_cells", "branch_precursor_number", "organ_precursor_number",
                    "organ_total_cells")


def validate_layered_config(cfg: Mapping[str, Any]) -> None:
    for sec in ("topology", "simulation", "output"):
        if sec not in cfg:
            raise ValueError(f"layered config: missing section '{sec}'")
    if "topology_json" not in cfg["topology"]:
        raise ValueError("layered config: topology.topology_json is required")
    p = cfg["simulation"].get("parameters", {})
    miss = [k for k in _REQUIRED_PARAMS if k not in p]
    if miss:
        raise ValueError(f"layered config: simulation.parameters missing {miss}")
    if not cfg["simulation"].get("layers"):
        raise ValueError("layered config: simulation.layers must name at least one layer with mu_unit")
    if "observation" in cfg:
        ob = cfg["observation"]
        if "layer_contributions" not in ob:
            raise ValueError("layered config: observation.layer_contributions is required")
        plantsoma_obs.normalize_config(ob.get("model"))


def run_layered_config(config_path: str | Path) -> Dict[str, Any]:
    config_path = Path(config_path)
    cfg = json.loads(config_path.read_text())
    validate_layered_config(cfg)
    base = config_path.parent
    tcfg = cfg["topology"]
    import run_from_config as _rfc          # same path rules as `simsoma run`
    tpath = _rfc._resolve_topology_json_path(str(tcfg["topology_json"]), base_dir=base)
    if not tpath.is_file():
        raise FileNotFoundError(f"topology.topology_json not found: {tpath}\n"
                                "Relative paths are read relative to the folder of the config file.")
    unit = str(tcfg.get("mapping_unit", "steps"))
    rate = float(tcfg.get("mapping_rate", 1.0))
    mapping = None if unit == "steps" else {"unit": unit, "rate": rate, "mode": "deterministic"}
    phyl = tcfg.get("phyllotaxy") or {"mode": "off"}
    sim = cfg["simulation"]; P = dict(sim["parameters"])
    seed, n_rep = int(sim.get("seed", 0)), int(sim.get("n_replicates", 1))
    topo = topology_io.load_topology_auto(tpath, mapping=mapping, phyllotaxy_config=phyl,
                                          phyllotaxy_seed=layer_seed(seed, 0, "phyllotaxy"))
    out = _rfc._resolve_outdir_root(str(cfg["output"]["dir"]), base_dir=base)
    out.mkdir(parents=True, exist_ok=True)
    kappa = rate if unit != "steps" else float(P.get("kappa_sr", 1.0))
    written = []
    for rep in range(n_rep):
        rdir = out / f"replicate_{rep:04d}"; rdir.mkdir(exist_ok=True)
        car = simulate_layers(topo, sim["layers"], m=P["m"], rho=P["rho"], kappa_sr=kappa,
                              sam_boundary_cells=P["sam_boundary_cells"],
                              branch_precursor_number=P["branch_precursor_number"],
                              organ_precursor_number=P["organ_precursor_number"],
                              organ_total_cells=P["organ_total_cells"], sequenced_cells=P.get("sequenced_cells"),
                              organ_mu_multiplier=P.get("organ_mu_multiplier", 1.0),
                              victim_locality=P.get("victim_locality", 0.0), seed=seed, replicate=rep)
        car.to_csv(rdir / "layer_carriers.csv.gz", index=False); written.append(str(rdir / "layer_carriers.csv.gz"))
        if "observation" in cfg:
            ob = cfg["observation"]
            ev, obs, site_id, organs = observe_layered(car, ob["layer_contributions"], ob.get("model", {}),
                                                       phase=ob.get("phase", "unphased"),
                                                       seed=layer_seed(seed, rep, "observation"),
                                                       layers=list(sim["layers"]))
            ev.to_csv(rdir / "read_evidence.csv.gz", index=False); written.append(str(rdir / "read_evidence.csv.gz"))
            vt = ob.get("vafsoma_tables")
            if vt:
                mean_depth = float(obs["config"]["depth"]["mean"])
                cut0 = int(vt.get("first_cut", round(float(vt.get("first_cut_fraction", 0.5)) * mean_depth)))
                tabs = plantsoma_obs.depth_tier_tables(obs, organs, site_id, first_cut=cut0,
                                                       n_tiers=int(vt.get("n_tiers", 4)),
                                                       site_filter_all_samples=bool(vt.get("site_filter_all_samples", True)))
                for cut, df in tabs.items():
                    f = rdir / f"vafsoma_dp_{cut}_vaf.csv"; df.to_csv(f, index=False); written.append(str(f))
    settings = {"layered_api_version": LAYERED_API_VERSION, "config": cfg,
                "plantsoma_obs_version": plantsoma_obs.OBSERVATION_MODEL_VERSION,
                "observation_config_sha256": (plantsoma_obs.config_sha256(cfg["observation"].get("model", {}))
                                              if "observation" in cfg else None)}
    try:
        import simsoma
        settings["simsoma_version"] = simsoma.__version__
    except Exception:
        pass
    (out / "layered_settings.json").write_text(json.dumps(settings, indent=2))
    return {"output_dir": str(out), "files": written}
