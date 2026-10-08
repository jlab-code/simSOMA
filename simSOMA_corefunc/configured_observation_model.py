#!/usr/bin/env python3
"""Config-driven observation-model and fitSOMA handoff for simSOMA.

This module intentionally leaves the developmental simulator unchanged.  It
reads the ordinary ``grid_parameter`` outputs, applies the existing deterministic
layer/phasing transform, and optionally exports complete mutation-by-sample
files for fitSOMA.  Two public modes are supported:

``deterministic``
    Exact layer/phasing transformation without stochastic read noise.
``read_counts``
    The same deterministic assay transformation followed by depth/read/error
    sampling and caller emulation through ``transform_vaf_observations``.

If the top-level ``observation_model`` config block is absent, this module is not
called and the legacy simSOMA run is byte-for-byte unaffected apart from normal
runtime metadata such as timestamps.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd

from transform_observation_model import (
    run_observation_model_transform,
    transform_vaf_observations,
)

VALID_MODES = {"deterministic", "read_counts"}
VALID_SAMPLING = {"bulk", "layer_specific"}
VALID_PHASE = {"phased", "unphased"}


def _safe_name(value: Any) -> str:
    import re
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())
    text = re.sub(r"_+", "_", text).strip("_")
    if text in {"", ".", ".."}:
        raise ValueError("observation_model.scenario_name is invalid")
    return text


def _as_list(value: Any, *, name: str) -> list[Any]:
    if isinstance(value, str):
        out = [x.strip() for x in value.split(",") if x.strip()]
    elif isinstance(value, (list, tuple)):
        out = list(value)
    else:
        raise ValueError(f"{name} must be a list or comma-separated string")
    if not out:
        raise ValueError(f"{name} must not be empty")
    return out


def normalize_observation_model_config(block: Mapping[str, Any]) -> dict[str, Any]:
    if not isinstance(block, Mapping):
        raise ValueError("observation_model must be a JSON object")
    mode = str(block.get("mode", "deterministic")).strip().lower()
    if mode not in VALID_MODES:
        raise ValueError(f"observation_model.mode must be one of {sorted(VALID_MODES)}")

    layers = [str(x) for x in _as_list(block.get("layers", ["layer_equivalent"]), name="observation_model.layers")]
    weights = [float(x) for x in _as_list(block.get("layer_weights", [1.0]), name="observation_model.layer_weights")]
    if len(layers) != len(weights):
        raise ValueError("observation_model.layers and layer_weights must have the same length")
    if len(set(layers)) != len(layers):
        raise ValueError("observation_model.layers must be unique")
    if any((not math.isfinite(x) or x < 0.0) for x in weights):
        raise ValueError("observation_model.layer_weights must be finite and non-negative")
    if not math.isclose(sum(weights), 1.0, rel_tol=1e-9, abs_tol=1e-9):
        raise ValueError("observation_model.layer_weights must sum to 1.0")

    sampling = str(block.get("sampling", block.get("sampling_mode", "layer_specific"))).strip().lower()
    phase = str(block.get("phase", block.get("phase_mode", "phased"))).strip().lower()
    if sampling not in VALID_SAMPLING:
        raise ValueError(f"observation_model.sampling must be one of {sorted(VALID_SAMPLING)}")
    if phase not in VALID_PHASE:
        raise ValueError(f"observation_model.phase must be one of {sorted(VALID_PHASE)}")
    target = block.get("target_layer")
    if sampling == "layer_specific":
        target = str(target if target is not None else layers[0])
        if target not in layers:
            raise ValueError("observation_model.target_layer must be present in layers")
    else:
        target = None

    scenario = _safe_name(block.get("scenario_name", f"{mode}_{sampling}_{phase}"))
    deterministic_depth = int(block.get("deterministic_depth", 1_000_000))
    if deterministic_depth < 1:
        raise ValueError("observation_model.deterministic_depth must be >=1")

    read_cfg = dict(block.get("read_counts", {}))
    # Also accept the more explicit top-level subblocks proposed for the public config.
    if "depth" in block:
        read_cfg.setdefault("depth", block["depth"])
    if "read_sampling" in block:
        read_cfg.setdefault("read_sampling", block["read_sampling"])
    if "sequencing_error" in block:
        read_cfg.setdefault("sequencing_error", block["sequencing_error"])
    if "caller" in block:
        read_cfg.setdefault("caller", block["caller"])
    if "background" in block:
        read_cfg.setdefault("background", block["background"])
    # Ascertainment may be written in any of the public forms below.  The
    # nested read_counts value remains authoritative when multiple aliases are
    # present, preserving backward compatibility.
    if "ascertainment" in block and isinstance(block["ascertainment"], Mapping):
        ascertainment = dict(block["ascertainment"])
        if "retain_called_any" in ascertainment:
            read_cfg.setdefault("retain_called_any", ascertainment["retain_called_any"])
    if "retain_called_any" in block:
        read_cfg.setdefault("retain_called_any", block["retain_called_any"])

    return {
        "mode": mode,
        "scenario_name": scenario,
        "layers": layers,
        "layer_weights": weights,
        "sampling": sampling,
        "target_layer": target,
        "phase": phase,
        "export_fitsoma": bool(block.get("export_fitsoma", False)),
        "overwrite": bool(block.get("overwrite", False)),
        "deterministic_depth": deterministic_depth,
        "seed": int(block.get("seed", 910_000_000)),
        "read_counts": read_cfg,
    }


def observation_requires_raw_vafs(raw_config: Mapping[str, Any]) -> bool:
    block = raw_config.get("observation_model") if isinstance(raw_config, Mapping) else None
    if block is None:
        return False
    return bool(normalize_observation_model_config(block)["export_fitsoma"])


def _sample_info(organs: list[str], cfg: Mapping[str, Any]) -> pd.DataFrame:
    rows = []
    for organ in organs:
        row: dict[str, Any] = {
            "sample_id": str(organ),
            "topology_sample_id": str(organ),
            "sampling_mode": cfg["sampling"],
            "phase_mode": cfg["phase"],
            "target_layer": cfg["target_layer"] or "",
        }
        for layer, weight in zip(cfg["layers"], cfg["layer_weights"]):
            row[f"{layer}_weight"] = float(weight)
        rows.append(row)
    return pd.DataFrame(rows)


def _contribution(cfg: Mapping[str, Any], layer: str) -> float:
    if cfg["sampling"] == "layer_specific":
        return 1.0 if layer == cfg["target_layer"] else 0.0
    return float(dict(zip(cfg["layers"], cfg["layer_weights"]))[layer])


def _phase_factor(phase: str) -> float:
    return 1.0 if phase == "phased" else 0.5


def _latent_for_layer(raw_group: pd.DataFrame, organs: list[str], layer: str) -> pd.DataFrame:
    if raw_group.empty:
        return pd.DataFrame(columns=["mutation_id", "topology_sample_id", "latent_vaf"])
    x = raw_group[["mutation_id", "organ_id", "allele_count", "sequenced_cells"]].copy()
    x["mutation_id"] = layer + "::" + x["mutation_id"].astype(str)
    x["topology_sample_id"] = x["organ_id"].astype(str)
    x["latent_vaf"] = x["allele_count"].astype(float) / x["sequenced_cells"].astype(float)
    return x[["mutation_id", "topology_sample_id", "latent_vaf"]]


def _deterministic_evidence(
    raw_group: pd.DataFrame,
    sample_info: pd.DataFrame,
    cfg: Mapping[str, Any],
) -> pd.DataFrame:
    depth = int(cfg["deterministic_depth"])
    phase = _phase_factor(str(cfg["phase"]))
    parts: list[pd.DataFrame] = []
    for layer in cfg["layers"]:
        contribution = _contribution(cfg, layer)
        if contribution <= 0.0:
            continue
        latent = _latent_for_layer(raw_group, sample_info["topology_sample_id"].tolist(), layer)
        if latent.empty:
            continue
        mids = sorted(latent["mutation_id"].unique())
        pivot = latent.pivot_table(index="mutation_id", columns="topology_sample_id", values="latent_vaf", aggfunc="max", fill_value=0.0)
        pivot = pivot.reindex(index=mids, columns=sample_info["topology_sample_id"].tolist(), fill_value=0.0)
        assay = np.clip(pivot.to_numpy(float) * contribution * phase, 0.0, 1.0)
        alt = np.rint(assay * depth).astype(int)
        rows: list[dict[str, Any]] = []
        for i, mid in enumerate(mids):
            for j, sample in enumerate(sample_info.itertuples(index=False)):
                v = float(assay[i, j])
                rows.append({
                    "mutation_id": str(mid),
                    "sample_id": str(sample.sample_id),
                    "topology_sample_id": str(sample.topology_sample_id),
                    "alt_count": int(alt[i, j]),
                    "ref_count": int(depth - alt[i, j]),
                    "depth": int(depth),
                    "observed_vaf": float(alt[i, j]) / float(depth),
                    "exact_vaf": v,
                    "callable": True,
                    "caller_call": bool(v > 0.0),
                    "caller_status": "CALLED" if v > 0.0 else "ABSENT",
                    "source_layer": str(layer),
                    "observation_mode": "deterministic",
                })
        parts.append(pd.DataFrame(rows))
    if not parts:
        return pd.DataFrame(columns=["mutation_id", "sample_id", "alt_count", "depth", "callable", "caller_call", "caller_status"])
    return pd.concat(parts, ignore_index=True)


def _read_count_config(cfg: Mapping[str, Any], layer: str) -> dict[str, Any]:
    rc = dict(cfg.get("read_counts", {}))
    depth = dict(rc.get("depth", {}))
    read = dict(rc.get("read_sampling", {}))
    error = rc.get("sequencing_error", 0.001)
    caller = dict(rc.get("caller", {}))
    if isinstance(error, Mapping):
        e1 = float(error.get("reference_to_alternate", error.get("rate", 0.001)))
        e2 = float(error.get("alternate_to_reference", e1))
        if not math.isclose(e1, e2, rel_tol=0.0, abs_tol=0.0):
            raise ValueError("Current simSOMA read_counts mode requires symmetric sequencing-error rates")
        error = e1
    mode = str(depth.get("mode", "fixed"))
    value = float(depth.get("value", depth.get("default_mean", 100.0)))
    sd = float(depth.get("sd", depth.get("default_sd", 0.0)))
    depth_model = {
        # historical mapping: anything other than "fixed" meant per-sample normal depth
        "mode": mode if mode in ("fixed", "poisson", "lognormal_site_sample") else "sample_design",
        "default_mean": value,
        "default_sd": sd,
        "minimum": int(depth.get("minimum", 1)),
    }
    for k in ("site_sdlog", "sample_sdlog", "sample_factor_sdlog", "max_site_factor"):
        if k in depth:
            depth_model[k] = float(depth[k])
    out = {
        "assay": {"source_layer": layer},
        "read_model": {
            "type": str(read.get("distribution", read.get("type", "beta_binomial"))),
            "sequencing_error": float(error),
            "concentration": float(read.get("concentration", 200.0)),
        },
        "depth_model": depth_model,
        "caller": {
            "min_depth": int(caller.get("minimum_depth", caller.get("min_depth", 20))),
            "min_alt_reads": int(caller.get("minimum_alt_reads", caller.get("min_alt_reads", 3))),
            "min_observed_vaf": float(caller.get("minimum_vaf", caller.get("min_observed_vaf", 0.0))),
        },
        "ascertainment": {"retain_called_any": bool(rc.get("retain_called_any", True))},
    }
    # Background artefact sites belong to the sample, not to a layer: generated once, with the
    # first contributing layer (see _read_count_evidence).
    bg = rc.get("background")
    if bg and int(dict(bg).get("n_sites", 0)) > 0:
        first = next((l for l in cfg["layers"] if _contribution(cfg, l) > 0.0), None)
        if layer == first:
            out["background"] = dict(bg)
    return out


def _read_count_evidence(raw_group: pd.DataFrame, sample_info: pd.DataFrame, cfg: Mapping[str, Any], seed: int) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    parts: list[pd.DataFrame] = []
    metadata: list[dict[str, Any]] = []
    for li, layer in enumerate(cfg["layers"]):
        if _contribution(cfg, layer) <= 0.0:
            continue
        latent = _latent_for_layer(raw_group, sample_info["topology_sample_id"].tolist(), layer)
        if latent.empty:
            continue
        ev, meta = transform_vaf_observations(latent, sample_info, _read_count_config(cfg, layer), seed=int(seed) + li * 100003)
        if not ev.empty:
            ev["source_layer"] = str(layer)
            ev["observation_mode"] = "read_counts"
            parts.append(ev)
        metadata.append(meta)
    if not parts:
        return pd.DataFrame(columns=["mutation_id", "sample_id", "alt_count", "depth", "callable", "caller_call", "caller_status"]), metadata
    return pd.concat(parts, ignore_index=True), metadata


def _jsonable(value: Any) -> Any:
    if pd.isna(value):
        return None
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    return value


def _canonical_truth_aliases(configured: Mapping[str, Any], run_config: Mapping[str, Any], obs_cfg: Mapping[str, Any]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    mu = configured.get("mu_unit", configured.get("mu_year"))
    if mu is not None:
        out["mu_unit"] = float(mu)
    k = float(configured.get("effective_kappa_sr", configured.get("mapping_rate", 1.0)))
    out["kappa_sr"] = k
    out["K_sr"] = k
    if configured.get("sam_boundary_cells") is not None:
        out["C"] = int(float(configured["sam_boundary_cells"]))
    if configured.get("branch_precursor_number_realized_cells") is not None:
        out["P_b_eff"] = int(float(configured["branch_precursor_number_realized_cells"]))
    elif configured.get("branch_precursor_number") is not None:
        out["P_b_eff"] = min(int(float(configured["branch_precursor_number"])), int(float(configured.get("m", 1))))
    if "P_b_eff" in out:
        out["P_a_eff"] = out["P_b_eff"]  # deprecated alias (fitSOMA <= 0.3.19); remove later
    if configured.get("organ_precursor_number_realized_cells") is not None:
        out["P_o"] = int(float(configured["organ_precursor_number_realized_cells"]))
    elif configured.get("organ_precursor_number") is not None:
        out["P_o"] = int(float(configured["organ_precursor_number"]))
    # Expected clonal composition of founders (local_sector_v1; exact for uniform focal
    # placement, i.e. no phyllotaxy): phi = (P_eff - 1) m / C, pi = Pr(polyclonal founding).
    try:
        from founder_diversity import FOUNDER_DEFINITION_VERSION, polyclonal_founding_expectation
        out["founder_definition"] = FOUNDER_DEFINITION_VERSION
        m_i = int(float(configured.get("m", 0)))
        C_i = out.get("C")
        for key, P in (("B", out.get("P_b_eff")), ("O", out.get("P_o"))):
            if m_i >= 1 and C_i and P and int(P) <= int(C_i):
                e = polyclonal_founding_expectation(int(P), int(C_i), m_i)
                out[f"phi_{key}"] = e["phi"]
                out[f"pi_{key}_expected"] = e["probability"]
                out[f"sector_count_{key}_expected"] = e["expected_sector_count"]
    except ImportError:
        pass
    rho = float(configured.get("rho", 0.0))
    rho = min(max(rho, 0.0), 1.0 - 1e-15)
    lam = 0.0 if rho <= 0.0 else -k * math.log1p(-rho)
    out["rho_turn"] = rho
    out["lambda_turn"] = lam
    out["p_turn"] = 0.0 if lam <= 0.0 else 1.0 - math.exp(-lam)
    topo = run_config.get("topology", {}) if isinstance(run_config, Mapping) else {}
    phyl = topo.get("phyllotaxy", {}) if isinstance(topo, Mapping) else {}
    out["phyllotaxy_mode"] = str(phyl.get("mode", "off")) if isinstance(phyl, Mapping) else "off"
    if obs_cfg.get("mode") == "read_counts":
        rc = dict(obs_cfg.get("read_counts", {}))
        read = dict(rc.get("read_sampling", {}))
        err = rc.get("sequencing_error", 0.001)
        if isinstance(err, Mapping):
            err = err.get("reference_to_alternate", err.get("rate", 0.001))
        out["sequencing_error"] = float(err)
        out["read_concentration"] = float(read.get("concentration", 200.0))
    else:
        out["sequencing_error"] = 0.0
        out["read_concentration"] = None
    return out


def _realized_summary(events: pd.DataFrame) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for event_type, prefix in [("branch", "B"), ("organ", "O")]:
        d = events.loc[events["event_type"].eq(event_type)] if not events.empty else events
        if d.empty:
            out[f"pi_{prefix}"] = None
            out[f"d_{prefix}"] = None
            out[f"{event_type}_founder_sector_mean"] = None
            continue
        sectors = pd.to_numeric(d["founder_sector_count"], errors="coerce")
        diversity = pd.to_numeric(d["founder_diversity"], errors="coerce")
        # local_sector_v1: realized fraction of polyclonal foundings (see founder_diversity.py)
        out[f"pi_{prefix}"] = float((sectors > 1).mean())
        out[f"d_{prefix}"] = float(diversity.mean())
        out[f"{event_type}_founder_sector_mean"] = float(sectors.mean())
    return out


def _observation_contract(cfg: Mapping[str, Any]) -> dict[str, Any]:
    """Return the canonical generative observation contract used by fitSOMA."""
    active_layers = [layer for layer in cfg["layers"] if _contribution(cfg, layer) > 0.0]
    contract: dict[str, Any] = {
        "schema_version": "1.0",
        "observation_mode": str(cfg["mode"]),
        "sampling_modes": [str(cfg["sampling"])],
        "phase_modes": [str(cfg["phase"])],
        "layers": [],
    }
    for layer in active_layers:
        if cfg["mode"] == "read_counts":
            layer_cfg = _read_count_config(cfg, layer)
            contract["layers"].append({
                "source_layer": str(layer),
                "read_model": dict(layer_cfg["read_model"]),
                "depth_model": dict(layer_cfg["depth_model"]),
                "caller": dict(layer_cfg["caller"]),
                "ascertainment": dict(layer_cfg["ascertainment"]),
                **({"background": dict(layer_cfg["background"])} if layer_cfg.get("background") else {}),
            })
        else:
            contract["layers"].append({
                "source_layer": str(layer),
                "deterministic_depth": int(cfg["deterministic_depth"]),
            })
    return contract


def _contract_sha256(contract: Mapping[str, Any]) -> str:
    payload = json.dumps(dict(contract), sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def export_fitsoma_handoff(run_dir: Path, scenario_dir: Path, cfg: Mapping[str, Any]) -> dict[str, Any]:
    grid = run_dir / "grid_parameter"
    raw_path = grid / "raw_vafs.csv.gz"
    params_path = grid / "parameter_sets.csv"
    organs_path = grid / "organ_replicate_summaries.csv"
    events_path = grid / "realized_event_truth.csv.gz"
    if not raw_path.exists():
        raise FileNotFoundError("fitSOMA export requires grid_parameter/raw_vafs.csv.gz")
    raw = pd.read_csv(raw_path)
    params = pd.read_csv(params_path)
    organs_table = pd.read_csv(organs_path)
    events = pd.read_csv(events_path) if events_path.exists() else pd.DataFrame()
    config_used_path = run_dir / "config_used.json"
    run_config = json.loads(config_used_path.read_text(encoding="utf-8")) if config_used_path.exists() else {}
    out_root = scenario_dir / "fitSOMA"
    out_root.mkdir(parents=True, exist_ok=True)
    manifest_rows: list[dict[str, Any]] = []

    pairs = organs_table[["set_id", "rep"]].drop_duplicates().sort_values(["set_id", "rep"])
    for pair in pairs.itertuples(index=False):
        set_id, rep = str(pair.set_id), int(pair.rep)
        organs = sorted(organs_table.loc[(organs_table["set_id"].astype(str) == set_id) & (organs_table["rep"].astype(int) == rep), "organ_id"].astype(str).unique())
        sample_info = _sample_info(organs, cfg)
        rg = raw.loc[(raw["set_id"].astype(str) == set_id) & (raw["rep"].astype(int) == rep)].copy()
        obs_seed = int(cfg["seed"]) + int(set_id.split("_")[-1]) * 1_000_000 + rep
        if cfg["mode"] == "deterministic":
            evidence = _deterministic_evidence(rg, sample_info, cfg)
            read_meta: Any = None
        else:
            evidence, read_meta = _read_count_evidence(rg, sample_info, cfg, obs_seed)

        rep_dir = out_root / set_id / f"replicate_{rep:05d}"
        rep_dir.mkdir(parents=True, exist_ok=True)
        evidence.to_csv(rep_dir / "variants.tsv.gz", sep="\t", index=False, compression="gzip")
        sample_info.to_csv(rep_dir / "sample_info.tsv", sep="\t", index=False)
        prow = params.loc[params["set_id"].astype(str).eq(set_id)]
        configured = {k: _jsonable(v) for k, v in prow.iloc[0].to_dict().items()} if len(prow) else {"set_id": set_id}
        configured.update(_canonical_truth_aliases(configured, run_config, cfg))
        (rep_dir / "configured_truth.json").write_text(json.dumps(configured, indent=2, sort_keys=True), encoding="utf-8")
        er = events.loc[(events["set_id"].astype(str) == set_id) & (events["rep"].astype(int) == rep)].copy() if not events.empty else events
        realized = {**configured, **_realized_summary(er), "set_id": set_id, "rep": rep, "event_truth_file": "realized_event_truth.tsv" if not er.empty else None}
        (rep_dir / "realized_truth.json").write_text(json.dumps(realized, indent=2, sort_keys=True), encoding="utf-8")
        if not er.empty:
            er.to_csv(rep_dir / "realized_event_truth.tsv", sep="\t", index=False)
        observation_contract = _observation_contract(cfg)
        provenance = {
            "provenance_schema_version": "1.1",
            "generator": "simSOMA",
            "observation_mode": cfg["mode"],
            "scenario_name": cfg["scenario_name"],
            "observation_seed": obs_seed,
            "complete_mutation_by_sample_matrix": True,
            "ubiquitous_filter_applied": False,
            "observation_contract": observation_contract,
            "observation_contract_sha256": _contract_sha256(observation_contract),
            "read_metadata": read_meta,
        }
        (rep_dir / "provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True), encoding="utf-8")
        manifest_rows.append({"set_id": set_id, "rep": rep, "n_samples": len(sample_info), "n_candidates": int(evidence["mutation_id"].nunique()) if not evidence.empty else 0, "directory": str(rep_dir.relative_to(scenario_dir))})

    pd.DataFrame(manifest_rows).to_csv(out_root / "handoff_manifest.tsv", sep="\t", index=False)
    return {"n_handoffs": len(manifest_rows), "handoff_root": str(out_root)}


def run_configured_observation_model(run_dir: str | Path, raw_config: Mapping[str, Any], *, overwrite: bool | None = None) -> dict[str, Any] | None:
    block = raw_config.get("observation_model") if isinstance(raw_config, Mapping) else None
    if block is None:
        return None
    cfg = normalize_observation_model_config(block)
    if overwrite is not None:
        cfg["overwrite"] = bool(overwrite)
    run_dir = Path(run_dir).expanduser().resolve()
    result = run_observation_model_transform(
        run_or_grid_dir=run_dir,
        scenario_name=cfg["scenario_name"],
        layer_names=cfg["layers"],
        layer_weights=cfg["layer_weights"],
        sampling_mode=cfg["sampling"],
        phase_mode=cfg["phase"],
        target_layers=cfg["target_layer"] if cfg["sampling"] == "layer_specific" else "all",
        overwrite=cfg["overwrite"],
    )
    scenario_dir = run_dir / "observation_model_transforms" / cfg["scenario_name"]
    resolved = dict(cfg)
    (scenario_dir / "observation_model_config_resolved.json").write_text(json.dumps(resolved, indent=2, sort_keys=True), encoding="utf-8")
    handoff = export_fitsoma_handoff(run_dir, scenario_dir, cfg) if cfg["export_fitsoma"] else None
    report = {"scenario_dir": str(scenario_dir), "mode": cfg["mode"], "summary_transform": result, "fitSOMA_handoff": handoff}
    (scenario_dir / "configured_observation_report.json").write_text(json.dumps(report, indent=2, sort_keys=True, default=str), encoding="utf-8")
    return report
