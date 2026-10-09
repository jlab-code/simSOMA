#!/usr/bin/env python3
"""
transform_observation_model.py

Standalone post-processing tool for simSOMA outputs.

This tool transforms an existing layer-equivalent simSOMA VAF/count spectrum
into an observation-space spectrum under one user-specified observation scenario.
All observation-model hyperparameters are passed directly to the Python function
or through command-line arguments. No observation-model JSON file is used.

Default deterministic transform:

    observed_vaf = source_vaf * effective_layer_contribution * phase_factor

where:
  - effective_layer_contribution is the layer weight for bulk samples.
  - effective_layer_contribution is 1.0 for layer-specific samples.
  - phase_factor is 1.0 for phased observations and 0.5 for unphased observations.

Layer interpretation:
  The source simSOMA output is treated as one layer-equivalent developmental
  spectrum. For a multi-layer transform, the source spectrum is copied to each
  modeled layer, but copied rows represent distinct layer-specific mutation
  identities under an infinite-sites interpretation. A copied row is not the
  same mutation present in multiple layers.

Expected standard input files under <run-dir>/grid_parameter/:
  - vaf_count_spectrum_aggregated_summaries.csv
  - vaf_count_spectrum_replicate_summaries.csv.gz

Standard output folder:
  <run-dir>/observation_model_transforms/<scenario_name>/

No JSON input or output is written. Metadata and manifest files are TSV/TXT.
"""

from __future__ import annotations

import argparse
import gzip
import math
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Any, Mapping

import numpy as np
import pandas as pd

STD_AGG_FILE = "vaf_count_spectrum_aggregated_summaries.csv"
STD_REP_FILE = "vaf_count_spectrum_replicate_summaries.csv.gz"

VALID_SAMPLING_MODES = {"bulk", "layer_specific"}
VALID_PHASE_MODES = {"phased", "unphased"}
VALID_LAYER_WEIGHT_MODES = {"deterministic"}

COUNT_COLS_REPLICATE = ["n_variants"]
COUNT_COLS_AGGREGATED = ["n_variants_mean", "n_variants_sd"]


def _now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def _safe_name(s: str) -> str:
    s = str(s).strip()
    if not s:
        raise ValueError("scenario_name must not be empty")
    s = re.sub(r"[^A-Za-z0-9_.-]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    if s in {"", ".", ".."}:
        raise ValueError("scenario_name resolves to an invalid output folder name")
    return s


def _infer_sep(path: Path) -> str:
    name = path.name.lower()
    if name.endswith(".tsv") or name.endswith(".tsv.gz"):
        return "\t"
    return ","


def _read_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing input table: {path}")
    return pd.read_csv(path, sep=_infer_sep(path))


def _write_table(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    sep = "\t" if (path.name.endswith(".tsv") or path.name.endswith(".tsv.gz")) else ","
    compression = "gzip" if path.name.endswith(".gz") else None
    df.to_csv(path, sep=sep, index=False, compression=compression)


def _relative_to(path: Path, base: Path) -> str:
    try:
        return str(path.resolve().relative_to(base.resolve()))
    except Exception:
        return str(path.resolve())


def _resolve_run_and_grid_dir(path: Path) -> tuple[Path, Path]:
    """Return (run_dir, grid_dir). Accepts either experiment dir or grid_parameter dir."""
    p = path.expanduser().resolve()
    if not p.exists():
        raise FileNotFoundError(f"Run/grid directory does not exist: {p}")
    if (p / STD_AGG_FILE).exists() or (p / STD_REP_FILE).exists():
        return p.parent, p
    if (p / "grid_parameter" / STD_AGG_FILE).exists() or (p / "grid_parameter" / STD_REP_FILE).exists():
        return p, p / "grid_parameter"
    raise FileNotFoundError(
        f"Could not find standard simSOMA VAF tables in {p} or {p / 'grid_parameter'}"
    )


def _resolve_optional_table(value: str | Path | None, *, run_dir: Path, grid_dir: Path) -> Path | None:
    if value is None or str(value).strip() == "":
        return None
    raw = Path(str(value)).expanduser()
    if raw.is_absolute():
        return raw.resolve()

    candidates = [
        (run_dir / raw).resolve(),
        (grid_dir / raw).resolve(),
        raw.resolve(),
    ]
    for c in candidates:
        if c.exists():
            return c
    return candidates[0]


def _phase_factor(phase_mode: str) -> float:
    if phase_mode == "phased":
        return 1.0
    if phase_mode == "unphased":
        return 0.5
    raise ValueError(f"Invalid phase_mode: {phase_mode}")


def _parse_csv_str(value: str | Iterable[str], *, name: str) -> list[str]:
    if isinstance(value, str):
        vals = [x.strip() for x in value.split(",") if x.strip()]
    else:
        vals = [str(x).strip() for x in value if str(x).strip()]
    if not vals:
        raise ValueError(f"{name} must not be empty")
    return vals


def _parse_csv_float(value: str | Iterable[float], *, name: str, n: int | None = None) -> list[float]:
    if isinstance(value, str):
        raw_vals = [x.strip() for x in value.split(",") if x.strip()]
    else:
        raw_vals = list(value)
    if not raw_vals:
        raise ValueError(f"{name} must not be empty")
    vals = [float(x) for x in raw_vals]
    if n is not None and len(vals) != n:
        raise ValueError(f"{name} must have length {n}; observed length {len(vals)}")
    if any((not math.isfinite(v)) for v in vals):
        raise ValueError(f"{name} contains non-finite values")
    return vals


def _parse_target_layers(value: str | Iterable[str] | None, layer_names: list[str]) -> list[str]:
    if value is None or value == "" or value == "all":
        return list(layer_names)
    vals = _parse_csv_str(value, name="target_layers")
    unknown = [x for x in vals if x not in layer_names]
    if unknown:
        raise ValueError(f"Unknown target_layers {unknown}; available layers: {layer_names}")
    return vals


def _validate_model_args(
    *,
    scenario_name: str,
    layer_names: str | Iterable[str],
    layer_weights: str | Iterable[float],
    sampling_mode: str,
    phase_mode: str,
    layer_weight_mode: str = "deterministic",
    target_layers: str | Iterable[str] | None = "all",
    layer_mutation_rate_multipliers: str | Iterable[float] | None = None,
    normalize_layer_weights: bool = False,
) -> dict:
    scenario_name = _safe_name(scenario_name)

    layer_names_list = _parse_csv_str(layer_names, name="layer_names")
    if len(set(layer_names_list)) != len(layer_names_list):
        raise ValueError("layer_names must be unique")
    layer_count = len(layer_names_list)

    layer_weights_list = _parse_csv_float(layer_weights, name="layer_weights", n=layer_count)
    if any(v < 0 for v in layer_weights_list):
        raise ValueError("layer_weights must be non-negative")
    weight_sum = float(sum(layer_weights_list))
    if weight_sum <= 0:
        raise ValueError("layer_weights must sum to a positive value")

    if normalize_layer_weights:
        layer_weights_list = [v / weight_sum for v in layer_weights_list]
        weight_sum = 1.0
    elif not math.isclose(weight_sum, 1.0, rel_tol=1e-6, abs_tol=1e-6):
        raise ValueError(
            f"layer_weights must sum to 1.0 for deterministic bulk transforms; observed sum={weight_sum:g}. "
            "Use --normalize-layer-weights to normalize explicitly."
        )

    if sampling_mode not in VALID_SAMPLING_MODES:
        raise ValueError(f"sampling_mode must be one of {sorted(VALID_SAMPLING_MODES)}")
    if phase_mode not in VALID_PHASE_MODES:
        raise ValueError(f"phase_mode must be one of {sorted(VALID_PHASE_MODES)}")
    if layer_weight_mode not in VALID_LAYER_WEIGHT_MODES:
        raise NotImplementedError(
            f"layer_weight_mode='{layer_weight_mode}' is not implemented. "
            f"Implemented modes: {sorted(VALID_LAYER_WEIGHT_MODES)}"
        )

    if layer_mutation_rate_multipliers is None or layer_mutation_rate_multipliers == "":
        multipliers = [1.0] * layer_count
    else:
        multipliers = _parse_csv_float(
            layer_mutation_rate_multipliers,
            name="layer_mutation_rate_multipliers",
            n=layer_count,
        )
    if any(v < 0 for v in multipliers):
        raise ValueError("layer_mutation_rate_multipliers must be non-negative")

    selected_target_layers = _parse_target_layers(target_layers, layer_names_list)

    return {
        "scenario_name": scenario_name,
        "layer_count": layer_count,
        "layer_names": layer_names_list,
        "layer_weights": layer_weights_list,
        "layer_weight_sum": weight_sum,
        "sampling_mode": sampling_mode,
        "phase_mode": phase_mode,
        "phase_factor": _phase_factor(phase_mode),
        "layer_weight_mode": layer_weight_mode,
        "target_layers": selected_target_layers,
        "layer_mutation_rate_multipliers": multipliers,
        "normalize_layer_weights": bool(normalize_layer_weights),
    }


def _numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def _check_required_columns(df: pd.DataFrame, path: Path) -> None:
    required = {"n_sampled_cells", "allele_count", "sampled_vaf"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"Input table {path} is missing required columns: {missing}")


def _observed_vaf_bin(v: pd.Series, nbins: int) -> tuple[pd.Series, pd.Series, pd.Series, pd.Series]:
    vv = pd.to_numeric(v, errors="coerce")
    idx = np.floor(vv * nbins).astype("float")
    idx = idx.mask(vv <= 0, 0)
    idx = idx.mask(vv >= 1, nbins - 1)
    idx = idx.clip(lower=0, upper=nbins - 1)
    lower = idx / float(nbins)
    upper = (idx + 1) / float(nbins)
    label = lower.map(lambda x: f"{x:.6g}") + "-" + upper.map(lambda x: f"{x:.6g}")
    idx_int = idx.astype("Int64")
    return idx_int, lower, upper, label


def _layer_records(model: dict) -> list[dict]:
    names = model["layer_names"]
    weights = model["layer_weights"]
    multipliers = model["layer_mutation_rate_multipliers"]
    sampling = model["sampling_mode"]

    rows: list[dict] = []
    if sampling == "bulk":
        for i, name in enumerate(names):
            rows.append({
                "source_layer_index": i + 1,
                "source_layer": name,
                "sample_layer": "bulk",
                "source_layer_weight": weights[i],
                "effective_layer_contribution": weights[i],
                "layer_mutation_rate_multiplier": multipliers[i],
            })
    elif sampling == "layer_specific":
        target = set(model["target_layers"])
        for i, name in enumerate(names):
            if name not in target:
                continue
            rows.append({
                "source_layer_index": i + 1,
                "source_layer": name,
                "sample_layer": name,
                "source_layer_weight": weights[i],
                "effective_layer_contribution": 1.0,
                "layer_mutation_rate_multiplier": multipliers[i],
            })
    else:
        raise ValueError(f"Unhandled sampling mode: {sampling}")

    if not rows:
        raise ValueError("No layer records selected for transformation")
    return rows


def _transform_table(
    df: pd.DataFrame,
    *,
    source_path: Path,
    model: dict,
    table_kind: str,
    observed_nbins: int,
) -> pd.DataFrame:
    _check_required_columns(df, source_path)

    base = df.copy()
    base["source_n_sampled_cells"] = _numeric(base["n_sampled_cells"])
    base["source_allele_count"] = _numeric(base["allele_count"])
    base["source_vaf"] = _numeric(base["sampled_vaf"])
    base = base.drop(columns=[c for c in ["n_sampled_cells", "allele_count", "sampled_vaf"] if c in base.columns])

    transformed_pieces: list[pd.DataFrame] = []
    phase_factor = float(model["phase_factor"])
    phase_denominator_factor = 1.0 / phase_factor

    for rec in _layer_records(model):
        d = base.copy()
        for k, v in rec.items():
            d[k] = v

        d["observation_scenario_name"] = model["scenario_name"]
        d["sampling_mode"] = model["sampling_mode"]
        d["phase_mode"] = model["phase_mode"]
        d["phase_factor"] = phase_factor
        d["layer_count"] = int(model["layer_count"])
        d["layer_weight_mode"] = model["layer_weight_mode"]
        d["observation_table_kind"] = table_kind

        d["observed_n_sampled_cells"] = d["source_n_sampled_cells"]
        d["observed_vaf_denominator_expected"] = d["source_n_sampled_cells"] * phase_denominator_factor
        d["observed_alt_count_expected"] = (
            d["source_vaf"] * d["effective_layer_contribution"] * d["source_n_sampled_cells"]
        )
        d["observed_vaf"] = d["source_vaf"] * d["effective_layer_contribution"] * phase_factor

        bin_idx, bin_low, bin_high, bin_label = _observed_vaf_bin(d["observed_vaf"], observed_nbins)
        d["observed_vaf_bin_index"] = bin_idx
        d["observed_vaf_bin_lower"] = bin_low
        d["observed_vaf_bin_upper"] = bin_high
        d["observed_vaf_bin"] = bin_label

        multiplier = float(rec["layer_mutation_rate_multiplier"])
        for c in COUNT_COLS_REPLICATE + COUNT_COLS_AGGREGATED:
            if c in d.columns:
                d[f"source_{c}"] = _numeric(d[c])
                d[c] = _numeric(d[c]) * multiplier

        transformed_pieces.append(d)

    out = pd.concat(transformed_pieces, ignore_index=True)

    first_cols = [
        "observation_scenario_name",
        "sampling_mode",
        "phase_mode",
        "phase_factor",
        "layer_weight_mode",
        "layer_count",
        "source_layer_index",
        "source_layer",
        "sample_layer",
        "source_layer_weight",
        "effective_layer_contribution",
        "layer_mutation_rate_multiplier",
        "source_n_sampled_cells",
        "source_allele_count",
        "source_vaf",
        "observed_n_sampled_cells",
        "observed_vaf_denominator_expected",
        "observed_alt_count_expected",
        "observed_vaf",
        "observed_vaf_bin_index",
        "observed_vaf_bin_lower",
        "observed_vaf_bin_upper",
        "observed_vaf_bin",
    ]
    first_cols = [c for c in first_cols if c in out.columns]
    rest = [c for c in out.columns if c not in first_cols]
    return out[first_cols + rest]


def _count_rows(path: Path) -> int | None:
    if not path.exists():
        return None
    opener = gzip.open if path.name.endswith(".gz") else open
    try:
        with opener(path, "rt") as fh:
            return max(0, sum(1 for _ in fh) - 1)
    except Exception:
        return None


def _write_metadata_tsv(path: Path, metadata: dict) -> None:
    rows = []
    for key, value in metadata.items():
        if isinstance(value, (list, tuple)):
            value = ",".join(str(v) for v in value)
        elif isinstance(value, dict):
            value = "; ".join(f"{k}={v}" for k, v in value.items())
        rows.append({"key": key, "value": value})
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def _write_manifest_tsv(path: Path, rows: list[dict]) -> None:
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def run_observation_model_transform(
    *,
    run_or_grid_dir: str | Path,
    scenario_name: str,
    layer_names: str | Iterable[str],
    layer_weights: str | Iterable[float],
    sampling_mode: str,
    phase_mode: str,
    layer_weight_mode: str = "deterministic",
    target_layers: str | Iterable[str] | None = "all",
    layer_mutation_rate_multipliers: str | Iterable[float] | None = None,
    normalize_layer_weights: bool = False,
    out_base_dir: str | Path | None = None,
    overwrite: bool = False,
    observed_nbins: int = 50,
    skip_replicate: bool = False,
    skip_aggregated: bool = False,
    source_vaf_table: str | Path | None = None,
    source_vaf_replicate_table: str | Path | None = None,
) -> dict:
    """
    Apply one observation-model transform to an existing simSOMA run.

    Parameters are explicit function arguments. No observation-model JSON file is
    used or written.
    """
    if observed_nbins <= 0:
        raise ValueError("observed_nbins must be positive")
    if skip_replicate and skip_aggregated:
        raise ValueError("Cannot skip both replicate and aggregated tables")

    model = _validate_model_args(
        scenario_name=scenario_name,
        layer_names=layer_names,
        layer_weights=layer_weights,
        sampling_mode=sampling_mode,
        phase_mode=phase_mode,
        layer_weight_mode=layer_weight_mode,
        target_layers=target_layers,
        layer_mutation_rate_multipliers=layer_mutation_rate_multipliers,
        normalize_layer_weights=normalize_layer_weights,
    )

    run_dir, grid_dir = _resolve_run_and_grid_dir(Path(run_or_grid_dir))
    scenario = model["scenario_name"]
    base_out = Path(out_base_dir).expanduser().resolve() if out_base_dir else run_dir / "observation_model_transforms"
    scenario_dir = base_out / scenario

    if scenario_dir.exists():
        if not overwrite:
            raise FileExistsError(
                f"Observation-transform output folder already exists: {scenario_dir}. "
                "Use overwrite=True or --overwrite to replace it."
            )
        shutil.rmtree(scenario_dir)
    scenario_dir.mkdir(parents=True, exist_ok=False)

    agg_source = _resolve_optional_table(source_vaf_table, run_dir=run_dir, grid_dir=grid_dir) or (grid_dir / STD_AGG_FILE)
    rep_source = _resolve_optional_table(source_vaf_replicate_table, run_dir=run_dir, grid_dir=grid_dir) or (grid_dir / STD_REP_FILE)

    outputs: dict[str, str] = {}
    manifest_rows: list[dict] = []

    if not skip_aggregated:
        agg_df = _read_table(agg_source)
        obs_agg = _transform_table(
            agg_df,
            source_path=agg_source,
            model=model,
            table_kind="aggregated",
            observed_nbins=observed_nbins,
        )
        agg_out = scenario_dir / "observed_vaf_count_spectrum_aggregated_summaries.tsv"
        _write_table(obs_agg, agg_out)
        outputs["observed_vaf_count_spectrum_aggregated_summaries"] = str(agg_out)
        manifest_rows.append({
            "table": "observed_vaf_count_spectrum_aggregated_summaries",
            "source_path": str(agg_source),
            "output_path": str(agg_out),
            "source_rows": len(agg_df),
            "output_rows": len(obs_agg),
        })

    if not skip_replicate:
        rep_df = _read_table(rep_source)
        obs_rep = _transform_table(
            rep_df,
            source_path=rep_source,
            model=model,
            table_kind="replicate",
            observed_nbins=observed_nbins,
        )
        rep_out = scenario_dir / "observed_vaf_count_spectrum_replicate_summaries.tsv.gz"
        _write_table(obs_rep, rep_out)
        outputs["observed_vaf_count_spectrum_replicate_summaries"] = str(rep_out)
        manifest_rows.append({
            "table": "observed_vaf_count_spectrum_replicate_summaries",
            "source_path": str(rep_source),
            "output_path": str(rep_out),
            "source_rows": len(rep_df),
            "output_rows": len(obs_rep),
        })

    metadata = {
        "transform_type": "simSOMA_observation_model_transform",
        "created_utc": _now_iso(),
        "run_dir": str(run_dir),
        "grid_dir": str(grid_dir),
        "scenario_dir": str(scenario_dir),
        "relative_scenario_dir_from_run": _relative_to(scenario_dir, run_dir),
        "scenario_name": scenario,
        "interpretation": (
            "Layer spectra are deterministic copies of the source layer-equivalent simSOMA spectrum. "
            "Copied rows are treated as distinct layer-specific mutation identities under an infinite-sites model."
        ),
        "transform_formula": "observed_vaf = source_vaf * effective_layer_contribution * phase_factor",
        "layer_count": model["layer_count"],
        "layer_names": model["layer_names"],
        "layer_weights": model["layer_weights"],
        "layer_weight_sum": model["layer_weight_sum"],
        "sampling_mode": model["sampling_mode"],
        "phase_mode": model["phase_mode"],
        "phase_factor": model["phase_factor"],
        "layer_weight_mode": model["layer_weight_mode"],
        "target_layers": model["target_layers"],
        "layer_mutation_rate_multipliers": model["layer_mutation_rate_multipliers"],
        "normalize_layer_weights": model["normalize_layer_weights"],
        "observed_nbins": int(observed_nbins),
        "source_vaf_table": str(agg_source) if not skip_aggregated else "SKIPPED",
        "source_vaf_replicate_table": str(rep_source) if not skip_replicate else "SKIPPED",
    }

    metadata_path = scenario_dir / "observation_model_metadata.tsv"
    manifest_path = scenario_dir / "manifest.tsv"
    readme_path = scenario_dir / "README.txt"
    _write_metadata_tsv(metadata_path, metadata)
    _write_manifest_tsv(manifest_path, manifest_rows)
    readme_path.write_text(
        "simSOMA observation-model transform\n"
        "====================================\n\n"
        "This folder contains derived observation-space VAF spectra generated from an existing simSOMA run.\n"
        "No observation-model JSON file was used. All hyperparameters were supplied as function/CLI arguments.\n\n"
        f"Scenario: {scenario}\n"
        f"Formula: {metadata['transform_formula']}\n"
        f"Metadata: {metadata_path.name}\n"
        f"Manifest: {manifest_path.name}\n",
        encoding="utf-8",
    )

    print(f"Run dir:       {run_dir}")
    print(f"Grid dir:      {grid_dir}")
    print(f"Scenario:      {scenario}")
    print(f"Output folder: {scenario_dir}")
    for name, path in outputs.items():
        print(f"{name}: {path}")
    print(f"Metadata:      {metadata_path}")
    print(f"Manifest:      {manifest_path}")

    return {
        "scenario_dir": str(scenario_dir),
        "metadata_path": str(metadata_path),
        "manifest_path": str(manifest_path),
        "outputs": outputs,
        "metadata": metadata,
    }



# -----------------------------------------------------------------------------
# Public mutation-level VAF transformation API
# -----------------------------------------------------------------------------
# This API extends the original aggregated-spectrum post-processing functions
# above.  The original run_observation_model_transform() function and CLI are
# intentionally preserved for backward compatibility.

VAF_TRANSFORMATION_API_VERSION = "0.2.0"
VALID_MUTATION_SAMPLING_MODES = {"bulk", "layer_enriched", "layer_specific"}
VALID_READ_MODELS = {"binomial", "beta_binomial"}


def get_vaf_transformation_capabilities() -> dict[str, Any]:
    """Return the stable capabilities exposed to external clients such as fitSOMA."""
    return {
        "api_version": VAF_TRANSFORMATION_API_VERSION,
        "mutation_level": True,
        "sampling_modes": sorted(VALID_MUTATION_SAMPLING_MODES),
        "phase_modes": sorted(VALID_PHASE_MODES),
        "read_models": sorted(VALID_READ_MODELS),
        "caller_emulation": True,
        "called_any_ascertainment": True,
        "ubiquitous_filtering": False,
        "note": (
            "Ubiquitous-variant filtering and topology-aware summary construction are "
            "analysis operations and remain outside this generative transformation API."
        ),
    }


def _parameter_value(value: Any) -> float:
    """Accept either a scalar or an inference-ready {value: ...} parameter block."""
    if isinstance(value, Mapping):
        if "value" not in value:
            raise ValueError("Parameter dictionaries must contain a 'value' field")
        value = value["value"]
    return float(value)


def _normalize_mutation_observation_config(config: Mapping[str, Any] | None) -> dict[str, Any]:
    cfg = dict(config or {})
    assay = dict(cfg.get("assay", {}))
    assay.setdefault("source_layer", "L2")

    read_model = dict(cfg.get("read_model", {}))
    read_model.setdefault("type", "beta_binomial")
    read_model.setdefault("sequencing_error", 0.001)
    read_model.setdefault("concentration", 200.0)

    depth_model = dict(cfg.get("depth_model", {}))
    depth_model.setdefault("mode", "sample_design")
    depth_model.setdefault("default_mean", 100.0)
    depth_model.setdefault("default_sd", 15.0)
    depth_model.setdefault("minimum", 1)

    caller = dict(cfg.get("caller", {}))
    caller.setdefault("min_depth", 20)
    caller.setdefault("min_alt_reads", 3)
    caller.setdefault("min_observed_vaf", 0.0)

    ascertainment = dict(cfg.get("ascertainment", {}))
    ascertainment.setdefault("retain_called_any", True)

    if read_model["type"] not in VALID_READ_MODELS:
        raise ValueError(f"read_model.type must be one of {sorted(VALID_READ_MODELS)}")
    eps = _parameter_value(read_model["sequencing_error"])
    if not (0.0 <= eps < 0.5):
        raise ValueError("read_model.sequencing_error must be in [0,0.5)")
    concentration = _parameter_value(read_model["concentration"])
    if concentration <= 0:
        raise ValueError("read_model.concentration must be >0")
    if int(depth_model["minimum"]) < 0:
        raise ValueError("depth_model.minimum must be >=0")
    if int(caller["min_depth"]) < 0 or int(caller["min_alt_reads"]) < 0:
        raise ValueError("caller depth/read thresholds must be >=0")
    if not (0.0 <= float(caller["min_observed_vaf"]) <= 1.0):
        raise ValueError("caller.min_observed_vaf must be in [0,1]")

    return {
        "assay": assay,
        "read_model": read_model,
        "depth_model": depth_model,
        "caller": caller,
        "ascertainment": ascertainment,
    }


def _source_layer_weight(row: pd.Series, source_layer: str) -> float | None:
    candidates = [
        f"{source_layer}_weight",
        f"weight_{source_layer}",
        f"{source_layer.lower()}_weight",
        f"weight_{source_layer.lower()}",
    ]
    for col in candidates:
        if col in row.index and pd.notna(row[col]):
            return float(row[col])
    return None


def _prepare_mutation_sample_design(
    sample_design: pd.DataFrame,
    *,
    source_layer: str,
) -> pd.DataFrame:
    required = {"sample_id", "topology_sample_id", "sampling_mode", "phase_mode"}
    missing = sorted(required - set(sample_design.columns))
    if missing:
        raise ValueError(f"sample_design missing columns: {missing}")
    design = sample_design.copy()
    design["sample_id"] = design["sample_id"].astype(str)
    design["topology_sample_id"] = design["topology_sample_id"].astype(str)
    if design["sample_id"].duplicated().any():
        raise ValueError("sample_design.sample_id values must be unique")

    contributions: list[float] = []
    for _, row in design.iterrows():
        sampling = str(row["sampling_mode"]).strip().lower()
        phase = str(row["phase_mode"]).strip().lower()
        if sampling not in VALID_MUTATION_SAMPLING_MODES:
            raise ValueError(
                f"Unsupported sampling_mode {sampling!r}; expected one of "
                f"{sorted(VALID_MUTATION_SAMPLING_MODES)}"
            )
        if phase not in VALID_PHASE_MODES:
            raise ValueError(f"Unsupported phase_mode {phase!r}")

        if "effective_layer_contribution" in design.columns and pd.notna(row.get("effective_layer_contribution")):
            contribution = float(row["effective_layer_contribution"])
        elif sampling == "layer_specific":
            target = str(row.get("target_layer", source_layer) or source_layer)
            contribution = 1.0 if target == source_layer else 0.0
        else:
            weight = _source_layer_weight(row, source_layer)
            if weight is None:
                raise ValueError(
                    f"Sample {row['sample_id']!r} uses {sampling!r} sampling but has no "
                    f"effective_layer_contribution or {source_layer}_weight column"
                )
            contribution = weight
        if not (0.0 <= contribution <= 1.0):
            raise ValueError("effective layer contributions must be in [0,1]")
        contributions.append(contribution)

    design["sampling_mode"] = design["sampling_mode"].astype(str).str.lower()
    design["phase_mode"] = design["phase_mode"].astype(str).str.lower()
    design["effective_layer_contribution"] = contributions
    return design


def _latent_mutation_matrix(
    latent_mutations: pd.DataFrame,
    sample_design: pd.DataFrame,
) -> tuple[list[str], np.ndarray]:
    required = {"mutation_id", "topology_sample_id", "latent_vaf"}
    missing = sorted(required - set(latent_mutations.columns))
    if missing:
        raise ValueError(f"latent_mutations missing columns: {missing}")
    latent = latent_mutations.copy()
    latent["mutation_id"] = latent["mutation_id"].astype(str)
    latent["topology_sample_id"] = latent["topology_sample_id"].astype(str)
    latent["latent_vaf"] = pd.to_numeric(latent["latent_vaf"], errors="raise")
    if not latent["latent_vaf"].between(0.0, 1.0).all():
        raise ValueError("latent_vaf values must be in [0,1]")
    mutation_ids = sorted(latent["mutation_id"].unique())
    topo_samples = sample_design["topology_sample_id"].tolist()
    pivot = latent.pivot_table(
        index="mutation_id",
        columns="topology_sample_id",
        values="latent_vaf",
        aggfunc="max",
        fill_value=0.0,
    ).reindex(index=mutation_ids, columns=topo_samples, fill_value=0.0)
    return mutation_ids, pivot.to_numpy(dtype=float)


def _draw_mutation_depths(
    design: pd.DataFrame,
    config: Mapping[str, Any],
    n_mutations: int,
    rng: np.random.Generator,
) -> np.ndarray:
    dm = config["depth_model"]
    mode = str(dm.get("mode", "sample_design"))
    minimum = int(dm.get("minimum", 1))
    out = np.empty((n_mutations, len(design)), dtype=int)
    for col, row in enumerate(design.itertuples(index=False)):
        row_mean = getattr(row, "depth_mean", np.nan)
        row_sd = getattr(row, "depth_sd", np.nan)
        if mode == "fixed":
            mean = float(dm.get("default_mean", 100.0) if pd.isna(row_mean) else row_mean)
            values = np.repeat(int(round(mean)), n_mutations)
        elif mode == "sample_design":
            mean = float(dm.get("default_mean", 100.0) if pd.isna(row_mean) else row_mean)
            sd = float(dm.get("default_sd", 15.0) if pd.isna(row_sd) else row_sd)
            if sd <= 0:
                values = np.repeat(int(round(mean)), n_mutations)
            else:
                values = np.rint(rng.normal(mean, sd, size=n_mutations)).astype(int)
        else:
            raise ValueError(f"Unsupported depth_model.mode: {mode!r}")
        out[:, col] = np.maximum(values, minimum)
    return out


def transform_vaf_observations(
    latent_mutations: pd.DataFrame,
    sample_design: pd.DataFrame,
    observation_config: Mapping[str, Any] | None = None,
    *,
    seed: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Transform latent mutation VAFs into empirical-style read evidence.

    This is the mutation-level extension of simSOMA's existing VAF observation
    transform.  It combines assay sampling, phasing, depth/read sampling, caller
    emulation, and called-in-any-sample ascertainment.  It deliberately does not
    remove ubiquitous variants or calculate fitSOMA summaries.
    """
    cfg = _normalize_mutation_observation_config(observation_config)
    design = _prepare_mutation_sample_design(
        sample_design, source_layer=str(cfg["assay"]["source_layer"])
    )
    mutation_ids, latent_matrix = _latent_mutation_matrix(latent_mutations, design)
    n_mut, n_samples = latent_matrix.shape
    rng = np.random.default_rng(int(seed))

    contribution = design["effective_layer_contribution"].to_numpy(dtype=float)
    phase = np.asarray([_phase_factor(x) for x in design["phase_mode"]], dtype=float)
    assay_matrix = np.clip(latent_matrix * contribution[None, :] * phase[None, :], 0.0, 1.0)

    eps = _parameter_value(cfg["read_model"]["sequencing_error"])
    read_probability = np.clip(
        assay_matrix * (1.0 - eps) + (1.0 - assay_matrix) * eps, 0.0, 1.0
    )
    depth = _draw_mutation_depths(design, cfg, n_mut, rng)

    read_type = str(cfg["read_model"]["type"])
    if read_type == "binomial":
        alt = rng.binomial(depth, read_probability)
    else:
        concentration = _parameter_value(cfg["read_model"]["concentration"])
        safe_p = np.clip(read_probability, 1e-10, 1.0 - 1e-10)
        locus_p = rng.beta(safe_p * concentration, (1.0 - safe_p) * concentration)
        alt = rng.binomial(depth, locus_p)

    observed_vaf = np.divide(
        alt, depth, out=np.zeros_like(alt, dtype=float), where=depth > 0
    )
    caller = cfg["caller"]
    callable_mask = depth >= int(caller["min_depth"])
    called = (
        callable_mask
        & (alt >= int(caller["min_alt_reads"]))
        & (observed_vaf >= float(caller["min_observed_vaf"]))
    )
    status = np.full(called.shape, "BELOW_CALL_THRESHOLD", dtype=object)
    status[~callable_mask] = "NO_COVERAGE"
    status[called] = "CALLED"

    rows: list[dict[str, Any]] = []
    for i, mutation_id in enumerate(mutation_ids):
        for s, sample in enumerate(design.itertuples(index=False)):
            rows.append({
                "mutation_id": str(mutation_id),
                "sample_id": str(sample.sample_id),
                "topology_sample_id": str(sample.topology_sample_id),
                "sampling_mode": str(sample.sampling_mode),
                "phase_mode": str(sample.phase_mode),
                "effective_layer_contribution": float(sample.effective_layer_contribution),
                "latent_vaf": float(latent_matrix[i, s]),
                "assay_vaf": float(assay_matrix[i, s]),
                "read_probability": float(read_probability[i, s]),
                "alt_count": int(alt[i, s]),
                "ref_count": int(depth[i, s] - alt[i, s]),
                "depth": int(depth[i, s]),
                "observed_vaf": float(observed_vaf[i, s]) if depth[i, s] > 0 else np.nan,
                "callable": bool(callable_mask[i, s]),
                "caller_call": bool(called[i, s]),
                "caller_status": str(status[i, s]),
            })
    evidence = pd.DataFrame(rows)
    n_before = int(evidence["mutation_id"].nunique()) if not evidence.empty else 0
    if bool(cfg["ascertainment"].get("retain_called_any", True)) and not evidence.empty:
        called_any = evidence.groupby("mutation_id")["caller_call"].transform("any")
        evidence = evidence.loc[called_any].reset_index(drop=True)
    n_after = int(evidence["mutation_id"].nunique()) if not evidence.empty else 0

    metadata = {
        "transform_type": "simSOMA_mutation_level_vaf_transformation",
        "api_version": VAF_TRANSFORMATION_API_VERSION,
        "created_utc": _now_iso(),
        "seed": int(seed),
        "source_layer": str(cfg["assay"]["source_layer"]),
        "sampling_modes": sorted(design["sampling_mode"].unique().tolist()),
        "phase_modes": sorted(design["phase_mode"].unique().tolist()),
        "read_model": str(cfg["read_model"]["type"]),
        "sequencing_error": eps,
        "read_concentration": _parameter_value(cfg["read_model"]["concentration"]),
        "depth_model": dict(cfg["depth_model"]),
        "caller": dict(cfg["caller"]),
        "ascertainment": dict(cfg["ascertainment"]),
        "n_samples": int(n_samples),
        "n_latent_mutations": n_before,
        "n_candidates_called_any": n_after,
        "complete_mutation_by_sample_matrix": True,
        "ubiquitous_filter_applied": False,
    }
    return evidence, metadata

def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Transform simSOMA VAF spectra under one explicit observation scenario. No JSON input is used."
    )
    ap.add_argument("--run-dir", required=True, help="simSOMA experiment directory, or its grid_parameter directory.")
    ap.add_argument("--scenario-name", required=True, help="Name for the output scenario folder.")
    ap.add_argument("--layer-names", required=True, help="Comma-separated layer names, e.g. L1,L2,L3.")
    ap.add_argument("--layer-weights", required=True, help="Comma-separated layer weights, e.g. 0.10,0.70,0.20.")
    ap.add_argument("--sampling-mode", required=True, choices=sorted(VALID_SAMPLING_MODES), help="bulk or layer_specific.")
    ap.add_argument("--phase-mode", required=True, choices=sorted(VALID_PHASE_MODES), help="phased or unphased.")
    ap.add_argument("--layer-weight-mode", default="deterministic", choices=sorted(VALID_LAYER_WEIGHT_MODES))
    ap.add_argument("--target-layers", default="all", help="For layer_specific sampling: all or comma-separated target layers, e.g. L2.")
    ap.add_argument("--layer-mutation-rate-multipliers", default="", help="Optional comma-separated count multipliers per layer.")
    ap.add_argument("--normalize-layer-weights", action="store_true", help="Normalize supplied layer weights to sum to 1.")
    ap.add_argument("--out-base-dir", default="", help="Optional base output directory. Default: <run-dir>/observation_model_transforms.")
    ap.add_argument("--observed-nbins", type=int, default=50, help="Number of observed VAF bins. Default: 50.")
    ap.add_argument("--overwrite", action="store_true", help="Replace an existing scenario output folder.")
    ap.add_argument("--skip-replicate", action="store_true", help=f"Do not transform {STD_REP_FILE}.")
    ap.add_argument("--skip-aggregated", action="store_true", help=f"Do not transform {STD_AGG_FILE}.")
    ap.add_argument("--source-vaf-table", default="", help="Optional explicit aggregated source VAF table.")
    ap.add_argument("--source-vaf-replicate-table", default="", help="Optional explicit replicate source VAF table.")
    args = ap.parse_args(argv)

    run_observation_model_transform(
        run_or_grid_dir=args.run_dir,
        scenario_name=args.scenario_name,
        layer_names=args.layer_names,
        layer_weights=args.layer_weights,
        sampling_mode=args.sampling_mode,
        phase_mode=args.phase_mode,
        layer_weight_mode=args.layer_weight_mode,
        target_layers=args.target_layers,
        layer_mutation_rate_multipliers=args.layer_mutation_rate_multipliers or None,
        normalize_layer_weights=bool(args.normalize_layer_weights),
        out_base_dir=args.out_base_dir or None,
        overwrite=bool(args.overwrite),
        observed_nbins=int(args.observed_nbins),
        skip_replicate=bool(args.skip_replicate),
        skip_aggregated=bool(args.skip_aggregated),
        source_vaf_table=args.source_vaf_table or None,
        source_vaf_replicate_table=args.source_vaf_replicate_table or None,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
