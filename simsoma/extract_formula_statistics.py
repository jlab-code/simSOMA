#!/usr/bin/env python3
"""
Extract manuscript-style formula statistics from simSOMA outputs.

This is an official post-processing tool. It does not run simulations and does
not modify primary simSOMA outputs. Point it at one experiment folder and it will
look for standardized primary and observation-transformed VAF spectrum tables.

Detected input tables
---------------------
Primary simSOMA output:
  <experiment>/grid_parameter/vaf_count_spectrum_aggregated_summaries.csv

Observation-model transformed output(s):
  <experiment>/observation_model_transforms/<scenario>/observed_vaf_count_spectrum_aggregated_summaries.tsv

Output folders
--------------
For each detected source table, this tool overwrites/writes:
  <source-folder>/formula_statistics_tables/

Statistics
----------
For each organ and parameter set, compute:
  fixed_fraction
  intermediate_fraction
  private_fraction
  normalized_sharedness

The across-organ summary is grouped only by set_id plus data-space labels.
Combinatorial input parameters are displayed by merging grid_parameter/parameter_sets.csv.
Derived organ descriptors such as n_bottlenecks_root_to_organ are never treated
as parameter columns in the across-organ summary.

For primary outputs, the formula VAF is sampled_vaf.
For observation-transformed outputs, the formula VAF is scenario-normalized:
  formula_vaf = observed_vaf / expected_fixed_observed_vaf
where expected_fixed_observed_vaf is detected from the transformed table or
computed as effective_layer_contribution * phase_factor.
"""

from __future__ import annotations

import argparse
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

PRIMARY_VAF_FILE = "vaf_count_spectrum_aggregated_summaries.csv"
OBSERVED_VAF_FILE = "observed_vaf_count_spectrum_aggregated_summaries.tsv"
FORMULA_DIR_NAME = "formula_statistics_tables"

PARAM_ALIASES = {
    "bottlenecks": "n_bottlenecks_root_to_organ",
    "m": "m",
    "rho": "rho",
    "branch_precursor_number": "branch_precursor_number_realized_cells",
    "branch_precursor_number_realized_cells": "branch_precursor_number_realized_cells",
    "organ_precursor_number": "organ_precursor_number_realized_cells",
    "organ_precursor_number_realized_cells": "organ_precursor_number_realized_cells",
    "sam_boundary_cells": "sam_boundary_cells",
}

# Columns that may appear in the VAF table and are useful in the per-organ
# diagnostic table. These are not automatically treated as simulation input
# parameters in the across-organ summary.
ORGAN_DESCRIPTOR_COLUMNS = [
    "n_bottlenecks_root_to_organ",
]

# Columns that should never be copied from parameter_sets.csv into the
# parameter-set-level summary. The summary should display only true
# combinatorial input parameters from the config/grid.
NON_CONFIG_PARAMETER_COLUMNS = {
    "organ_id",
    "organ",
    "organ_name",
    "sampled_organ",
    "branch_id",
    "branch",
    "n_bottlenecks_root_to_organ",
    "n_bottlenecks",
    "root_to_organ_length",
    "root_to_organ_distance",
    "root_to_organ_age",
    "developmental_depth",
    "tree_depth",
    "topological_depth",
}

# Backward-compatible name used internally for optional per-organ descriptors.
FIGURE_PARAMETER_COLUMNS = ORGAN_DESCRIPTOR_COLUMNS

STAT_COLUMNS = [
    "fixed_fraction",
    "intermediate_fraction",
    "private_fraction",
    "normalized_sharedness",
]


@dataclass(frozen=True)
class SourceTable:
    data_space: str  # primary_simSOMA or observation_model_transform
    source_folder: Path
    table_path: Path
    scenario_name: str
    sep: str


def _num(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def _fmt(v) -> str:
    if pd.isna(v):
        return "NA"
    try:
        f = float(v)
        if abs(f - round(f)) < 1e-9:
            return str(int(round(f)))
        return f"{f:g}"
    except Exception:
        return str(v)


def _sort_key(x):
    try:
        return (0, float(x))
    except Exception:
        return (1, str(x))


def _parse_list(s: str | None) -> list[str]:
    if s is None or str(s).strip() == "":
        return []
    return [x.strip() for x in str(s).split(",") if x.strip()]


def _resolve_col(name: str, cols: Iterable[str]) -> str:
    c = PARAM_ALIASES.get(name, name)
    if c not in set(cols):
        raise ValueError(f"Parameter column '{name}' -> '{c}' is not available in the input table.")
    return c


def apply_fixed_filters(df: pd.DataFrame, fixed: list[str]) -> pd.DataFrame:
    """Apply optional plotter-style filters, e.g. --fixed rho=0 --fixed m=min,max."""
    out = df.copy()

    for spec in fixed:
        if "=" not in spec:
            raise ValueError(f"Invalid --fixed specification '{spec}'. Use column=value or alias=value.")

        key, value_string = spec.split("=", 1)
        col = _resolve_col(key.strip(), out.columns)
        tokens = _parse_list(value_string)
        vals_avail = sorted(out[col].dropna().unique(), key=_sort_key)
        if not vals_avail:
            raise ValueError(f"Column {col} has no non-missing values.")

        selected = []
        for tok in tokens:
            tl = tok.lower()
            if tl == "min":
                selected.append(vals_avail[0])
                continue
            if tl == "max":
                selected.append(vals_avail[-1])
                continue

            match = None
            for a in vals_avail:
                if str(a) == tok:
                    match = a
                    break
                try:
                    if float(a) == float(tok):
                        match = a
                        break
                except Exception:
                    pass
            if match is None:
                raise ValueError(
                    f"Value '{tok}' not found for {col}. Available: {','.join(_fmt(x) for x in vals_avail)}"
                )
            selected.append(match)

        out = out[out[col].isin(selected)].copy()

    if out.empty:
        raise ValueError("No rows remain after applying fixed filters.")
    return out


def discover_sources(result_dir: Path, target: str = "all") -> list[SourceTable]:
    """Discover primary and observation-transformed spectrum tables."""
    root = result_dir.expanduser().resolve()
    sources: list[SourceTable] = []

    if target not in {"all", "primary", "observed"}:
        raise ValueError("target must be one of: all, primary, observed")

    # Direct grid directory.
    if target in {"all", "primary"} and (root / PRIMARY_VAF_FILE).exists():
        sources.append(SourceTable("primary_simSOMA", root, root / PRIMARY_VAF_FILE, "primary", ","))

    # Experiment directory with grid_parameter.
    if target in {"all", "primary"} and (root / "grid_parameter" / PRIMARY_VAF_FILE).exists():
        sources.append(SourceTable("primary_simSOMA", root / "grid_parameter", root / "grid_parameter" / PRIMARY_VAF_FILE, "primary", ","))

    # Direct observation scenario directory.
    if target in {"all", "observed"} and (root / OBSERVED_VAF_FILE).exists():
        sources.append(SourceTable("observation_model_transform", root, root / OBSERVED_VAF_FILE, root.name, "\t"))

    # Experiment directory with observation_model_transforms/*.
    obs_root = root / "observation_model_transforms"
    if target in {"all", "observed"} and obs_root.exists():
        for p in sorted(obs_root.glob(f"*/{OBSERVED_VAF_FILE}")):
            sources.append(SourceTable("observation_model_transform", p.parent, p, p.parent.name, "\t"))

    # De-duplicate by path while preserving order.
    seen = set()
    unique: list[SourceTable] = []
    for s in sources:
        key = s.table_path.resolve()
        if key not in seen:
            unique.append(s)
            seen.add(key)

    return unique


def load_and_prepare_vaf(source: SourceTable) -> pd.DataFrame:
    df = pd.read_csv(source.table_path, sep=source.sep)

    if source.data_space == "primary_simSOMA":
        if "sampled_vaf" not in df.columns:
            raise ValueError(f"Primary input table is missing sampled_vaf: {source.table_path}")
        df["formula_vaf"] = _num(df["sampled_vaf"])
        df["expected_fixed_observed_vaf"] = 1.0
        df["formula_vaf_definition"] = "sampled_vaf"
    else:
        if "relative_observed_vaf" in df.columns:
            df["formula_vaf"] = _num(df["relative_observed_vaf"])
            df["formula_vaf_definition"] = "relative_observed_vaf"
            if "expected_fixed_observed_vaf" not in df.columns:
                if {"effective_layer_contribution", "phase_factor"}.issubset(df.columns):
                    df["expected_fixed_observed_vaf"] = (
                        _num(df["effective_layer_contribution"]) * _num(df["phase_factor"])
                    )
        else:
            if "observed_vaf" not in df.columns:
                raise ValueError(f"Observed input table is missing observed_vaf: {source.table_path}")
            if "expected_fixed_observed_vaf" in df.columns:
                fixed_vaf = _num(df["expected_fixed_observed_vaf"])
                df["formula_vaf_definition"] = "observed_vaf / expected_fixed_observed_vaf"
            elif {"effective_layer_contribution", "phase_factor"}.issubset(df.columns):
                fixed_vaf = _num(df["effective_layer_contribution"]) * _num(df["phase_factor"])
                df["expected_fixed_observed_vaf"] = fixed_vaf
                df["formula_vaf_definition"] = "observed_vaf / (effective_layer_contribution * phase_factor)"
            else:
                raise ValueError(
                    "Observed input table must contain either relative_observed_vaf, "
                    "expected_fixed_observed_vaf, or both effective_layer_contribution and phase_factor. "
                    f"Missing from: {source.table_path}"
                )
            obs = _num(df["observed_vaf"])
            df["formula_vaf"] = np.where(fixed_vaf > 0, obs / fixed_vaf, np.nan)

    for c in ["formula_vaf", "n_variants_mean", "sharing_degree", "expected_fixed_observed_vaf"]:
        if c in df.columns:
            df[c] = _num(df[c])

    df["data_space"] = source.data_space
    df["observation_scenario_name"] = source.scenario_name
    return df


def organ_all(vaf: pd.DataFrame) -> pd.DataFrame:
    required = {"level", "variant_class", "formula_vaf", "n_variants_mean", "set_id", "organ_id"}
    missing = sorted(required - set(vaf.columns))
    if missing:
        raise ValueError(f"Input VAF table is missing required columns: {missing}")

    d = vaf[(vaf["level"].astype(str) == "organ") & (vaf["variant_class"].astype(str) == "all")].copy()
    if d.empty:
        raise ValueError("No organ-level variant_class='all' rows found in VAF table.")

    d["formula_vaf"] = _num(d["formula_vaf"])
    d["n_variants_mean"] = _num(d["n_variants_mean"])
    d = d.dropna(subset=["formula_vaf", "n_variants_mean"])
    return d


def _param_cols(d: pd.DataFrame) -> list[str]:
    return [c for c in FIGURE_PARAMETER_COLUMNS if c in d.columns]


def _infer_experiment_dir(source: SourceTable) -> Path:
    """Infer the experiment directory that contains grid_parameter/."""
    folder = source.source_folder.resolve()

    # Primary source may be either the experiment folder or the grid_parameter folder.
    if folder.name == "grid_parameter":
        return folder.parent
    if (folder / "grid_parameter").exists():
        return folder

    # Observation transform scenario folder: <experiment>/observation_model_transforms/<scenario>/
    parts = list(folder.parts)
    if "observation_model_transforms" in parts:
        idx = parts.index("observation_model_transforms")
        if idx > 0:
            return Path(*parts[:idx])

    # Fallback: source folder itself.
    return folder


def load_parameter_sets_for_source(source: SourceTable) -> pd.DataFrame:
    """Load the authoritative parameter_sets.csv for the experiment, if available.

    These columns define the combinatorial parameter setting for each set_id. They
    are merged into the parameter-set summary so the summary is directly usable
    for plotting/filtering.
    """
    exp_dir = _infer_experiment_dir(source)
    candidates = [
        exp_dir / "grid_parameter" / "parameter_sets.csv",
        source.source_folder / "parameter_sets.csv",
    ]
    for path in candidates:
        if path.exists():
            d = pd.read_csv(path)
            if "set_id" in d.columns:
                return d.drop_duplicates(subset=["set_id"]).copy()
    return pd.DataFrame()


def compute_fixed_intermediate(vaf: pd.DataFrame, low_threshold: float) -> pd.DataFrame:
    d = organ_all(vaf)
    param_cols = _param_cols(d)
    for c in param_cols:
        d[c] = _num(d[c])

    group_cols = ["set_id", "organ_id"] + param_cols
    rows = []

    for keys, g in d.groupby(group_cols, dropna=False):
        total = float(g["n_variants_mean"].sum())
        if total <= 0 or not np.isfinite(total):
            continue

        v = g["formula_vaf"].to_numpy(float)
        n = g["n_variants_mean"].to_numpy(float)

        fixed_mask = np.isclose(v, 1.0, rtol=1e-7, atol=1e-9)
        intermediate_mask = (v >= low_threshold) & (v < 1.0) & (~fixed_mask)

        fixed_n = float(n[fixed_mask].sum())
        intermediate_n = float(n[intermediate_mask].sum())

        row = dict(zip(group_cols, keys if isinstance(keys, tuple) else (keys,)))
        row.update({
            "total_observed_variants_mean_vaf": total,
            "fixed_variants_mean": fixed_n,
            "fixed_fraction": fixed_n / total,
            "intermediate_variants_mean": intermediate_n,
            "intermediate_fraction": intermediate_n / total,
            "low_threshold": low_threshold,
        })
        rows.append(row)

    return pd.DataFrame(rows)


def compute_private_sharedness(vaf: pd.DataFrame) -> pd.DataFrame:
    required = {"level", "variant_class", "sharing_degree", "n_variants_mean", "set_id", "organ_id"}
    missing = sorted(required - set(vaf.columns))
    if missing:
        raise ValueError(f"Input VAF table is missing required columns: {missing}")

    d = vaf[(vaf["level"].astype(str) == "organ") & (vaf["variant_class"].astype(str) == "degree")].copy()
    if d.empty:
        raise ValueError(
            "No organ-level variant_class='degree' rows found. Private fraction and normalized sharedness cannot be computed exactly."
        )

    param_cols = _param_cols(d)
    for c in ["sharing_degree", "n_variants_mean"] + param_cols:
        d[c] = _num(d[c])
    d = d.dropna(subset=["sharing_degree", "n_variants_mean"])

    norg = organ_all(vaf).groupby("set_id", dropna=False)["organ_id"].nunique().rename("n_organs").reset_index()

    group_cols = ["set_id", "organ_id"] + param_cols
    rows = []

    for keys, g in d.groupby(group_cols, dropna=False):
        total = float(g["n_variants_mean"].sum())
        if total <= 0 or not np.isfinite(total):
            continue

        private_n = float(g.loc[np.isclose(g["sharing_degree"], 1.0), "n_variants_mean"].sum())
        mean_degree = float((g["sharing_degree"] * g["n_variants_mean"]).sum() / total)

        row = dict(zip(group_cols, keys if isinstance(keys, tuple) else (keys,)))
        row.update({
            "total_observed_variants_mean_degree": total,
            "private_variants_mean": private_n,
            "private_fraction": private_n / total,
            "mean_sharing_degree": mean_degree,
        })
        rows.append(row)

    out = pd.DataFrame(rows)
    if out.empty:
        return out

    out = out.merge(norg, on="set_id", how="left")
    denom = out["n_organs"] - 1
    out["normalized_sharedness"] = np.where(
        denom > 0,
        ((out["mean_sharing_degree"] - 1) / denom).clip(0, 1),
        np.nan,
    )
    return out


def merge_statistics(vaf_stats: pd.DataFrame, sharing_stats: pd.DataFrame, source: SourceTable, formula_definition: str) -> pd.DataFrame:
    key_cols = [
        c for c in ["set_id", "organ_id"] + FIGURE_PARAMETER_COLUMNS
        if c in vaf_stats.columns and c in sharing_stats.columns
    ]
    if not key_cols:
        raise ValueError("Could not identify shared key columns for merging VAF and sharing statistics.")

    out = vaf_stats.merge(sharing_stats, on=key_cols, how="outer")
    out.insert(0, "data_space", source.data_space)
    out.insert(1, "observation_scenario_name", source.scenario_name)
    out.insert(2, "formula_vaf_definition", formula_definition)

    first_cols = [
        c for c in [
            "data_space", "observation_scenario_name", "formula_vaf_definition",
            "set_id", "organ_id",
        ] + FIGURE_PARAMETER_COLUMNS if c in out.columns
    ]
    formula_cols = [
        "total_observed_variants_mean_vaf",
        "fixed_variants_mean",
        "fixed_fraction",
        "intermediate_variants_mean",
        "intermediate_fraction",
        "low_threshold",
        "total_observed_variants_mean_degree",
        "private_variants_mean",
        "private_fraction",
        "mean_sharing_degree",
        "n_organs",
        "normalized_sharedness",
    ]
    formula_cols = [c for c in formula_cols if c in out.columns]
    rest = [c for c in out.columns if c not in first_cols + formula_cols]

    return out[first_cols + formula_cols + rest].sort_values(first_cols).reset_index(drop=True)


def _organ_varying_columns(organ_stats: pd.DataFrame) -> list[str]:
    """Return known/candidate columns that vary across organs within set_id."""
    base_cols = [c for c in ["data_space", "observation_scenario_name", "set_id"] if c in organ_stats.columns]
    candidate_cols = [c for c in FIGURE_PARAMETER_COLUMNS if c in organ_stats.columns]
    varying: list[str] = []
    if not base_cols:
        return candidate_cols
    for col in candidate_cols:
        nunique = organ_stats.groupby(base_cols, dropna=False)[col].nunique(dropna=False)
        if len(nunique) > 0 and (nunique > 1).any():
            varying.append(col)
    return varying


def _prepare_parameter_set_columns(parameter_sets: pd.DataFrame, set_ids: pd.Series) -> pd.DataFrame:
    """Prepare true config/grid input parameter columns for summary output.

    The only authoritative source for combinatorial parameter settings is
    grid_parameter/parameter_sets.csv. Columns derived from organ placement or
    downstream summaries are excluded even if they accidentally appear in that
    file. If parameter_sets.csv is unavailable, return set_id only rather than
    inferring parameter columns from organ-level VAF tables.
    """
    set_id_values = sorted(pd.Series(set_ids).dropna().unique())
    if parameter_sets.empty or "set_id" not in parameter_sets.columns:
        return pd.DataFrame({"set_id": set_id_values})

    keep = parameter_sets.copy()
    keep = keep[keep["set_id"].isin(set_id_values)].copy()
    if keep.empty:
        return pd.DataFrame({"set_id": set_id_values})

    excluded = []
    allowed_cols = []
    for c in keep.columns:
        if c == "set_id":
            allowed_cols.append(c)
            continue
        if c in NON_CONFIG_PARAMETER_COLUMNS:
            excluded.append(c)
            continue
        if keep[c].isna().all():
            excluded.append(c)
            continue
        allowed_cols.append(c)

    keep = keep[allowed_cols].drop_duplicates(subset=["set_id"]).copy()
    keep.attrs["excluded_non_config_parameter_columns"] = ",".join(excluded)
    return keep

def make_parameter_set_summary(organ_stats: pd.DataFrame, parameter_sets: pd.DataFrame | None = None) -> pd.DataFrame:
    """
    Build the figure-ready parameter-set summary from the organ-level statistics.

    This function deliberately uses organ_level_formula_statistics as the only
    statistical input. It melts the four formula statistics into long format,
    groups only by data_space, observation_scenario_name, set_id, and statistic,
    and then summarizes across organ_id. The combinatorial parameter setting is
    added only after this aggregation by merging grid_parameter/parameter_sets.csv
    on set_id.

    Therefore, organ-derived descriptors such as n_bottlenecks_root_to_organ can
    remain in organ_level_formula_statistics.tsv for diagnostics, but they cannot
    split the across-organ summary table.
    """
    required = {"set_id", "organ_id"}
    missing = sorted(required - set(organ_stats.columns))
    if missing:
        raise ValueError(f"Cannot build parameter-set summary; organ_stats is missing columns: {missing}")

    id_cols = [c for c in ["data_space", "observation_scenario_name", "set_id", "organ_id"] if c in organ_stats.columns]
    value_cols = [c for c in STAT_COLUMNS if c in organ_stats.columns]
    if not value_cols:
        return pd.DataFrame()

    # Explicitly remove duplicate organ-stat rows before averaging. Duplicate rows
    # for the same set_id x organ_id x statistic would otherwise overweight one organ.
    long = organ_stats[id_cols + value_cols].melt(
        id_vars=id_cols,
        value_vars=value_cols,
        var_name="statistic",
        value_name="organ_stat_value",
    )
    long["organ_stat_value"] = pd.to_numeric(long["organ_stat_value"], errors="coerce")
    long = long.dropna(subset=["organ_stat_value"]).copy()
    if long.empty:
        return pd.DataFrame()

    dedup_cols = [c for c in ["data_space", "observation_scenario_name", "set_id", "organ_id", "statistic"] if c in long.columns]
    before = len(long)
    long = long.drop_duplicates(subset=dedup_cols, keep="first").copy()
    duplicate_organ_stat_rows_removed = before - len(long)

    group_cols = [c for c in ["data_space", "observation_scenario_name", "set_id", "statistic"] if c in long.columns]
    out = (
        long.groupby(group_cols, dropna=False)["organ_stat_value"]
        .agg(
            organ_mean="mean",
            organ_sd="std",
            organ_min="min",
            organ_median="median",
            organ_max="max",
            n_organs_with_stat="count",
        )
        .reset_index()
    )

    # Track organ-derived columns that exist in the organ-level table but are
    # intentionally excluded from this summary.
    excluded_cols = sorted([
        c for c in organ_stats.columns
        if c not in set(id_cols + value_cols)
        and (c in NON_CONFIG_PARAMETER_COLUMNS or c in ORGAN_DESCRIPTOR_COLUMNS or "organ" in c.lower() or "bottleneck" in c.lower())
    ])
    if "n_bottlenecks_root_to_organ" in organ_stats.columns and "n_bottlenecks_root_to_organ" not in excluded_cols:
        excluded_cols.append("n_bottlenecks_root_to_organ")

    # Add the true config/grid parameter combination after the across-organ
    # aggregation. This is the only source of parameter columns in the summary.
    param = _prepare_parameter_set_columns(parameter_sets if parameter_sets is not None else pd.DataFrame(), out["set_id"])
    excluded_non_config = param.attrs.get("excluded_non_config_parameter_columns", "") if hasattr(param, "attrs") else ""
    if not param.empty:
        out = out.merge(param, on="set_id", how="left")

    out["summary_level"] = "parameter_set_across_organs"
    out["variability_unit"] = "sampled_organs"
    out["summary_source_table"] = "organ_level_formula_statistics.tsv"
    out["summary_grouping_definition"] = "grouped_by_data_space_observation_scenario_name_set_id_statistic_then_averaged_across_organ_id"
    out["parameter_columns_source"] = "grid_parameter/parameter_sets.csv" if parameter_sets is not None and not parameter_sets.empty else "set_id_only_no_parameter_sets_csv_found"
    out["parameter_columns_definition"] = "true_config_grid_input_parameters_merged_by_set_id_after_across_organ_aggregation"
    out["excluded_organ_varying_columns"] = ",".join(sorted(set(excluded_cols)))
    out["excluded_non_config_parameter_columns"] = excluded_non_config
    out["duplicate_organ_stat_rows_removed"] = duplicate_organ_stat_rows_removed

    # Put the full combinatorial parameter setting next to set_id and before statistic.
    id_cols_out = [c for c in ["data_space", "observation_scenario_name", "set_id"] if c in out.columns]
    param_cols = [c for c in param.columns if c != "set_id" and c in out.columns]
    descriptor_cols = [
        "statistic",
        "summary_level",
        "variability_unit",
        "summary_source_table",
        "summary_grouping_definition",
        "parameter_columns_source",
        "parameter_columns_definition",
        "excluded_organ_varying_columns",
        "excluded_non_config_parameter_columns",
        "duplicate_organ_stat_rows_removed",
    ]
    stat_cols = [
        "organ_mean",
        "organ_sd",
        "organ_min",
        "organ_median",
        "organ_max",
        "n_organs_with_stat",
    ]
    first_cols = [c for c in id_cols_out + param_cols + descriptor_cols if c in out.columns]
    stat_cols = [c for c in stat_cols if c in out.columns]
    rest = [c for c in out.columns if c not in first_cols + stat_cols]
    out = out[first_cols + stat_cols + rest]
    sort_cols = [c for c in ["data_space", "observation_scenario_name", "set_id", "statistic"] if c in out.columns]
    return out.sort_values(sort_cols).reset_index(drop=True)


def write_definitions(out_dir: Path, low_threshold: float, data_space: str, formula_definition: str) -> None:
    if data_space == "primary_simSOMA":
        fixed_formula = "Fix_o = |{i in M_o: sampled_vaf_io = 1}| / |M_o|"
        int_formula = f"Int_o = |{{i in M_o: {low_threshold:g} <= sampled_vaf_io < 1}}| / |M_o|"
    else:
        fixed_formula = "Fix_o = |{i in M_o: v_formula_io = 1}| / |M_o|, where v_formula = observed_vaf / expected_fixed_observed_vaf"
        int_formula = f"Int_o = |{{i in M_o: {low_threshold:g} <= v_formula_io < 1}}| / |M_o|"

    definitions = pd.DataFrame([
        {
            "statistic": "fixed_fraction",
            "formula": fixed_formula,
            "implementation": "fixed_variants_mean / total_observed_variants_mean_vaf; counts are n_variants_mean from organ-level variant_class='all' rows",
            "formula_vaf_definition": formula_definition,
        },
        {
            "statistic": "intermediate_fraction",
            "formula": int_formula,
            "implementation": "intermediate_variants_mean / total_observed_variants_mean_vaf; counts are n_variants_mean from organ-level variant_class='all' rows",
            "formula_vaf_definition": formula_definition,
        },
        {
            "statistic": "private_fraction",
            "formula": "Priv_o = |{i in M_o: D_i = 1}| / |M_o|",
            "implementation": "private_variants_mean / total_observed_variants_mean_degree; counts are n_variants_mean from organ-level variant_class='degree' rows with sharing_degree == 1",
            "formula_vaf_definition": "not VAF-based",
        },
        {
            "statistic": "normalized_sharedness",
            "formula": "Shr_o = mean_{i in M_o}((D_i - 1) / (N_organs - 1))",
            "implementation": "weighted mean of (sharing_degree - 1)/(n_organs - 1), using n_variants_mean from organ-level variant_class='degree' rows; equivalent to (mean_sharing_degree - 1)/(n_organs - 1)",
            "formula_vaf_definition": "not VAF-based",
        },
    ])
    definitions.to_csv(out_dir / "statistic_definitions.tsv", sep="\t", index=False)


def extract_one_source(
    source: SourceTable,
    *,
    low_threshold: float,
    fixed_filters: list[str],
    overwrite: bool = True,
) -> dict:
    if not source.table_path.exists():
        raise FileNotFoundError(source.table_path)

    out_dir = source.source_folder / FORMULA_DIR_NAME
    if out_dir.exists() and overwrite:
        shutil.rmtree(out_dir)
    elif out_dir.exists() and not overwrite:
        raise FileExistsError(f"Output folder exists: {out_dir}")
    out_dir.mkdir(parents=True, exist_ok=True)

    vaf = load_and_prepare_vaf(source)
    formula_def = str(vaf["formula_vaf_definition"].dropna().iloc[0]) if "formula_vaf_definition" in vaf.columns and not vaf.empty else "unknown"

    if fixed_filters:
        vaf = apply_fixed_filters(vaf, fixed_filters)

    vaf_stats = compute_fixed_intermediate(vaf, low_threshold)
    sharing_stats = compute_private_sharedness(vaf)
    organ_stats = merge_statistics(vaf_stats, sharing_stats, source, formula_def)

    organ_stats.to_csv(out_dir / "organ_level_formula_statistics.tsv", sep="\t", index=False)

    parameter_sets = load_parameter_sets_for_source(source)
    parameter_summary = make_parameter_set_summary(organ_stats, parameter_sets=parameter_sets)
    if not parameter_summary.empty:
        parameter_summary.to_csv(out_dir / "parameter_set_formula_statistics_summary.tsv", sep="\t", index=False)

    write_definitions(out_dir, low_threshold, source.data_space, formula_def)

    summary = pd.DataFrame([
        {
            "key": "data_space",
            "value": source.data_space,
        },
        {
            "key": "observation_scenario_name",
            "value": source.scenario_name,
        },
        {
            "key": "source_folder",
            "value": str(source.source_folder),
        },
        {
            "key": "source_table",
            "value": str(source.table_path),
        },
        {
            "key": "output_folder",
            "value": str(out_dir),
        },
        {
            "key": "low_threshold",
            "value": low_threshold,
        },
        {
            "key": "fixed_filters",
            "value": ";".join(fixed_filters),
        },
        {
            "key": "formula_vaf_definition",
            "value": formula_def,
        },
        {
            "key": "n_organ_stat_rows",
            "value": len(organ_stats),
        },
        {
            "key": "organ_level_table_definition",
            "value": "one row per set_id x organ_id; includes organ-varying covariates such as n_bottlenecks_root_to_organ when present",
        },
        {
            "key": "parameter_set_summary_definition",
            "value": "one row per set_id x statistic; built from organ_level_formula_statistics.tsv by averaging across organ_id, then true config/grid input parameters are merged from parameter_sets.csv",
        },
        {
            "key": "parameter_columns_source",
            "value": "grid_parameter/parameter_sets.csv, merged by set_id after across-organ aggregation; no organ-derived columns are used for grouping",
        },
        {
            "key": "outputs",
            "value": "organ_level_formula_statistics.tsv;parameter_set_formula_statistics_summary.tsv;statistic_definitions.tsv;extraction_run_summary.tsv",
        },
    ])
    summary.to_csv(out_dir / "extraction_run_summary.tsv", sep="\t", index=False)

    return {
        "data_space": source.data_space,
        "observation_scenario_name": source.scenario_name,
        "source_table": str(source.table_path),
        "output_folder": str(out_dir),
        "n_organ_stat_rows": len(organ_stats),
        "formula_vaf_definition": formula_def,
    }


def run_formula_statistics_extraction(
    *,
    result_dir: str | Path,
    target: str = "all",
    low_threshold: float = 0.05,
    fixed_filters: list[str] | None = None,
    overwrite: bool = True,
) -> list[dict]:
    fixed_filters = fixed_filters or []
    sources = discover_sources(Path(result_dir), target=target)
    if not sources:
        raise FileNotFoundError(
            f"No standardized VAF summary tables found under {Path(result_dir).expanduser().resolve()}. "
            f"Expected {PRIMARY_VAF_FILE} or {OBSERVED_VAF_FILE}."
        )

    rows = []
    for source in sources:
        rows.append(
            extract_one_source(
                source,
                low_threshold=low_threshold,
                fixed_filters=fixed_filters,
                overwrite=overwrite,
            )
        )

    # If result_dir is an experiment folder, write a compact manifest there.
    root = Path(result_dir).expanduser().resolve()
    try:
        pd.DataFrame(rows).to_csv(root / "formula_statistics_extraction_manifest.tsv", sep="\t", index=False)
    except Exception:
        pass
    return rows


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=(
            "Extract fixed, intermediate, private, and normalized sharedness statistics from "
            "primary simSOMA output and/or observation-model transformed output."
        )
    )
    ap.add_argument("--result-dir", required=True, help="Experiment folder, grid_parameter folder, or observation scenario folder.")
    ap.add_argument("--target", default="all", choices=["all", "primary", "observed"], help="Which standardized outputs to analyze. Default: all.")
    ap.add_argument("--low-threshold", type=float, default=0.05, help="Lower formula-VAF threshold for intermediate variants. Default: 0.05.")
    ap.add_argument("--fixed", action="append", default=[], help="Optional filter before extraction. Example: --fixed rho=0 --fixed m=min,max")
    ap.add_argument("--no-overwrite", action="store_true", help="Do not overwrite existing formula_statistics_tables folders.")
    args = ap.parse_args(argv)

    rows = run_formula_statistics_extraction(
        result_dir=args.result_dir,
        target=args.target,
        low_threshold=args.low_threshold,
        fixed_filters=args.fixed,
        overwrite=not args.no_overwrite,
    )

    print(f"Discovered and processed {len(rows)} source table(s).")
    for r in rows:
        print(f"- {r['data_space']} | {r['observation_scenario_name']} | rows={r['n_organ_stat_rows']}")
        print(f"  source: {r['source_table']}")
        print(f"  output: {r['output_folder']}")
        print(f"  VAF formula: {r['formula_vaf_definition']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
