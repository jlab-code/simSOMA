from __future__ import annotations

import argparse
import copy
import csv
import itertools
import gzip
import gc
import json
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import summaries

BIAS_FIELDS = [
    "bias_mode",
    "branch_bias_value",
    "branch_bias_mean",
    "branch_bias_kappa",
]

CONFIG_PARAMS = [
    "m",
    "rho",
    "mu_year",
    "victim_locality",
    *BIAS_FIELDS,
    "sam_boundary_cells",
    "branch_precursor_number",
    "organ_precursor_number",
    "organ_total_cells",
    "seq_fraction",
]

SELF_RENEWAL_FIELDS = [
    "m",
    "rho",
    "victim_locality",
    "mu_year",
    *BIAS_FIELDS,
]
PRE_BRANCHING_FIELDS = [
    "sam_boundary_cells",
]
BRANCHING_FIELDS = [
    "branch_precursor_number",
]
ORGAN_FIELDS = [
    "organ_precursor_number",
    "organ_total_cells",
    "seq_fraction",
]
DERIVED_OUTPUT_FIELDS = [
    "branch_precursor_number_realized_cells",
    "organ_precursor_number_realized_cells",
    "organ_total_cells_realized",
    "sequenced_cells_used",
]
OUTPUT_PARAM_FIELDS = SELF_RENEWAL_FIELDS + PRE_BRANCHING_FIELDS + BRANCHING_FIELDS + ORGAN_FIELDS + DERIVED_OUTPUT_FIELDS
TOPOLOGY_OUTPUT_FIELDS = [
    "mapping_unit",
    "mapping_mode",
    "mapping_rate",
    "effective_kappa_sr",
]

MODULE_PARAM_FIELDS = {
    "self_renewal": SELF_RENEWAL_FIELDS,
    "pre_branching": PRE_BRANCHING_FIELDS,
    "branching": BRANCHING_FIELDS,
    "organ": ORGAN_FIELDS,
}

DEFAULT_TOPOLOGY_PLOT = {
    # Single cleaned topology-plotting backend used for GitHub release.
    # It draws an axis-preserving developmental topology and prunes branches
    # that do not lead to sampled terminal organs.
    "enabled": True,
    "style": "developmental_pruned",
    "fig_width": 4.2,
    "fig_height": 2.6,
    "dpi": 300,
    "y_ticks": "auto",
    "y_label": "auto",
    "tick_size": 8,
    "label_size": 9,
    "branch_linewidth": 1.0,
    "connector_linewidth": 0.8,
    "branchpoint_size": 9,
    "organ_size": 16,
    "tip_size": 6,
    "min_branch_offset": 0.35,
    "show_nonterminal_organs": False,
    "organ_stub_length": 0.28,
    "nonterminal_organ_style": "dotted",
    "nonterminal_organ_size": 13,
    "branch_color_mode": "order",
    "write_pdf": True,
    "write_layout_csv": True,
    "title": "",
}


def _resolve_path(value: str, *, base_dir: Path) -> Path:
    p = Path(value).expanduser()
    if not p.is_absolute():
        p = (base_dir / p).resolve()
    return p


def _project_root_from_context(*, base_dir: Path) -> Path:
    env_home = os.environ.get("SIMSOMA_HOME")
    if env_home:
        return Path(env_home).expanduser().resolve()

    # Standard package layout: config files are in simSOMA_configs/.
    if base_dir.name == "simSOMA_configs" and (base_dir.parent / "simSOMA_corefunc").exists():
        return base_dir.parent.resolve()

    # Split configs are written below simSOMA_output/<experiment>/grid_splits/subconfigs/.
    # Walk upward until the project root is found.
    for parent in [base_dir, *base_dir.parents]:
        if (parent / "simSOMA_corefunc").exists() and (parent / "simSOMA_configs").exists():
            return parent.resolve()

    # Last-resort fallback keeps legacy behavior predictable.
    return base_dir.resolve()


def _resolve_project_prefixed_path(value: str, *, base_dir: Path, prefix: str, env_var: str | None = None) -> Path:
    p = Path(value).expanduser()
    if p.is_absolute():
        return p.resolve()

    if env_var and os.environ.get(env_var):
        env_root = Path(os.environ[env_var]).expanduser().resolve()
        if p.parts and p.parts[0] == prefix:
            return (env_root / Path(*p.parts[1:])).resolve()

    if p.parts and p.parts[0] == prefix:
        project_root = _project_root_from_context(base_dir=base_dir)
        return (project_root / p).resolve()

    return (base_dir / p).resolve()


def _resolve_outdir_root(value: str, *, base_dir: Path) -> Path:
    return _resolve_project_prefixed_path(
        value,
        base_dir=base_dir,
        prefix="simSOMA_output",
        env_var="SIMSOMA_OUTPUT_DIR",
    )


def _resolve_topology_json_path(value: str, *, base_dir: Path) -> Path:
    return _resolve_project_prefixed_path(
        value,
        base_dir=base_dir,
        prefix="simSOMA_inputs",
        env_var="SIMSOMA_INPUT_DIR",
    )


def _relative_path_str(target: Path, *, start: Path) -> str:
    try:
        return Path(os.path.relpath(str(target), str(start))).as_posix()
    except Exception:
        return str(target)


def _path_str_for_output(value: str, *, config_base_dir: Path, output_base_dir: Path) -> str:
    resolved = _resolve_topology_json_path(value, base_dir=config_base_dir)
    return _relative_path_str(resolved, start=output_base_dir)


def _json_dump(obj: Any, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, sort_keys=True), encoding="utf-8")


def _sha256_file(path: Path) -> str:
    h = __import__("hashlib").sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _software_version_info(pipeline_dir: Path) -> Dict[str, Any]:
    project_root = pipeline_dir.parent
    version_file = project_root / "VERSION"
    version = version_file.read_text(encoding="utf-8").strip() if version_file.exists() else "unversioned"
    tracked = [
        "run_from_config.py",
        "pipeline_wrapper.py",
        "self_renewal.py",
        "pre_branching.py",
        "branching.py",
        "organ.py",
        "summaries.py",
        "branch_bias.py",
        "topology_io.py",
        "inspect_topology.py",
        "launch_grid_splits.py",
    ]
    files: Dict[str, Any] = {}
    for name in tracked:
        fp = pipeline_dir / name
        files[name] = {
            "present": fp.exists(),
            "sha256": _sha256_file(fp) if fp.exists() else None,
            "size_bytes": fp.stat().st_size if fp.exists() else None,
        }
    return {
        "software_version": version,
        "project_root": str(project_root),
        "pipeline_dir": str(pipeline_dir),
        "tracked_files": files,
    }


def _write_software_version(run_dir: Path, pipeline_dir: Path) -> Path:
    path = run_dir / "software_version.json"
    _json_dump(_software_version_info(pipeline_dir), path)
    return path


def _write_csv(rows: List[Dict[str, Any]], path: Path, *, fieldnames: Optional[List[str]] = None) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    keys = fieldnames if fieldnames is not None else sorted({k for r in rows for k in r.keys()})
    if path.suffix == ".gz":
        with gzip.open(path, "wt", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            for r in rows:
                w.writerow({k: r.get(k, "") for k in keys})
        return
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in keys})



class _StreamingCsvWriter:
    """Lazy CSV/gzip writer with a fixed, pre-declared schema.

    This avoids holding all output rows in memory while still keeping CSV
    headers stable. For known output kinds, the schema includes optional columns
    that may appear only in later rows, such as ``sharing_degree`` or VAF-bin
    columns. Unexpected extra columns still raise an error because widening a
    CSV schema after writing the header would silently corrupt the file.
    """

    def __init__(self, path: Path, *, kind: str, fieldnames: Optional[List[str]] = None) -> None:
        self.path = path
        self.kind = kind
        self.fieldnames: List[str] | None = list(fieldnames) if fieldnames is not None else None
        self._fh: Any = None
        self._writer: csv.DictWriter | None = None
        self.n_rows = 0

    def _open(self, row: Dict[str, Any]) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        if self.fieldnames is None:
            self.fieldnames = _ordered_fieldnames(self.kind, [row], include_all_ordered=True)
        else:
            schema_row = {k: "" for k in self.fieldnames}
            schema_row.update(row)
            self.fieldnames = _ordered_fieldnames(self.kind, [schema_row], include_all_ordered=True)
        if self.path.suffix == ".gz":
            self._fh = gzip.open(self.path, "wt", newline="", encoding="utf-8")
        else:
            self._fh = self.path.open("w", newline="", encoding="utf-8")
        self._writer = csv.DictWriter(self._fh, fieldnames=self.fieldnames)
        self._writer.writeheader()

    def write_row(self, row: Dict[str, Any]) -> None:
        if self._writer is None:
            self._open(row)
        assert self.fieldnames is not None
        extra = set(row.keys()).difference(self.fieldnames)
        if extra:
            raise ValueError(
                f"Cannot write row with new column(s) after header was written for {self.path}: {sorted(extra)}"
            )
        assert self._writer is not None
        self._writer.writerow({k: row.get(k, "") for k in self.fieldnames})
        self.n_rows += 1

    def write_rows(self, rows: Iterable[Dict[str, Any]]) -> None:
        for row in rows:
            self.write_row(row)

    def close(self) -> None:
        if self._fh is not None:
            self._fh.close()
            self._fh = None
            self._writer = None


def _close_streaming_writers(writers: Iterable[_StreamingCsvWriter]) -> None:
    for writer in writers:
        writer.close()

def _natural_key(text: str) -> List[Any]:
    parts = re.split(r"(\d+)", text)
    out: List[Any] = []
    for part in parts:
        if part.isdigit():
            out.append(int(part))
        else:
            out.append(part)
    return out


def _numeric_bin_columns(keys: Iterable[str]) -> List[str]:
    bins = [k for k in keys if re.fullmatch(r"bin_\d+", k)]
    return sorted(bins, key=lambda x: int(x.split("_", 1)[1]))


def _topology_row(cfg: Dict[str, Any]) -> Dict[str, Any]:
    return {
        "mapping_unit": cfg.get("mapping_unit"),
        "mapping_mode": cfg.get("mapping_mode"),
        "mapping_rate": _effective_kappa_sr(cfg),
        "effective_kappa_sr": _effective_kappa_sr(cfg),
    }


def _load_topology_json(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _organ_topology_metadata(topo: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    branches = topo.get("branches", {})
    if not isinstance(branches, dict):
        return {}

    depth_by_branch: Dict[str, int] = {}

    def branch_depth(bid: str) -> int:
        if bid in depth_by_branch:
            return depth_by_branch[bid]
        b = branches.get(bid, {})
        parent_id = b.get("parent_id")
        if parent_id is None:
            depth_by_branch[bid] = 0
        else:
            depth_by_branch[bid] = branch_depth(str(parent_id)) + 1
        return depth_by_branch[bid]

    organ_meta: Dict[str, Dict[str, Any]] = {}
    for bid in sorted(branches.keys(), key=_natural_key):
        depth = branch_depth(str(bid))
        events = branches[bid].get("events", [])
        if not isinstance(events, list):
            continue
        for ev in events:
            if str(ev.get("type", "")).upper() != "ORGAN":
                continue
            organ_id = str(ev.get("target_id", ev.get("target", "")))
            if not organ_id:
                continue
            organ_meta[organ_id] = {
                "n_bottlenecks_root_to_organ": int(depth) + 1,
            }
    return organ_meta


def _ordered_fieldnames(kind: str, rows: List[Dict[str, Any]], *, include_all_ordered: bool = False) -> List[str]:
    keys = {k for r in rows for k in r.keys()}
    bins = _numeric_bin_columns(keys)
    remainder = keys.difference(bins)

    if kind == "parameter_sets":
        ordered: List[str] = [
            "set_id",
            "set_label",
            "n_sim",
            *TOPOLOGY_OUTPUT_FIELDS,
            *SELF_RENEWAL_FIELDS,
            *PRE_BRANCHING_FIELDS,
            *BRANCHING_FIELDS,
            *ORGAN_FIELDS,
            *DERIVED_OUTPUT_FIELDS,
        ]
    elif kind == "replicate_summaries":
        ordered = [
            "set_id",
            "set_label",
            "rep",
            "seed",
            *TOPOLOGY_OUTPUT_FIELDS,
            *SELF_RENEWAL_FIELDS,
            *PRE_BRANCHING_FIELDS,
            *BRANCHING_FIELDS,
            *ORGAN_FIELDS,
            *DERIVED_OUTPUT_FIELDS,
            "total_variants",
            "total_unique_variants",
            "total_variant_occurrences",
            "n_shared_mutations",
            "n_private_mutations_total",
        ]
    elif kind == "aggregated_summaries":
        ordered = [
            "set_id",
            "set_label",
            "n_replicates",
            *TOPOLOGY_OUTPUT_FIELDS,
            *SELF_RENEWAL_FIELDS,
            *PRE_BRANCHING_FIELDS,
            *BRANCHING_FIELDS,
            *ORGAN_FIELDS,
            *DERIVED_OUTPUT_FIELDS,
            "private_shared_available",
            "private_shared_reason",
            "total_variants_mean",
            "total_variants_sd",
            "total_unique_variants_mean",
            "total_unique_variants_sd",
            "total_variant_occurrences_mean",
            "total_variant_occurrences_sd",
            "n_shared_mutations_mean",
            "n_shared_mutations_sd",
            "n_private_mutations_total_mean",
            "n_private_mutations_total_sd",
        ]
    elif kind == "organ_replicate_summaries":
        ordered = [
            "set_id",
            "set_label",
            "rep",
            "seed",
            "organ_id",
            "n_bottlenecks_root_to_organ",
            "sequenced_cells",
            *TOPOLOGY_OUTPUT_FIELDS,
            *SELF_RENEWAL_FIELDS,
            *PRE_BRANCHING_FIELDS,
            *BRANCHING_FIELDS,
            *ORGAN_FIELDS,
            *DERIVED_OUTPUT_FIELDS,
            "total_variants",
            "singletons",
            "fixed",
            "n_lineages_sampled",
            "dominant_lineage_id",
            "dominant_fraction",
            "lineage_entropy",
            "boundary_cells",
            "organ_total_cells",
            "organ_precursor_number",
        ]
    elif kind == "organ_aggregated_summaries":
        ordered = [
            "set_id",
            "set_label",
            "n_replicates",
            "organ_id",
            "n_bottlenecks_root_to_organ",
            "sequenced_cells",
            *TOPOLOGY_OUTPUT_FIELDS,
            *SELF_RENEWAL_FIELDS,
            *PRE_BRANCHING_FIELDS,
            *BRANCHING_FIELDS,
            *ORGAN_FIELDS,
            *DERIVED_OUTPUT_FIELDS,
            "total_variants_mean",
            "total_variants_sd",
            "singletons_mean",
            "singletons_sd",
            "fixed_mean",
            "fixed_sd",
            "n_lineages_sampled_mean",
            "n_lineages_sampled_sd",
            "dominant_lineage_id_mean",
            "dominant_lineage_id_sd",
            "dominant_fraction_mean",
            "dominant_fraction_sd",
            "lineage_entropy_mean",
            "lineage_entropy_sd",
            "boundary_cells_mean",
            "boundary_cells_sd",
            "organ_total_cells_mean",
            "organ_total_cells_sd",
            "organ_precursor_number_mean",
            "organ_precursor_number_sd",
        ]
    elif kind == "vaf_class_summaries":
        ordered = [
            "set_id",
            "set_label",
            "rep",
            "seed",
            "n_replicates",
            "level",
            "organ_id",
            "n_bottlenecks_root_to_organ",
            "variant_class",
            "bin_index",
            "vaf_left",
            "vaf_right",
            *TOPOLOGY_OUTPUT_FIELDS,
            *SELF_RENEWAL_FIELDS,
            *PRE_BRANCHING_FIELDS,
            *BRANCHING_FIELDS,
            *ORGAN_FIELDS,
            *DERIVED_OUTPUT_FIELDS,
            "n_variants",
            "n_variants_mean",
            "n_variants_sd",
        ]
    elif kind == "sharing_summaries":
        ordered = [
            "set_id",
            "set_label",
            "rep",
            "seed",
            "n_replicates",
            "statistic",
            "variant_class",
            "sharing_degree",
            "organ_id",
            "organ_a",
            "organ_b",
            *TOPOLOGY_OUTPUT_FIELDS,
            *SELF_RENEWAL_FIELDS,
            *PRE_BRANCHING_FIELDS,
            *BRANCHING_FIELDS,
            *ORGAN_FIELDS,
            *DERIVED_OUTPUT_FIELDS,
            "n_variants",
            "n_variants_mean",
            "n_variants_sd",
        ]
    elif kind == "vaf_count_spectrum_summaries":
        ordered = [
            "set_id",
            "set_label",
            "rep",
            "seed",
            "n_replicates",
            "level",
            "organ_id",
            "n_bottlenecks_root_to_organ",
            "variant_class",
            "sharing_degree",
            "n_sampled_cells",
            "allele_count",
            "sampled_vaf",
            *TOPOLOGY_OUTPUT_FIELDS,
            *SELF_RENEWAL_FIELDS,
            *PRE_BRANCHING_FIELDS,
            *BRANCHING_FIELDS,
            *ORGAN_FIELDS,
            *DERIVED_OUTPUT_FIELDS,
            "n_variants",
            "n_variants_mean",
            "n_variants_sd",
        ]
    elif kind == "raw_vafs":
        ordered = [
            "set_id",
            "set_label",
            "rep",
            "seed",
            "organ_id",
            "n_bottlenecks_root_to_organ",
            "seq_fraction",
            "organ_total_cells",
            "sequenced_cells_used",
            "mutation_id",
            "allele_count",
            "sequenced_cells",
            "vaf",
        ]
    else:
        ordered = []

    fieldnames: List[str] = []
    for k in ordered:
        if include_all_ordered or k in remainder:
            if k not in fieldnames:
                fieldnames.append(k)
            remainder.discard(k)

    other_cols = sorted(remainder, key=_natural_key)
    fieldnames.extend(other_cols)
    fieldnames.extend(bins)
    return fieldnames


def _streaming_schema_hint_row(
    kind: str,
    *,
    cfg: Optional[Dict[str, Any]] = None,
    organ_meta: Optional[Dict[str, Dict[str, Any]]] = None,
) -> Dict[str, Any]:
    """Return optional columns that should be present in streaming CSV headers.

    Legacy non-streaming output derived schemas from the union of all rows. A
    streaming writer must instead declare columns before seeing later rows. This
    helper adds predictable optional columns that can appear after the first row
    for a file type.
    """
    hint: Dict[str, Any] = {}
    nbins = int((cfg or {}).get("vaf_nbins", 0) or 0)

    if kind in {"replicate_summaries", "organ_replicate_summaries"}:
        for i in range(nbins):
            hint[f"bin_{i}"] = ""

    if kind == "organ_aggregated_summaries":
        for i in range(nbins):
            hint[f"bin_{i}_mean"] = ""
            hint[f"bin_{i}_sd"] = ""

    if kind == "aggregated_summaries" and organ_meta:
        for oid in sorted(organ_meta.keys(), key=_natural_key):
            hint[f"{oid}_private_mean"] = ""
            hint[f"{oid}_private_sd"] = ""

    return hint


def _streaming_fieldnames(
    kind: str,
    *,
    cfg: Optional[Dict[str, Any]] = None,
    organ_meta: Optional[Dict[str, Dict[str, Any]]] = None,
) -> List[str]:
    return _ordered_fieldnames(
        kind,
        [_streaming_schema_hint_row(kind, cfg=cfg, organ_meta=organ_meta)],
        include_all_ordered=True,
    )


def _load_config(config_path: Path) -> Dict[str, Any]:
    cfg = json.loads(config_path.read_text(encoding="utf-8"))
    if not isinstance(cfg, dict):
        raise ValueError("Config must be a JSON object.")
    return cfg


def _is_grouped_config(cfg: Dict[str, Any]) -> bool:
    return any(k in cfg for k in ["run", "topology", "check", "simulation", "plot"])




def _normalize_config(cfg: Dict[str, Any]) -> Dict[str, Any]:
    if not _is_grouped_config(cfg):
        raise ValueError(
            "Legacy flat configs are no longer supported. Use the grouped config layout with run/topology/check/simulation blocks."
        )

    run = cfg.get("run", {})
    topology = cfg.get("topology", {})
    check = cfg.get("check", {})
    simulation = cfg.get("simulation", {})

    return {
        "experiment_name": run.get("experiment_name"),
        "outdir_root": run.get("outdir_root"),
        "seed": run.get("seed", 1),
        "topology_json": topology.get("topology_json"),
        "mapping_unit": topology.get("mapping_unit", "steps"),
        "mapping_mode": topology.get("mapping_mode", "deterministic"),
        "mapping_rate": topology.get("mapping_rate"),
        "phyllotaxy": topology.get("phyllotaxy"),
        "topology_plot": check.get("topology_plot", {}),
        "simulation_mode": simulation.get("mode", "grid_parameter"),
        "n_sim": simulation.get("n_sim", 1),
        "summaries": simulation.get("summaries", True),
        "vaf_nbins": simulation.get("vaf_nbins", 20),
        "summary_private_shared": simulation.get("summary_private_shared", True),
        "store_full_results": simulation.get("store_full_results", True),
        "export_raw_vafs": simulation.get("export_raw_vafs", False),
        # Current compact output profile:
        # - sharing summaries are retained by default because they are central diagnostics;
        # - exact allele-count VAF spectra are retained by default because they are re-binnable;
        # - legacy binned VAF-class summaries are opt-in only.
        "export_sharing_summaries": simulation.get("export_sharing_summaries", simulation.get("export_summary_feature_tables", True)),
        "export_replicate_sharing_summaries": simulation.get("export_replicate_sharing_summaries", simulation.get("export_replicate_feature_tables", False)),
        "export_vaf_count_spectra": simulation.get("export_vaf_count_spectra", True),
        "export_legacy_binned_vaf_summaries": simulation.get("export_legacy_binned_vaf_summaries", False),
        "export_legacy_binned_vaf_replicates": simulation.get("export_legacy_binned_vaf_replicates", False),
        "modules": simulation.get("modules"),
        "grid_subset": simulation.get("grid_subset"),
        "replicate_subset": simulation.get("replicate_subset"),
    }

def _grouped_config_uses_simulation_kappa(raw_cfg: Dict[str, Any]) -> bool:
    if not _is_grouped_config(raw_cfg):
        return False
    modules = raw_cfg.get("simulation", {}).get("modules")
    if not isinstance(modules, dict):
        return False
    sr = modules.get("self_renewal")
    return isinstance(sr, dict) and "kappa_sr" in sr


def _effective_kappa_sr(cfg: Dict[str, Any]) -> float:
    mapping_unit = str(cfg["mapping_unit"]).lower()
    rate = cfg.get("mapping_rate")
    if mapping_unit == "steps":
        return 1.0 if rate is None else float(rate)
    if rate is None:
        raise ValueError("topology.mapping_rate must be provided when topology.mapping_unit is years or meters.")
    return float(rate)



def _coerce_requested_cell_count(value: Any, *, name: str) -> int:
    try:
        fv = float(value)
    except Exception as exc:
        raise ValueError(f"{name} must be numeric.") from exc
    if not math.isfinite(fv):
        raise ValueError(f"{name} must be finite.")
    iv = int(round(fv))
    if abs(fv - iv) > 1e-9:
        raise ValueError(f"{name} must be an integer cell count.")
    if iv < 1:
        raise ValueError(f"{name} must be >= 1.")
    return iv


def _realize_branch_precursor_count(requested: Any, *, m: Any, name: str) -> Dict[str, Any]:
    req = _coerce_requested_cell_count(requested, name=name)
    m_int = int(m)
    if m_int < 1:
        raise ValueError("m must be >= 1")
    realized = min(req, m_int)
    return {
        "requested_cells": int(req),
        "realized_cells": int(realized),
        "capped": bool(realized < req),
        "realized_fraction_of_m": float(realized) / float(m_int),
    }


def _realize_organ_precursor_count(requested: Any, *, sam_boundary_cells: Any, organ_total_cells: Any, name: str) -> Dict[str, Any]:
    req = _coerce_requested_cell_count(requested, name=name)
    B = int(sam_boundary_cells)
    O = int(organ_total_cells)
    if B < 1:
        raise ValueError("sam_boundary_cells must be >= 1")
    if O < 1:
        raise ValueError("organ_total_cells must be >= 1")
    realized = min(req, B, O)
    return {
        "requested_cells": int(req),
        "realized_cells": int(realized),
        "capped": bool(realized < req),
        "realized_fraction_of_boundary": float(realized) / float(B),
    }


def _derive_n_cells_and_sequenced_cells(*, organ_total_cells: Any, seq_fraction: Any) -> Tuple[int, int]:
    O = int(organ_total_cells)
    frac = float(seq_fraction)
    if O <= 0:
        raise ValueError("organ_total_cells must be >= 1")
    if not math.isfinite(frac) or frac <= 0.0 or frac > 1.0:
        raise ValueError("seq_fraction must satisfy 0 < seq_fraction <= 1")
    sequenced_cells = int(math.ceil(frac * O))
    sequenced_cells = max(1, min(sequenced_cells, O))
    return O, sequenced_cells


def _params_with_effective_kappa(cfg: Dict[str, Any], params: Dict[str, Any]) -> Dict[str, Any]:
    out = dict(params)
    out["kappa_sr"] = _effective_kappa_sr(cfg)
    m_int = int(out["m"])
    B_int = int(out["sam_boundary_cells"])
    O_int = int(out["organ_total_cells"])
    if B_int < m_int:
        raise ValueError(f"sam_boundary_cells ({B_int}) must be >= m ({m_int}).")

    branch_rec = _realize_branch_precursor_count(out["branch_precursor_number"], m=m_int, name="branch_precursor_number")
    organ_rec = _realize_organ_precursor_count(
        out["organ_precursor_number"],
        sam_boundary_cells=B_int,
        organ_total_cells=O_int,
        name="organ_precursor_number",
    )
    n_cells, sequenced_cells = _derive_n_cells_and_sequenced_cells(organ_total_cells=O_int, seq_fraction=out["seq_fraction"])

    out["branch_precursor_number_realized_cells"] = int(branch_rec["realized_cells"])
    out["organ_precursor_number_realized_cells"] = int(organ_rec["realized_cells"])
    out["organ_total_cells_realized"] = int(n_cells)
    out["sequenced_cells_used"] = int(sequenced_cells)
    return out

def _resolve_param_spec(spec: Any) -> List[Any]:
    if isinstance(spec, dict):
        if "values" in spec:
            values = spec["values"]
            if not isinstance(values, list) or len(values) == 0:
                raise ValueError("Parameter spec with 'values' must contain a non-empty list.")
            return list(values)

        spacing = spec.get("spacing")
        if spacing is None:
            raise ValueError(
                "Parameter spec must be either a fixed scalar value, a {'values': [...]} block, "
                "or a spacing block with spacing/start/stop/num."
            )

        if "start" not in spec or "stop" not in spec or "num" not in spec:
            raise ValueError("Spacing-based parameter specs require 'start', 'stop', and 'num'.")

        start = spec["start"]
        stop = spec["stop"]
        num = int(spec["num"])
        if num < 1:
            raise ValueError("Parameter spec 'num' must be >= 1.")

        spacing = str(spacing).lower()
        if spacing == "linear":
            if num == 1:
                return [start]
            start_f = float(start)
            stop_f = float(stop)
            step = (stop_f - start_f) / float(num - 1)
            return [start_f + i * step for i in range(num)]

        if spacing == "log":
            start_f = float(start)
            stop_f = float(stop)
            if start_f <= 0 or stop_f <= 0:
                raise ValueError("Log-spaced parameter specs require start > 0 and stop > 0.")
            if num == 1:
                return [start_f]
            log_start = math.log(start_f)
            log_stop = math.log(stop_f)
            step = (log_stop - log_start) / float(num - 1)
            return [math.exp(log_start + i * step) for i in range(num)]

        if spacing == "integer":
            if num == 1:
                return [int(round(float(start)))]
            start_f = float(start)
            stop_f = float(stop)
            step = (stop_f - start_f) / float(num - 1)
            vals: List[int] = []
            for i in range(num):
                v = int(round(start_f + i * step))
                if v not in vals:
                    vals.append(v)
            return vals

        raise ValueError("Unsupported spacing. Use one of: linear, log, integer.")

    if isinstance(spec, list):
        if len(spec) == 0:
            raise ValueError("Parameter spec lists must be non-empty.")
        return list(spec)

    return [spec]


def _normalize_module_param_locations(modules: Dict[str, Any]) -> Dict[str, Any]:
    """Return a copy of simulation.modules using the canonical parameter layout.

    Canonical layout places branch-specific self-renewal bias parameters under
    simulation.modules.self_renewal.  For backward compatibility, old configs
    that still place these fields under simulation.modules.branching are accepted
    and migrated in-memory before validation.
    """
    if not isinstance(modules, dict):
        return modules

    out = copy.deepcopy(modules)
    sr = out.get("self_renewal")
    br = out.get("branching")
    if not isinstance(sr, dict) or not isinstance(br, dict):
        return out

    for field in BIAS_FIELDS:
        in_sr = field in sr
        in_br = field in br
        if in_br and not in_sr:
            sr[field] = br.pop(field)
        elif in_br and in_sr:
            if br[field] != sr[field]:
                raise ValueError(
                    f"Parameter {field!r} is present in both simulation.modules.self_renewal "
                    "and simulation.modules.branching with different values. "
                    "Keep it only under self_renewal."
                )
            br.pop(field)
    return out


def _extract_module_param_specs(modules: Dict[str, Any]) -> Dict[str, Any]:
    modules = _normalize_module_param_locations(modules)
    if not isinstance(modules, dict) or len(modules) == 0:
        raise ValueError("simulation.modules must be a non-empty JSON object.")

    unknown_modules = sorted(set(modules.keys()).difference(MODULE_PARAM_FIELDS.keys()), key=_natural_key)
    if unknown_modules:
        raise ValueError(
            f"Unknown simulation.modules sections: {unknown_modules}. "
            f"Allowed sections are: {sorted(MODULE_PARAM_FIELDS.keys())}."
        )

    missing_modules = [m for m in MODULE_PARAM_FIELDS if m not in modules]
    if missing_modules:
        raise ValueError(f"Missing required simulation.modules sections: {missing_modules}")

    out: Dict[str, Any] = {}
    for module_name, allowed_params in MODULE_PARAM_FIELDS.items():
        module_params = modules[module_name]
        if not isinstance(module_params, dict):
            raise ValueError(f"simulation.modules.{module_name} must be a JSON object.")
        extra_params = sorted(set(module_params.keys()).difference(allowed_params), key=_natural_key)
        if extra_params:
            raise ValueError(
                f"Unknown parameters in simulation.modules.{module_name}: {extra_params}. "
                f"Allowed parameters are: {list(allowed_params)}."
            )
        missing_params = [p for p in allowed_params if p not in module_params]
        if missing_params:
            raise ValueError(
                f"Missing required parameters in simulation.modules.{module_name}: {missing_params}"
            )
        for param_name in allowed_params:
            out[param_name] = module_params[param_name]
    return out



def _grid_axes_from_modules(modules: Dict[str, Any]) -> Tuple[List[str], Dict[str, List[Any]]]:
    specs = _extract_module_param_specs(modules)
    axes: Dict[str, List[Any]] = {}
    param_names: List[str] = []
    for k in CONFIG_PARAMS:
        axes[k] = _resolve_param_spec(specs[k])
        param_names.append(k)
    return param_names, axes



def _normalize_grid_subset(grid_subset: Any, *, total_sets: int) -> Optional[Dict[str, int]]:
    if grid_subset is None:
        return None
    if not isinstance(grid_subset, dict):
        raise ValueError("simulation.grid_subset must be a JSON object if provided.")
    if total_sets < 1:
        raise ValueError("simulation.grid_subset requires at least one parameter set.")

    start = grid_subset.get("start_index", grid_subset.get("start", 1))
    end = grid_subset.get("end_index", grid_subset.get("end", total_sets))
    try:
        start_i = int(start)
        end_i = int(end)
    except Exception as exc:
        raise ValueError("simulation.grid_subset.start_index and end_index must be integers.") from exc
    if start_i < 1 or end_i < 1:
        raise ValueError("simulation.grid_subset indices must be >= 1.")
    if end_i < start_i:
        raise ValueError("simulation.grid_subset.end_index must be >= start_index.")
    if start_i > total_sets:
        raise ValueError(
            f"simulation.grid_subset.start_index={start_i} exceeds the full grid size ({total_sets})."
        )
    end_i = min(end_i, total_sets)
    return {"start_index": int(start_i), "end_index": int(end_i)}


def _normalize_replicate_subset(replicate_subset: Any, *, total_replicates: int) -> Optional[Dict[str, int]]:
    if replicate_subset is None:
        return None
    if not isinstance(replicate_subset, dict):
        raise ValueError("simulation.replicate_subset must be a JSON object if provided.")
    if total_replicates < 1:
        raise ValueError("simulation.replicate_subset requires simulation.n_sim >= 1.")

    start = replicate_subset.get("start_index", replicate_subset.get("start_rep", replicate_subset.get("start", 1)))
    end = replicate_subset.get("end_index", replicate_subset.get("end_rep", replicate_subset.get("end", total_replicates)))
    try:
        start_i = int(start)
        end_i = int(end)
    except Exception as exc:
        raise ValueError("simulation.replicate_subset.start_index and end_index must be integers.") from exc
    if start_i < 1 or end_i < 1:
        raise ValueError("simulation.replicate_subset indices must be >= 1.")
    if end_i < start_i:
        raise ValueError("simulation.replicate_subset.end_index must be >= start_index.")
    if start_i > total_replicates:
        raise ValueError(
            f"simulation.replicate_subset.start_index={start_i} exceeds simulation.n_sim ({total_replicates})."
        )
    end_i = min(end_i, total_replicates)
    return {"start_index": int(start_i), "end_index": int(end_i)}


def _iter_grid_parameter_sets_with_indices(
    modules: Dict[str, Any],
    *,
    grid_subset: Any = None,
) -> Iterable[Tuple[int, Dict[str, Any]]]:
    total_sets = _grid_total_n_sets(modules)
    subset = _normalize_grid_subset(grid_subset, total_sets=total_sets)
    start_i = 1 if subset is None else int(subset["start_index"])
    end_i = total_sets if subset is None else int(subset["end_index"])
    for idx, params in enumerate(_iter_grid_parameter_sets(modules), start=1):
        if idx < start_i:
            continue
        if idx > end_i:
            break
        yield idx, params


def _grid_subset_n_sets(modules: Dict[str, Any], *, grid_subset: Any = None) -> int:
    total_sets = _grid_total_n_sets(modules)
    subset = _normalize_grid_subset(grid_subset, total_sets=total_sets)
    if subset is None:
        return total_sets
    return int(subset["end_index"]) - int(subset["start_index"]) + 1



def _validate_config(cfg: Dict[str, Any]) -> None:
    required = [
        "experiment_name",
        "topology_json",
        "outdir_root",
        "mapping_unit",
        "mapping_mode",
        "n_sim",
        "seed",
        "summaries",
        "vaf_nbins",
        "summary_private_shared",
        "simulation_mode",
    ]
    missing = [k for k in required if k not in cfg or cfg[k] is None]
    if missing:
        raise ValueError(f"Missing required config keys after normalization: {missing}")

    if int(cfg["n_sim"]) < 1:
        raise ValueError("simulation.n_sim must be >= 1")

    mapping_unit = str(cfg["mapping_unit"]).lower()
    if mapping_unit not in {"steps", "years", "meters"}:
        raise ValueError("topology.mapping_unit must be one of: steps, years, meters")

    if str(cfg["mapping_mode"]).lower() not in {"deterministic", "poisson"}:
        raise ValueError("topology.mapping_mode must be one of: deterministic, poisson")

    if mapping_unit != "steps" and cfg.get("mapping_rate") is None:
        raise ValueError(
            "For topology.mapping_unit = years or meters, provide topology.mapping_rate explicitly. "
            "This value defines the one-time conversion from observed topology units to internal self-renewal steps."
        )

    plot_cfg = cfg.get("topology_plot", {})
    if plot_cfg is not None and not isinstance(plot_cfg, dict):
        raise ValueError("check.topology_plot must be a JSON object if provided.")
    if isinstance(plot_cfg, dict):
        plot_style = str(plot_cfg.get("style", DEFAULT_TOPOLOGY_PLOT["style"])).lower()
        allowed_plot_styles = {"developmental_pruned"}
        if plot_style not in allowed_plot_styles:
            raise ValueError(
                "check.topology_plot.style must be 'developmental_pruned' in the cleaned GitHub release. "
                "Older topology plot styles were removed to keep the package minimal."
            )

        if str(plot_cfg.get("nonterminal_organ_style", DEFAULT_TOPOLOGY_PLOT["nonterminal_organ_style"])).lower() not in {"dotted", "dashed", "solid"}:
            raise ValueError("check.topology_plot.nonterminal_organ_style must be one of: dotted, dashed, solid")
        if str(plot_cfg.get("branch_color_mode", DEFAULT_TOPOLOGY_PLOT["branch_color_mode"])).lower() not in {"order", "black"}:
            raise ValueError("check.topology_plot.branch_color_mode must be one of: order, black")

    phyllotaxy_cfg = cfg.get("phyllotaxy")
    if phyllotaxy_cfg is not None and not isinstance(phyllotaxy_cfg, dict):
        raise ValueError("topology.phyllotaxy must be a JSON object if provided.")

    mode = str(cfg["simulation_mode"]).lower()
    if mode != "grid_parameter":
        raise ValueError(
            "simulation.mode='grid_parameter' is now the only supported workflow. "
            "Use singleton value lists for parameters you want to keep fixed."
        )

    if not isinstance(cfg.get("modules"), dict):
        raise ValueError("grid_parameter mode requires simulation.modules.")
    if not bool(cfg["summaries"]):
        raise ValueError("grid_parameter mode currently requires summaries = true.")

    _, axes = _grid_axes_from_modules(cfg["modules"])
    if any(str(v).lower() not in {"fixed", "draw", "map"} for v in axes["bias_mode"]):
        raise ValueError("bias_mode values must be one of: fixed, draw, map")

    total_sets = _grid_total_n_sets(cfg["modules"])
    _normalize_grid_subset(cfg.get("grid_subset"), total_sets=total_sets)
    _normalize_replicate_subset(cfg.get("replicate_subset"), total_replicates=int(cfg["n_sim"]))

    for _, params in _iter_grid_parameter_sets_with_indices(cfg["modules"], grid_subset=cfg.get("grid_subset")):
        _params_with_effective_kappa(cfg, params)

def _merged_topology_plot_cfg(cfg: Dict[str, Any]) -> Dict[str, Any]:
    out = dict(DEFAULT_TOPOLOGY_PLOT)
    out.update(cfg.get("topology_plot", {}))
    return out


def _same_json_file(path: Path, cfg: Dict[str, Any]) -> bool:
    if not path.exists():
        return False
    try:
        old = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return False
    return old == cfg


def _ensure_run_dir_and_config(*, run_dir: Path, raw_cfg: Dict[str, Any]) -> Path:
    config_used_json = run_dir / "config_used.json"

    if not run_dir.exists():
        run_dir.mkdir(parents=True, exist_ok=False)
        config_used_json.write_text(json.dumps(raw_cfg, indent=2, sort_keys=True), encoding="utf-8")
        return config_used_json

    if config_used_json.exists() and not _same_json_file(config_used_json, raw_cfg):
        raise FileExistsError(
            f"Run directory already exists with a different config: {run_dir}\n"
            "Use a new experiment_name or remove/rename the old output folder first."
        )

    if not config_used_json.exists():
        config_used_json.write_text(json.dumps(raw_cfg, indent=2, sort_keys=True), encoding="utf-8")

    return config_used_json


def _checked_topology_json_path(run_dir: Path) -> Path:
    return run_dir / "topology_check" / "topology_internal_steps.json"


def _display_param_value(value: Any) -> Any:
    if value is None:
        return "<none>"
    return value


def _format_label_value(value: Any) -> str:
    if value is None:
        return "none"
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if math.isfinite(value):
            if value.is_integer():
                return str(int(value))
            return format(value, ".6g")
        return str(value)
    return str(value).replace(" ", "")


def _param_id(params: Dict[str, Any]) -> str:
    parts = [
        f"m{_format_label_value(params['m'])}",
        f"rho{_format_label_value(params['rho'])}",
        f"mu{_format_label_value(params['mu_year'])}",
        f"mr{_format_label_value(params['kappa_sr'])}",
        f"vl{_format_label_value(params['victim_locality'])}",
        f"sb{_format_label_value(params['sam_boundary_cells'])}",
        f"bp{_format_label_value(params['branch_precursor_number'])}",
        f"op{_format_label_value(params['organ_precursor_number'])}",
        f"oc{_format_label_value(params['organ_total_cells'])}",
        f"sf{_format_label_value(params['seq_fraction'])}",
        f"mode{_format_label_value(params['bias_mode'])}",
        f"wb{_format_label_value(params['branch_bias_value'])}",
        f"wm{_format_label_value(params['branch_bias_mean'])}",
        f"wk{_format_label_value(params['branch_bias_kappa'])}",
    ]
    return "_".join(parts)




def _build_run_command(
    *,
    cfg: Dict[str, Any],
    pipeline_dir: Path,
    result_json: Path,
    topology_json: Path,
    params: Dict[str, Any],
    n_sim: int,
    seed: int,
) -> List[str]:
    cmd = [
        sys.executable,
        str(pipeline_dir / "pipeline_wrapper.py"),
        "--topology_json", str(topology_json),
        "--out_json", str(result_json),
        "--topology_unit", "steps",
        "--topology_mapping_mode", "deterministic",
        "--m", str(params["m"]),
        "--rho", str(params["rho"]),
        "--mu_year", str(params["mu_year"]),
        "--kappa_sr", str(params["kappa_sr"]),
        "--victim_locality", str(params["victim_locality"]),
        "--sam_boundary_cells", str(params["sam_boundary_cells"]),
        "--branch_precursor_number", str(params["branch_precursor_number_realized_cells"]),
        "--organ_precursor_number", str(params["organ_precursor_number_realized_cells"]),
        "--organ_total_cells", str(params["organ_total_cells_realized"]),
        "--bias_mode", str(params["bias_mode"]),
        "--branch_bias_value", str(params["branch_bias_value"]),
        "--branch_bias_mean", str(params["branch_bias_mean"]),
        "--branch_bias_kappa", str(params["branch_bias_kappa"]),
        "--summaries", "yes" if bool(cfg["summaries"]) else "no",
        "--vaf_nbins", str(cfg["vaf_nbins"]),
        "--summary_private_shared", "yes" if bool(cfg["summary_private_shared"]) else "no",
        "--n_sim", str(n_sim),
        "--seed", str(seed),
        "--sequenced_cells", str(params["sequenced_cells_used"]),
    ]

    phyllotaxy_cfg = cfg.get("phyllotaxy") or {}
    if phyllotaxy_cfg:
        cmd.extend(["--phyllotaxy_mode", str(phyllotaxy_cfg.get("mode", "off"))])
        if phyllotaxy_cfg.get("divergence_angle_deg") is not None:
            cmd.extend(["--phyllotaxy_divergence_deg", str(phyllotaxy_cfg["divergence_angle_deg"])])
        if phyllotaxy_cfg.get("tie_tol") is not None:
            cmd.extend(["--phyllotaxy_tie_tol", str(phyllotaxy_cfg["tie_tol"])])

    return cmd

def _build_check_commands(
    *,
    cfg: Dict[str, Any],
    pipeline_dir: Path,
    run_dir: Path,
    config_path: Path,
) -> Tuple[List[str], Path, Path, Path, Path, Dict[str, Any], Path]:
    topology_json = _resolve_topology_json_path(str(cfg["topology_json"]), base_dir=config_path.parent)
    plot_cfg = _merged_topology_plot_cfg(cfg)

    check_dir = run_dir / "topology_check"
    report_json = check_dir / "topology_report.json"
    internal_topology_json = check_dir / "topology_internal_steps.json"

    # Keep the historical PNG filename so existing downstream scripts and
    # habits do not break when the plotting backend changes.
    topology_png = check_dir / "topology_phylogram_original_units.png"

    inspect_cmd = [
        sys.executable,
        str(pipeline_dir / "inspect_topology.py"),
        "--topology_json", str(topology_json),
        "--topology_unit", str(cfg["mapping_unit"]),
        "--topology_mapping_mode", str(cfg["mapping_mode"]),
        "--export_topology_json", str(internal_topology_json),
        "--report_json", str(report_json),
    ]
    if str(cfg["mapping_unit"]).lower() != "steps":
        inspect_cmd.extend(["--kappa_sr", str(cfg["mapping_rate"])])

    return inspect_cmd, check_dir, report_json, internal_topology_json, topology_png, plot_cfg, topology_json


def _run_topology_plot(
    *,
    plot_cfg: Dict[str, Any],
    topology_json: Path,
    topology_png: Path,
    check_dir: Path,
    pipeline_dir: Path,
) -> Dict[str, Optional[Path]]:
    """Run the cleaned developmental topology plotter for step 'check'."""
    if not bool(plot_cfg.get("enabled", True)):
        return {"png": None, "pdf": None, "layout_csv": None}

    plot_style = str(plot_cfg.get("style", "developmental_pruned")).lower()
    if plot_style != "developmental_pruned":
        raise ValueError("Only check.topology_plot.style='developmental_pruned' is available in this cleaned release.")

    project_root = pipeline_dir.parent
    plotter = project_root / "simSOMA_corefunc" / "plot_topology_json.py"
    if not plotter.exists():
        raise FileNotFoundError(f"Topology plotter not found: {plotter}")

    out_name = topology_png.stem
    topology_pdf = check_dir / f"{out_name}.pdf"
    layout_csv = check_dir / f"{out_name}_layout.csv"
    report_json = check_dir / f"{out_name}_report.json"

    cmd = [
        sys.executable,
        str(plotter),
        "--topology_json", str(topology_json),
        "--outdir", str(check_dir),
        "--out_name", out_name,
        "--fig_width", str(plot_cfg.get("fig_width", DEFAULT_TOPOLOGY_PLOT["fig_width"])),
        "--fig_height", str(plot_cfg.get("fig_height", DEFAULT_TOPOLOGY_PLOT["fig_height"])),
        "--dpi", str(plot_cfg.get("dpi", DEFAULT_TOPOLOGY_PLOT["dpi"])),
        "--y_ticks", str(plot_cfg.get("y_ticks", DEFAULT_TOPOLOGY_PLOT["y_ticks"])),
        "--y_label", str(plot_cfg.get("y_label", DEFAULT_TOPOLOGY_PLOT["y_label"])),
        "--tick_size", str(plot_cfg.get("tick_size", DEFAULT_TOPOLOGY_PLOT["tick_size"])),
        "--label_size", str(plot_cfg.get("label_size", DEFAULT_TOPOLOGY_PLOT["label_size"])),
        "--branch_linewidth", str(plot_cfg.get("branch_linewidth", DEFAULT_TOPOLOGY_PLOT["branch_linewidth"])),
        "--connector_linewidth", str(plot_cfg.get("connector_linewidth", DEFAULT_TOPOLOGY_PLOT["connector_linewidth"])),
        "--branchpoint_size", str(plot_cfg.get("branchpoint_size", DEFAULT_TOPOLOGY_PLOT["branchpoint_size"])),
        "--organ_size", str(plot_cfg.get("organ_size", DEFAULT_TOPOLOGY_PLOT["organ_size"])),
        "--tip_size", str(plot_cfg.get("tip_size", DEFAULT_TOPOLOGY_PLOT["tip_size"])),
        "--min_branch_offset", str(plot_cfg.get("min_branch_offset", DEFAULT_TOPOLOGY_PLOT["min_branch_offset"])),
        "--organ_stub_length", str(plot_cfg.get("organ_stub_length", DEFAULT_TOPOLOGY_PLOT["organ_stub_length"])),
        "--nonterminal_organ_style", str(plot_cfg.get("nonterminal_organ_style", DEFAULT_TOPOLOGY_PLOT["nonterminal_organ_style"])),
        "--nonterminal_organ_size", str(plot_cfg.get("nonterminal_organ_size", DEFAULT_TOPOLOGY_PLOT["nonterminal_organ_size"])),
        "--branch_color_mode", str(plot_cfg.get("branch_color_mode", DEFAULT_TOPOLOGY_PLOT["branch_color_mode"])),
        "--title", str(plot_cfg.get("title", DEFAULT_TOPOLOGY_PLOT["title"])),
    ]
    if bool(plot_cfg.get("show_nonterminal_organs", DEFAULT_TOPOLOGY_PLOT["show_nonterminal_organs"])):
        cmd.append("--show_nonterminal_organs")
    if bool(plot_cfg.get("show_tip_labels", False)):
        cmd.append("--show_tip_labels")

    subprocess.run(cmd, check=True, cwd=str(project_root))

    if not bool(plot_cfg.get("write_pdf", DEFAULT_TOPOLOGY_PLOT["write_pdf"])) and topology_pdf.exists():
        topology_pdf.unlink()
        topology_pdf = None
    if not bool(plot_cfg.get("write_layout_csv", DEFAULT_TOPOLOGY_PLOT["write_layout_csv"])) and layout_csv.exists():
        layout_csv.unlink()
        layout_csv = None

    return {
        "png": topology_png,
        "pdf": topology_pdf if topology_pdf is not None and topology_pdf.exists() else None,
        "layout_csv": layout_csv if layout_csv is not None and layout_csv.exists() else None,
        "report_json": report_json if report_json.exists() else None,
    }


def _run_check_step(*, cfg: Dict[str, Any], raw_cfg: Dict[str, Any], config_path: Path, pipeline_dir: Path) -> None:
    outdir_root = _resolve_outdir_root(str(cfg["outdir_root"]), base_dir=config_path.parent)
    run_dir = outdir_root / str(cfg["experiment_name"])
    config_used_json = _ensure_run_dir_and_config(run_dir=run_dir, raw_cfg=raw_cfg)
    software_version_json = _write_software_version(run_dir, pipeline_dir)

    inspect_cmd, check_dir, report_json, internal_topology_json, topology_png, plot_cfg, topology_json = _build_check_commands(
        cfg=cfg,
        pipeline_dir=pipeline_dir,
        run_dir=run_dir,
        config_path=config_path,
    )

    if check_dir.exists():
        shutil.rmtree(check_dir)

    check_dir.mkdir(parents=True, exist_ok=True)

    try:
        subprocess.run(inspect_cmd, check=True, cwd=str(pipeline_dir))
        topology_plot_outputs = _run_topology_plot(
            plot_cfg=plot_cfg,
            topology_json=topology_json,
            topology_png=topology_png,
            check_dir=check_dir,
            pipeline_dir=pipeline_dir,
        )
    except Exception as exc:
        shutil.rmtree(check_dir, ignore_errors=True)
        if run_dir.exists() and not any(run_dir.iterdir()):
            shutil.rmtree(run_dir, ignore_errors=True)
        raise RuntimeError("Topology check failed. Partial topology_check outputs were removed.") from exc

    print("Completed step: check")
    print(f"run_dir: {run_dir}")
    print(f"config_used_json: {config_used_json}")
    print(f"software_version_json: {software_version_json}")
    print(f"topology_report_json: {report_json}")
    print(f"topology_internal_steps_json: {internal_topology_json}")
    if topology_plot_outputs.get("png") is not None:
        print(f"topology_plot_png: {topology_plot_outputs['png']}")
    else:
        print("topology_plot: disabled")
    if topology_plot_outputs.get("pdf") is not None:
        print(f"topology_plot_pdf: {topology_plot_outputs['pdf']}")
    if topology_plot_outputs.get("layout_csv") is not None:
        print(f"topology_plot_layout_csv: {topology_plot_outputs['layout_csv']}")

def _summarize_replicate_row(
    *,
    set_id: str,
    set_label: str,
    rep: int,
    seed: int | None,
    cfg: Dict[str, Any],
    params: Dict[str, Any],
    summ: Dict[str, Any],
    organ_meta: Dict[str, Dict[str, Any]] | None = None,
) -> Dict[str, Any]:
    agg = summ.get("aggregate", {})
    row: Dict[str, Any] = {
        "set_id": set_id,
        "set_label": set_label,
        "rep": int(rep),
        "seed": seed,
        "total_variants": int(agg.get("total_variants", 0)),
        "total_unique_variants": int(agg.get("total_unique_variants", agg.get("total_variants", 0))),
        "total_variant_occurrences": int(agg.get("total_variant_occurrences", 0)),
        **_topology_row(cfg),
    }
    for k in OUTPUT_PARAM_FIELDS:
        row[k] = _display_param_value(params.get(k))
    bins = agg.get("vaf_bins", [])
    for i, v in enumerate(bins):
        row[f"bin_{i}"] = int(v)
    ps = summ.get("private_shared")
    if isinstance(ps, dict) and ps.get("available"):
        row["n_shared_mutations"] = int(ps.get("n_shared_mutations", 0))
        row["n_private_mutations_total"] = int(ps.get("n_private_mutations_total", 0))
    return row


def _summarize_organ_replicate_rows(
    *,
    set_id: str,
    set_label: str,
    rep: int,
    seed: int | None,
    cfg: Dict[str, Any],
    params: Dict[str, Any],
    summ: Dict[str, Any],
    organ_meta: Dict[str, Dict[str, Any]],
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    per_organ = summ.get("per_organ", {})
    for oid in sorted(per_organ.keys(), key=_natural_key):
        organ_summary = per_organ[oid]
        row: Dict[str, Any] = {
            "set_id": set_id,
            "set_label": set_label,
            "rep": int(rep),
            "seed": seed,
            "organ_id": oid,
            "n_bottlenecks_root_to_organ": organ_meta.get(oid, {}).get("n_bottlenecks_root_to_organ"),
            "sequenced_cells": int(organ_summary.get("sequenced_cells", 0)),
            "total_variants": int(organ_summary.get("total_variants", 0)),
            "singletons": int(organ_summary.get("singletons", 0)),
            "fixed": int(organ_summary.get("fixed", 0)),
            "n_lineages_sampled": int(organ_summary.get("n_lineages_sampled", 0)),
            "dominant_lineage_id": int(organ_summary.get("dominant_lineage_id", -1)),
            "dominant_fraction": float(organ_summary.get("dominant_fraction", 0.0)),
            "lineage_entropy": float(organ_summary.get("lineage_entropy", 0.0)),
            "boundary_cells": int(organ_summary.get("boundary_cells", 0)),
            "organ_total_cells": int(organ_summary.get("organ_total_cells", 0)),
            "organ_precursor_number": int(organ_summary.get("organ_precursor_number", 0)),
            **_topology_row(cfg),
        }
        for k in OUTPUT_PARAM_FIELDS:
            row[k] = _display_param_value(params.get(k))
        bins = organ_summary.get("vaf_bins", [])
        for i, v in enumerate(bins):
            row[f"bin_{i}"] = int(v)
        rows.append(row)
    return rows


def _aggregate_organ_summary_rows(
    *,
    set_id: str,
    set_label: str,
    cfg: Dict[str, Any],
    params: Dict[str, Any],
    group_summary: Dict[str, Any],
    n_replicates: int,
    organ_meta: Dict[str, Dict[str, Any]],
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    per_organ = group_summary.get("per_organ", {})
    for oid in sorted(per_organ.keys(), key=_natural_key):
        organ_summary = per_organ[oid]
        row: Dict[str, Any] = {
            "set_id": set_id,
            "set_label": set_label,
            "n_replicates": int(n_replicates),
            "organ_id": oid,
            "n_bottlenecks_root_to_organ": organ_meta.get(oid, {}).get("n_bottlenecks_root_to_organ"),
            "sequenced_cells": int(organ_summary.get("sequenced_cells", 0)),
            "total_variants_mean": organ_summary.get("total_variants", {}).get("mean"),
            "total_variants_sd": organ_summary.get("total_variants", {}).get("sd"),
            "singletons_mean": organ_summary.get("singletons", {}).get("mean"),
            "singletons_sd": organ_summary.get("singletons", {}).get("sd"),
            "fixed_mean": organ_summary.get("fixed", {}).get("mean"),
            "fixed_sd": organ_summary.get("fixed", {}).get("sd"),
            "n_lineages_sampled_mean": organ_summary.get("n_lineages_sampled", {}).get("mean"),
            "n_lineages_sampled_sd": organ_summary.get("n_lineages_sampled", {}).get("sd"),
            "dominant_lineage_id_mean": organ_summary.get("dominant_lineage_id", {}).get("mean"),
            "dominant_lineage_id_sd": organ_summary.get("dominant_lineage_id", {}).get("sd"),
            "dominant_fraction_mean": organ_summary.get("dominant_fraction", {}).get("mean"),
            "dominant_fraction_sd": organ_summary.get("dominant_fraction", {}).get("sd"),
            "lineage_entropy_mean": organ_summary.get("lineage_entropy", {}).get("mean"),
            "lineage_entropy_sd": organ_summary.get("lineage_entropy", {}).get("sd"),
            "boundary_cells_mean": organ_summary.get("boundary_cells", {}).get("mean"),
            "boundary_cells_sd": organ_summary.get("boundary_cells", {}).get("sd"),
            "organ_total_cells_mean": organ_summary.get("organ_total_cells", {}).get("mean"),
            "organ_total_cells_sd": organ_summary.get("organ_total_cells", {}).get("sd"),
            "organ_precursor_number_mean": organ_summary.get("organ_precursor_number", {}).get("mean"),
            "organ_precursor_number_sd": organ_summary.get("organ_precursor_number", {}).get("sd"),
            **_topology_row(cfg),
        }
        for k in OUTPUT_PARAM_FIELDS:
            row[k] = _display_param_value(params.get(k))
        vaf_bins = organ_summary.get("vaf_bins", {})
        mean_bins = vaf_bins.get("mean", []) if isinstance(vaf_bins, dict) else []
        sd_bins = vaf_bins.get("sd", []) if isinstance(vaf_bins, dict) else []
        for i, v in enumerate(mean_bins):
            row[f"bin_{i}_mean"] = v
        for i, v in enumerate(sd_bins):
            row[f"bin_{i}_sd"] = v
        rows.append(row)
    return rows


def _raw_vaf_rows_for_sim(
    *,
    set_id: str,
    set_label: str,
    rep: int,
    seed: int | None,
    result: Dict[str, Any],
    organ_meta: Dict[str, Dict[str, Any]],
    params: Dict[str, Any],
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    organ_events = list(result.get("organ_events", []))
    organ_events.sort(key=lambda oe: _natural_key(str(oe.get("organ_id", ""))))
    for oe in organ_events:
        organ_id = str(oe.get("organ_id"))
        sequenced_cells = int(oe.get("sequenced_cells", 0))
        acm = oe.get("allele_counts_by_mutation", {})
        if not isinstance(acm, dict) or sequenced_cells <= 0:
            continue
        for mid, cnt in sorted(((int(k), int(v)) for k, v in acm.items()), key=lambda x: x[0]):
            if cnt <= 0:
                continue
            rows.append({
                "set_id": set_id,
                "set_label": set_label,
                "rep": int(rep),
                "seed": seed,
                "organ_id": organ_id,
                "n_bottlenecks_root_to_organ": organ_meta.get(organ_id, {}).get("n_bottlenecks_root_to_organ"),
                "seq_fraction": _display_param_value(params.get("seq_fraction")),
                "organ_total_cells": params.get("organ_total_cells"),
                "sequenced_cells_used": params.get("sequenced_cells_used"),
                "mutation_id": mid,
                "allele_count": int(cnt),
                "sequenced_cells": sequenced_cells,
                "vaf": float(cnt) / float(sequenced_cells),
            })
    return rows



def _bin_edge_row(bin_index: int, nbins: int) -> Dict[str, Any]:
    return {
        "bin_index": int(bin_index),
        "vaf_left": float(bin_index) / float(nbins),
        "vaf_right": float(bin_index + 1) / float(nbins),
    }


def _feature_base_row(*, set_id: str, set_label: str, cfg: Dict[str, Any], params: Dict[str, Any]) -> Dict[str, Any]:
    row: Dict[str, Any] = {
        "set_id": set_id,
        "set_label": set_label,
        **_topology_row(cfg),
    }
    for k in OUTPUT_PARAM_FIELDS:
        row[k] = _display_param_value(params.get(k))
    return row


def _append_vaf_bins_replicate(
    rows: List[Dict[str, Any]],
    *,
    base: Dict[str, Any],
    rep: int,
    seed: int | None,
    level: str,
    organ_id: str | None,
    organ_meta: Dict[str, Dict[str, Any]],
    variant_class: str,
    bins: List[Any],
    nbins: int,
) -> None:
    for i, value in enumerate(bins):
        row = dict(base)
        row.update({
            "rep": int(rep),
            "seed": seed,
            "level": level,
            "organ_id": organ_id or "",
            "n_bottlenecks_root_to_organ": organ_meta.get(str(organ_id), {}).get("n_bottlenecks_root_to_organ") if organ_id else "",
            "variant_class": variant_class,
            **_bin_edge_row(i, nbins),
            "n_variants": int(value),
        })
        rows.append(row)


def _append_vaf_bins_aggregate(
    rows: List[Dict[str, Any]],
    *,
    base: Dict[str, Any],
    n_replicates: int,
    level: str,
    organ_id: str | None,
    organ_meta: Dict[str, Dict[str, Any]],
    variant_class: str,
    bins_stats: Dict[str, Any],
    nbins: int,
) -> None:
    means = list(bins_stats.get("mean", [])) if isinstance(bins_stats, dict) else []
    sds = list(bins_stats.get("sd", [])) if isinstance(bins_stats, dict) else []
    for i, mean_value in enumerate(means):
        row = dict(base)
        row.update({
            "n_replicates": int(n_replicates),
            "level": level,
            "organ_id": organ_id or "",
            "n_bottlenecks_root_to_organ": organ_meta.get(str(organ_id), {}).get("n_bottlenecks_root_to_organ") if organ_id else "",
            "variant_class": variant_class,
            **_bin_edge_row(i, nbins),
            "n_variants_mean": mean_value,
            "n_variants_sd": sds[i] if i < len(sds) else None,
        })
        rows.append(row)


def _vaf_class_replicate_rows(
    *,
    set_id: str,
    set_label: str,
    rep: int,
    seed: int | None,
    cfg: Dict[str, Any],
    params: Dict[str, Any],
    summ: Dict[str, Any],
    organ_meta: Dict[str, Dict[str, Any]],
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    nbins = int(summ.get("nbins", cfg.get("vaf_nbins", 20)))
    base = _feature_base_row(set_id=set_id, set_label=set_label, cfg=cfg, params=params)
    agg = summ.get("aggregate", {})
    ps = summ.get("private_shared", {}) if isinstance(summ.get("private_shared"), dict) else {}

    _append_vaf_bins_replicate(rows, base=base, rep=rep, seed=seed, level="global", organ_id=None, organ_meta=organ_meta, variant_class="all", bins=list(agg.get("vaf_bins_organ_occurrences", agg.get("vaf_bins", []))), nbins=nbins)
    if ps.get("available"):
        _append_vaf_bins_replicate(rows, base=base, rep=rep, seed=seed, level="global", organ_id=None, organ_meta=organ_meta, variant_class="private", bins=list(ps.get("private_vaf_bins_global", [])), nbins=nbins)
        _append_vaf_bins_replicate(rows, base=base, rep=rep, seed=seed, level="global", organ_id=None, organ_meta=organ_meta, variant_class="shared", bins=list(ps.get("shared_vaf_bins_global", [])), nbins=nbins)

    per_organ = summ.get("per_organ", {})
    for oid in sorted(per_organ.keys(), key=_natural_key):
        os = per_organ[oid]
        _append_vaf_bins_replicate(rows, base=base, rep=rep, seed=seed, level="organ", organ_id=oid, organ_meta=organ_meta, variant_class="all", bins=list(os.get("vaf_bins", [])), nbins=nbins)
        if ps.get("available"):
            _append_vaf_bins_replicate(rows, base=base, rep=rep, seed=seed, level="organ", organ_id=oid, organ_meta=organ_meta, variant_class="private", bins=list(ps.get("private_vaf_bins_by_organ", {}).get(oid, [])), nbins=nbins)
            _append_vaf_bins_replicate(rows, base=base, rep=rep, seed=seed, level="organ", organ_id=oid, organ_meta=organ_meta, variant_class="shared", bins=list(ps.get("shared_vaf_bins_by_organ", {}).get(oid, [])), nbins=nbins)
    return rows


def _vaf_class_aggregate_rows(
    *,
    set_id: str,
    set_label: str,
    cfg: Dict[str, Any],
    params: Dict[str, Any],
    group_summary: Dict[str, Any],
    n_replicates: int,
    organ_meta: Dict[str, Dict[str, Any]],
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    nbins = int(group_summary.get("nbins", cfg.get("vaf_nbins", 20)))
    base = _feature_base_row(set_id=set_id, set_label=set_label, cfg=cfg, params=params)
    agg = group_summary.get("aggregate", {})
    ps = group_summary.get("private_shared", {}) if isinstance(group_summary.get("private_shared"), dict) else {}

    _append_vaf_bins_aggregate(rows, base=base, n_replicates=n_replicates, level="global", organ_id=None, organ_meta=organ_meta, variant_class="all", bins_stats=agg.get("vaf_bins_organ_occurrences", {}), nbins=nbins)
    if ps.get("available"):
        _append_vaf_bins_aggregate(rows, base=base, n_replicates=n_replicates, level="global", organ_id=None, organ_meta=organ_meta, variant_class="private", bins_stats=ps.get("private_vaf_bins_global", {}), nbins=nbins)
        _append_vaf_bins_aggregate(rows, base=base, n_replicates=n_replicates, level="global", organ_id=None, organ_meta=organ_meta, variant_class="shared", bins_stats=ps.get("shared_vaf_bins_global", {}), nbins=nbins)

    per_organ = group_summary.get("per_organ", {})
    for oid in sorted(per_organ.keys(), key=_natural_key):
        os = per_organ[oid]
        _append_vaf_bins_aggregate(rows, base=base, n_replicates=n_replicates, level="organ", organ_id=oid, organ_meta=organ_meta, variant_class="all", bins_stats=os.get("vaf_bins", {}), nbins=nbins)
        if ps.get("available"):
            _append_vaf_bins_aggregate(rows, base=base, n_replicates=n_replicates, level="organ", organ_id=oid, organ_meta=organ_meta, variant_class="private", bins_stats=ps.get("private_vaf_bins_by_organ", {}).get(oid, {}), nbins=nbins)
            _append_vaf_bins_aggregate(rows, base=base, n_replicates=n_replicates, level="organ", organ_id=oid, organ_meta=organ_meta, variant_class="shared", bins_stats=ps.get("shared_vaf_bins_by_organ", {}).get(oid, {}), nbins=nbins)
    return rows


def _vaf_count_spectrum_replicate_rows(
    *,
    set_id: str,
    set_label: str,
    rep: int,
    seed: int | None,
    cfg: Dict[str, Any],
    params: Dict[str, Any],
    summ: Dict[str, Any],
    organ_meta: Dict[str, Dict[str, Any]],
) -> List[Dict[str, Any]]:
    spec = summ.get("vaf_count_spectrum", {}) if isinstance(summ.get("vaf_count_spectrum"), dict) else {}
    if not spec.get("available"):
        return []
    base = _feature_base_row(set_id=set_id, set_label=set_label, cfg=cfg, params=params)
    rows: List[Dict[str, Any]] = []
    for r in spec.get("rows", []):
        row = dict(base)
        organ_id = str(r.get("organ_id", ""))
        row.update({
            "rep": int(rep),
            "seed": seed,
            "level": r.get("level", ""),
            "organ_id": organ_id,
            "n_bottlenecks_root_to_organ": organ_meta.get(organ_id, {}).get("n_bottlenecks_root_to_organ") if organ_id else "",
            "variant_class": r.get("variant_class", ""),
            "sharing_degree": r.get("sharing_degree", ""),
            "n_sampled_cells": r.get("n_sampled_cells", ""),
            "allele_count": r.get("allele_count", ""),
            "sampled_vaf": r.get("sampled_vaf", ""),
            "n_variants": int(r.get("n_variants", 0)),
        })
        rows.append(row)
    return rows


def _vaf_count_spectrum_aggregate_rows(
    *,
    set_id: str,
    set_label: str,
    cfg: Dict[str, Any],
    params: Dict[str, Any],
    group_summary: Dict[str, Any],
    n_replicates: int,
    organ_meta: Dict[str, Dict[str, Any]],
) -> List[Dict[str, Any]]:
    spec = group_summary.get("vaf_count_spectrum", {}) if isinstance(group_summary.get("vaf_count_spectrum"), dict) else {}
    if not spec.get("available"):
        return []
    base = _feature_base_row(set_id=set_id, set_label=set_label, cfg=cfg, params=params)
    rows: List[Dict[str, Any]] = []
    for r in spec.get("rows", []):
        row = dict(base)
        organ_id = str(r.get("organ_id", ""))
        row.update({
            "n_replicates": int(n_replicates),
            "level": r.get("level", ""),
            "organ_id": organ_id,
            "n_bottlenecks_root_to_organ": organ_meta.get(organ_id, {}).get("n_bottlenecks_root_to_organ") if organ_id else "",
            "variant_class": r.get("variant_class", ""),
            "sharing_degree": r.get("sharing_degree", ""),
            "n_sampled_cells": r.get("n_sampled_cells", ""),
            "allele_count": r.get("allele_count", ""),
            "sampled_vaf": r.get("sampled_vaf", ""),
            "n_variants_mean": r.get("n_variants_mean", ""),
            "n_variants_sd": r.get("n_variants_sd", ""),
        })
        rows.append(row)
    return rows


def _sharing_replicate_rows(
    *,
    set_id: str,
    set_label: str,
    rep: int,
    seed: int | None,
    cfg: Dict[str, Any],
    params: Dict[str, Any],
    summ: Dict[str, Any],
) -> List[Dict[str, Any]]:
    ps = summ.get("private_shared", {}) if isinstance(summ.get("private_shared"), dict) else {}
    if not ps.get("available"):
        return []
    rows: List[Dict[str, Any]] = []
    base = _feature_base_row(set_id=set_id, set_label=set_label, cfg=cfg, params=params)
    def add(**kw: Any) -> None:
        row = dict(base)
        row.update({"rep": int(rep), "seed": seed})
        row.update(kw)
        rows.append(row)
    add(statistic="global_count", variant_class="all_unique", n_variants=int(ps.get("n_mutations_total", 0)))
    add(statistic="global_count", variant_class="private", n_variants=int(ps.get("n_private_mutations_total", 0)))
    add(statistic="global_count", variant_class="shared", n_variants=int(ps.get("n_shared_mutations", 0)))
    for degree, value in sorted(ps.get("sharing_degree_counts", {}).items(), key=lambda kv: int(kv[0])):
        add(statistic="sharing_degree", sharing_degree=int(degree), variant_class="degree", n_variants=int(value))
    for oid, value in sorted(ps.get("n_private_by_organ", {}).items(), key=lambda kv: _natural_key(kv[0])):
        add(statistic="private_by_organ", variant_class="private", organ_id=oid, n_variants=int(value))
    for key, value in sorted(ps.get("pairwise_shared_counts", {}).items(), key=lambda kv: _natural_key(kv[0])):
        a, b = str(key).split("|", 1) if "|" in str(key) else (str(key), "")
        add(statistic="pairwise_shared", variant_class="shared", organ_a=a, organ_b=b, n_variants=int(value))
    return rows


def _sharing_aggregate_rows(
    *,
    set_id: str,
    set_label: str,
    cfg: Dict[str, Any],
    params: Dict[str, Any],
    group_summary: Dict[str, Any],
    n_replicates: int,
) -> List[Dict[str, Any]]:
    ps = group_summary.get("private_shared", {}) if isinstance(group_summary.get("private_shared"), dict) else {}
    if not ps.get("available"):
        return []
    rows: List[Dict[str, Any]] = []
    base = _feature_base_row(set_id=set_id, set_label=set_label, cfg=cfg, params=params)
    def stats_value(d: Dict[str, Any], field: str) -> Any:
        return d.get(field) if isinstance(d, dict) else None
    def add(**kw: Any) -> None:
        row = dict(base)
        row.update({"n_replicates": int(n_replicates)})
        row.update(kw)
        rows.append(row)
    agg = group_summary.get("aggregate", {}) if isinstance(group_summary.get("aggregate"), dict) else {}
    add(statistic="global_count", variant_class="all_unique", n_variants_mean=stats_value(agg.get("total_unique_variants", {}), "mean"), n_variants_sd=stats_value(agg.get("total_unique_variants", {}), "sd"))
    add(statistic="global_count", variant_class="private", n_variants_mean=stats_value(ps.get("n_private_mutations_total", {}), "mean"), n_variants_sd=stats_value(ps.get("n_private_mutations_total", {}), "sd"))
    add(statistic="global_count", variant_class="shared", n_variants_mean=stats_value(ps.get("n_shared_mutations", {}), "mean"), n_variants_sd=stats_value(ps.get("n_shared_mutations", {}), "sd"))
    for degree, st in sorted(ps.get("sharing_degree_counts", {}).items(), key=lambda kv: int(kv[0])):
        add(statistic="sharing_degree", sharing_degree=int(degree), variant_class="degree", n_variants_mean=stats_value(st, "mean"), n_variants_sd=stats_value(st, "sd"))
    for oid, st in sorted(ps.get("n_private_by_organ", {}).items(), key=lambda kv: _natural_key(kv[0])):
        add(statistic="private_by_organ", variant_class="private", organ_id=oid, n_variants_mean=stats_value(st, "mean"), n_variants_sd=stats_value(st, "sd"))
    for key, st in sorted(ps.get("pairwise_shared_counts", {}).items(), key=lambda kv: _natural_key(kv[0])):
        a, b = str(key).split("|", 1) if "|" in str(key) else (str(key), "")
        add(statistic="pairwise_shared", variant_class="shared", organ_a=a, organ_b=b, n_variants_mean=stats_value(st, "mean"), n_variants_sd=stats_value(st, "sd"))
    return rows


def _aggregate_summary_row(
    *,
    set_id: str,
    set_label: str,
    cfg: Dict[str, Any],
    params: Dict[str, Any],
    group_summary: Dict[str, Any],
    n_replicates: int,
) -> Dict[str, Any]:
    row: Dict[str, Any] = {
        "set_id": set_id,
        "set_label": set_label,
        "n_replicates": int(n_replicates),
        "private_shared_available": bool(group_summary.get("private_shared", {}).get("available", False)),
        **_topology_row(cfg),
    }
    for k in OUTPUT_PARAM_FIELDS:
        row[k] = _display_param_value(params.get(k))

    agg = group_summary.get("aggregate", {})
    if isinstance(agg, dict):
        row["total_variants_mean"] = agg.get("total_variants", {}).get("mean")
        row["total_variants_sd"] = agg.get("total_variants", {}).get("sd")
        row["total_unique_variants_mean"] = agg.get("total_unique_variants", {}).get("mean")
        row["total_unique_variants_sd"] = agg.get("total_unique_variants", {}).get("sd")
        row["total_variant_occurrences_mean"] = agg.get("total_variant_occurrences", {}).get("mean")
        row["total_variant_occurrences_sd"] = agg.get("total_variant_occurrences", {}).get("sd")

    ps = group_summary.get("private_shared", {})
    if isinstance(ps, dict) and ps.get("available"):
        row["n_shared_mutations_mean"] = ps.get("n_shared_mutations", {}).get("mean")
        row["n_shared_mutations_sd"] = ps.get("n_shared_mutations", {}).get("sd")
        row["n_private_mutations_total_mean"] = ps.get("n_private_mutations_total", {}).get("mean")
        row["n_private_mutations_total_sd"] = ps.get("n_private_mutations_total", {}).get("sd")
        for oid, stats in sorted(ps.get("n_private_by_organ", {}).items(), key=lambda kv: _natural_key(kv[0])):
            row[f"{oid}_private_mean"] = stats.get("mean")
            row[f"{oid}_private_sd"] = stats.get("sd")
    else:
        row["private_shared_reason"] = ps.get("reason") if isinstance(ps, dict) else None
    return row


def _iter_grid_parameter_sets(modules: Dict[str, Any]) -> Iterable[Dict[str, Any]]:
    param_names, axes = _grid_axes_from_modules(modules)
    value_lists = [axes[p] for p in param_names]
    for combo in itertools.product(*value_lists):
        yield {p: v for p, v in zip(param_names, combo)}


def _grid_total_n_sets(modules: Dict[str, Any]) -> int:
    _, axes = _grid_axes_from_modules(modules)
    total = 1
    for vals in axes.values():
        total *= len(vals)
    return total




def _run_grid_step(*, cfg: Dict[str, Any], raw_cfg: Dict[str, Any], config_path: Path, pipeline_dir: Path) -> None:
    outdir_root = _resolve_outdir_root(str(cfg["outdir_root"]), base_dir=config_path.parent)
    run_dir = outdir_root / str(cfg["experiment_name"])
    config_used_json = _ensure_run_dir_and_config(run_dir=run_dir, raw_cfg=raw_cfg)
    software_version_json = _write_software_version(run_dir, pipeline_dir)
    checked_topology_json = _checked_topology_json_path(run_dir)

    if not checked_topology_json.exists():
        raise FileNotFoundError(
            f"Checked topology not found: {checked_topology_json}\n"
            "Run step 'check' first, inspect the topology figure, and then rerun step 'run'."
        )

    grid_dir = run_dir / "grid_parameter"
    if grid_dir.exists():
        raise FileExistsError(
            f"Grid output directory already exists: {grid_dir}\n"
            "Use a new experiment_name or remove/rename that folder before rerunning step 'run'."
        )
    grid_dir.mkdir(parents=True, exist_ok=False)

    checked_topology = _load_topology_json(checked_topology_json)
    organ_meta = _organ_topology_metadata(checked_topology)

    store_full = bool(cfg.get("store_full_results", False))
    export_raw_vafs = bool(cfg.get("export_raw_vafs", False))
    export_sharing_summaries = bool(cfg.get("export_sharing_summaries", True))
    export_replicate_sharing_summaries = bool(cfg.get("export_replicate_sharing_summaries", False))
    export_vaf_count_spectra = bool(cfg.get("export_vaf_count_spectra", True))
    export_legacy_binned_vaf_summaries = bool(cfg.get("export_legacy_binned_vaf_summaries", False))
    export_legacy_binned_vaf_replicates = bool(cfg.get("export_legacy_binned_vaf_replicates", False))
    full_results_dir = grid_dir / "full_results"
    tmp_results_dir = grid_dir / "_tmp_results"
    if store_full:
        full_results_dir.mkdir(parents=True, exist_ok=True)
    else:
        tmp_results_dir.mkdir(parents=True, exist_ok=True)

    def make_writer(filename: str, kind: str) -> _StreamingCsvWriter:
        return _StreamingCsvWriter(
            grid_dir / filename,
            kind=kind,
            fieldnames=_streaming_fieldnames(kind, cfg=cfg, organ_meta=organ_meta),
        )

    writers: Dict[str, _StreamingCsvWriter] = {
        "parameter_sets": make_writer("parameter_sets.csv", "parameter_sets"),
        "replicate_summaries": make_writer("replicate_summaries.csv", "replicate_summaries"),
        "aggregated_summaries": make_writer("aggregated_summaries.csv", "aggregated_summaries"),
        "organ_replicate_summaries": make_writer("organ_replicate_summaries.csv", "organ_replicate_summaries"),
        "organ_aggregated_summaries": make_writer("organ_aggregated_summaries.csv", "organ_aggregated_summaries"),
    }
    if export_legacy_binned_vaf_summaries:
        writers["vaf_class_aggregated_summaries"] = make_writer("vaf_class_aggregated_summaries.csv", "vaf_class_summaries")
    if export_sharing_summaries:
        writers["sharing_aggregated_summaries"] = make_writer("sharing_aggregated_summaries.csv", "sharing_summaries")
    if export_legacy_binned_vaf_replicates:
        writers["vaf_class_replicate_summaries"] = make_writer("vaf_class_replicate_summaries.csv.gz", "vaf_class_summaries")
    if export_replicate_sharing_summaries:
        writers["sharing_replicate_summaries"] = make_writer("sharing_replicate_summaries.csv.gz", "sharing_summaries")
    if export_vaf_count_spectra:
        writers["vaf_count_spectrum_aggregated_summaries"] = make_writer("vaf_count_spectrum_aggregated_summaries.csv", "vaf_count_spectrum_summaries")
        writers["vaf_count_spectrum_replicate_summaries"] = make_writer("vaf_count_spectrum_replicate_summaries.csv.gz", "vaf_count_spectrum_summaries")
    if export_raw_vafs:
        writers["raw_vafs"] = make_writer("raw_vafs.csv.gz", "raw_vafs")

    total_sets_full = _grid_total_n_sets(cfg["modules"])
    subset_info = _normalize_grid_subset(cfg.get("grid_subset"), total_sets=total_sets_full)
    total_sets_selected = _grid_subset_n_sets(cfg["modules"], grid_subset=cfg.get("grid_subset"))
    total_replicates_full = int(cfg["n_sim"])
    replicate_subset_info = _normalize_replicate_subset(cfg.get("replicate_subset"), total_replicates=total_replicates_full)
    if replicate_subset_info is None:
        replicate_start_index = 1
        replicate_end_index = total_replicates_full
        n_sim_selected = total_replicates_full
    else:
        replicate_start_index = int(replicate_subset_info["start_index"])
        replicate_end_index = int(replicate_subset_info["end_index"])
        n_sim_selected = replicate_end_index - replicate_start_index + 1

    try:
        for idx, params in _iter_grid_parameter_sets_with_indices(cfg["modules"], grid_subset=cfg.get("grid_subset")):
            set_id = f"set_{idx:05d}"
            params_eff = _params_with_effective_kappa(cfg, params)
            set_label = _param_id(params_eff)
            parameter_row: Dict[str, Any] = {
                "set_id": set_id,
                "set_label": set_label,
                "n_sim": int(n_sim_selected),
                **_topology_row(cfg),
            }
            for k in OUTPUT_PARAM_FIELDS:
                parameter_row[k] = _display_param_value(params_eff.get(k))
            writers["parameter_sets"].write_row(parameter_row)

            summary_acc = summaries.SummaryAccumulator()

            for rep_index in range(int(replicate_start_index), int(replicate_end_index) + 1):
                rep_abs = int(rep_index) - 1
                result_json = (full_results_dir if store_full else tmp_results_dir) / f"{set_id}_rep_{rep_abs:05d}.json"
                rep_seed = int(cfg["seed"]) + idx * 1000000 + rep_abs
                cmd = _build_run_command(
                    cfg=cfg,
                    pipeline_dir=pipeline_dir,
                    result_json=result_json,
                    topology_json=checked_topology_json,
                    params=params_eff,
                    n_sim=1,
                    seed=rep_seed,
                )

                try:
                    subprocess.run(cmd, check=True, cwd=str(pipeline_dir))
                except subprocess.CalledProcessError as exc:
                    result_json.unlink(missing_ok=True)
                    raise RuntimeError(f"pipeline_wrapper.py failed for {set_id} ({set_label}), rep={rep_abs}.") from exc

                payload = json.loads(result_json.read_text(encoding="utf-8"))
                runs = payload.get("runs", [])
                if len(runs) != 1:
                    raise RuntimeError(f"Unexpected output structure for {set_id}, rep={rep_abs}: expected exactly one run entry.")
                run_entry = runs[0]
                sims = list(run_entry.get("sims", []))
                if len(sims) != 1:
                    raise RuntimeError(f"Unexpected output structure for {set_id}, rep={rep_abs}: expected exactly one simulation.")
                sim = sims[0]
                sim["rep"] = rep_abs
                seed = sim.get("seed")
                if not isinstance(sim.get("summary"), dict):
                    raise RuntimeError(f"Missing replicate summary for {set_id}, rep={rep_abs}. grid_parameter mode requires summaries=true.")

                if store_full:
                    payload["runs"][0]["sims"] = [sim]
                    payload["runs"][0]["summary_aggregate"] = None
                    result_json.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")

                summ = sim["summary"]
                summary_acc.update(summ)

                writers["replicate_summaries"].write_row(
                    _summarize_replicate_row(
                        set_id=set_id,
                        set_label=set_label,
                        rep=rep_abs,
                        seed=seed,
                        cfg=cfg,
                        params=params_eff,
                        summ=summ,
                        organ_meta=organ_meta,
                    )
                )
                writers["organ_replicate_summaries"].write_rows(
                    _summarize_organ_replicate_rows(
                        set_id=set_id,
                        set_label=set_label,
                        rep=rep_abs,
                        seed=seed,
                        cfg=cfg,
                        params=params_eff,
                        summ=summ,
                        organ_meta=organ_meta,
                    )
                )
                if export_legacy_binned_vaf_replicates:
                    writers["vaf_class_replicate_summaries"].write_rows(
                        _vaf_class_replicate_rows(
                            set_id=set_id,
                            set_label=set_label,
                            rep=rep_abs,
                            seed=seed,
                            cfg=cfg,
                            params=params_eff,
                            summ=summ,
                            organ_meta=organ_meta,
                        )
                    )
                if export_replicate_sharing_summaries:
                    writers["sharing_replicate_summaries"].write_rows(
                        _sharing_replicate_rows(
                            set_id=set_id,
                            set_label=set_label,
                            rep=rep_abs,
                            seed=seed,
                            cfg=cfg,
                            params=params_eff,
                            summ=summ,
                        )
                    )
                if export_vaf_count_spectra:
                    writers["vaf_count_spectrum_replicate_summaries"].write_rows(
                        _vaf_count_spectrum_replicate_rows(
                            set_id=set_id,
                            set_label=set_label,
                            rep=rep_abs,
                            seed=seed,
                            cfg=cfg,
                            params=params_eff,
                            summ=summ,
                            organ_meta=organ_meta,
                        )
                    )
                if export_raw_vafs and isinstance(sim.get("result"), dict):
                    writers["raw_vafs"].write_rows(
                        _raw_vaf_rows_for_sim(
                            set_id=set_id,
                            set_label=set_label,
                            rep=rep_abs,
                            seed=seed,
                            result=sim["result"],
                            organ_meta=organ_meta,
                            params=params_eff,
                        )
                    )

                if not store_full:
                    result_json.unlink(missing_ok=True)
                del payload, run_entry, sims, sim, summ
                gc.collect()

            group_summary = summary_acc.finalize()
            n_replicates_done = int(summary_acc.n)
            writers["aggregated_summaries"].write_row(
                _aggregate_summary_row(
                    set_id=set_id,
                    set_label=set_label,
                    cfg=cfg,
                    params=params_eff,
                    group_summary=group_summary,
                    n_replicates=n_replicates_done,
                )
            )
            writers["organ_aggregated_summaries"].write_rows(
                _aggregate_organ_summary_rows(
                    set_id=set_id,
                    set_label=set_label,
                    cfg=cfg,
                    params=params_eff,
                    group_summary=group_summary,
                    n_replicates=n_replicates_done,
                    organ_meta=organ_meta,
                )
            )
            if export_legacy_binned_vaf_summaries:
                writers["vaf_class_aggregated_summaries"].write_rows(
                    _vaf_class_aggregate_rows(
                        set_id=set_id,
                        set_label=set_label,
                        cfg=cfg,
                        params=params_eff,
                        group_summary=group_summary,
                        n_replicates=n_replicates_done,
                        organ_meta=organ_meta,
                    )
                )
            if export_sharing_summaries:
                writers["sharing_aggregated_summaries"].write_rows(
                    _sharing_aggregate_rows(
                        set_id=set_id,
                        set_label=set_label,
                        cfg=cfg,
                        params=params_eff,
                        group_summary=group_summary,
                        n_replicates=n_replicates_done,
                    )
                )
            if export_vaf_count_spectra:
                writers["vaf_count_spectrum_aggregated_summaries"].write_rows(
                    _vaf_count_spectrum_aggregate_rows(
                        set_id=set_id,
                        set_label=set_label,
                        cfg=cfg,
                        params=params_eff,
                        group_summary=group_summary,
                        n_replicates=n_replicates_done,
                        organ_meta=organ_meta,
                    )
                )
            del group_summary, summary_acc
            gc.collect()
    finally:
        _close_streaming_writers(writers.values())

    manifest = {
        "mode": "grid_parameter",
        "output_mode": "streaming_bounded_memory",
        "n_parameter_sets_full": total_sets_full,
        "n_parameter_sets": total_sets_selected,
        "n_sim_per_set_full": int(total_replicates_full),
        "n_sim_per_set": int(n_sim_selected),
        "store_full_results": store_full,
        "full_results_layout": "one_json_per_parameter_set_and_replicate" if store_full else None,
        "export_raw_vafs": export_raw_vafs,
        "export_sharing_summaries": export_sharing_summaries,
        "export_replicate_sharing_summaries": export_replicate_sharing_summaries,
        "export_vaf_count_spectra": export_vaf_count_spectra,
        "export_legacy_binned_vaf_summaries": export_legacy_binned_vaf_summaries,
        "export_legacy_binned_vaf_replicates": export_legacy_binned_vaf_replicates,
        "topology_json_input": _path_str_for_output(str(cfg["topology_json"]), config_base_dir=config_path.parent, output_base_dir=grid_dir),
        "topology_json_checked": _relative_path_str(checked_topology_json, start=grid_dir),
        "topology_json_used": _path_str_for_output(str(cfg["topology_json"]), config_base_dir=config_path.parent, output_base_dir=grid_dir),
        "software_version": _software_version_info(pipeline_dir)["software_version"],
        "software_version_json": _relative_path_str(software_version_json, start=grid_dir),
    }
    if subset_info is not None:
        manifest["grid_subset"] = dict(subset_info)
    if replicate_subset_info is not None:
        manifest["replicate_subset"] = dict(replicate_subset_info)
        manifest["replicate_subset_rep_ids"] = {
            "start_rep": int(replicate_subset_info["start_index"]) - 1,
            "end_rep": int(replicate_subset_info["end_index"]) - 1,
        }
    _json_dump(manifest, grid_dir / "manifest.json")

    if not store_full:
        shutil.rmtree(tmp_results_dir, ignore_errors=True)

    print("Completed step: run")
    print(f"run_dir: {run_dir}")
    print(f"config_used_json: {config_used_json}")
    print(f"software_version_json: {software_version_json}")
    print(f"topology_json_input: {_path_str_for_output(str(cfg['topology_json']), config_base_dir=config_path.parent, output_base_dir=grid_dir)}")
    print(f"topology_json_checked: {_relative_path_str(checked_topology_json, start=grid_dir)}")
    print(f"grid_output_dir: {grid_dir}")
    print(f"parameter_sets_csv: {grid_dir / 'parameter_sets.csv'}")
    print(f"replicate_summaries_csv: {grid_dir / 'replicate_summaries.csv'}")
    print(f"aggregated_summaries_csv: {grid_dir / 'aggregated_summaries.csv'}")
    print(f"organ_replicate_summaries_csv: {grid_dir / 'organ_replicate_summaries.csv'}")
    print(f"organ_aggregated_summaries_csv: {grid_dir / 'organ_aggregated_summaries.csv'}")
    if export_legacy_binned_vaf_summaries:
        print(f"vaf_class_aggregated_summaries_csv: {grid_dir / 'vaf_class_aggregated_summaries.csv'}")
    if export_sharing_summaries:
        print(f"sharing_aggregated_summaries_csv: {grid_dir / 'sharing_aggregated_summaries.csv'}")
    if export_legacy_binned_vaf_replicates:
        print(f"vaf_class_replicate_summaries_csv_gz: {grid_dir / 'vaf_class_replicate_summaries.csv.gz'}")
    if export_replicate_sharing_summaries:
        print(f"sharing_replicate_summaries_csv_gz: {grid_dir / 'sharing_replicate_summaries.csv.gz'}")
    if export_vaf_count_spectra:
        print(f"vaf_count_spectrum_aggregated_summaries_csv: {grid_dir / 'vaf_count_spectrum_aggregated_summaries.csv'}")
        print(f"vaf_count_spectrum_replicate_summaries_csv_gz: {grid_dir / 'vaf_count_spectrum_replicate_summaries.csv.gz'}")
    if export_raw_vafs:
        print(f"raw_vafs_csv_gz: {grid_dir / 'raw_vafs.csv.gz'}")
    if subset_info is not None:
        print(f"grid_subset: {subset_info['start_index']}..{subset_info['end_index']} of {total_sets_full}")
    print(f"n_parameter_sets: {total_sets_selected}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run simSOMA experiments from a JSON config.")
    parser.add_argument("--config", type=Path, required=True, help="Path to JSON config file.")
    parser.add_argument(
        "--step",
        type=str,
        default="run",
        choices=["check", "run"],
        help="Which workflow step to execute. 'check' = inspect/convert/plot topology, 'run' = run simulation from the checked topology.",
    )
    args = parser.parse_args()

    config_path = args.config.expanduser().resolve()
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    pipeline_dir = Path(__file__).resolve().parent
    raw_cfg = _load_config(config_path)
    if _grouped_config_uses_simulation_kappa(raw_cfg):
        raise ValueError(
            "In grouped configs, do not place kappa_sr under simulation. "
            "Set topology.mapping_rate instead; the driver will use that value both for topology conversion and as the effective kappa_sr during simulation."
        )
    cfg = _normalize_config(raw_cfg)
    _validate_config(cfg)

    if args.step == "check":
        _run_check_step(cfg=cfg, raw_cfg=raw_cfg, config_path=config_path, pipeline_dir=pipeline_dir)
        return

    _run_grid_step(cfg=cfg, raw_cfg=raw_cfg, config_path=config_path, pipeline_dir=pipeline_dir)


if __name__ == "__main__":
    main()
