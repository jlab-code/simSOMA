from __future__ import annotations

import argparse
import copy
import csv
import gzip
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import run_from_config as rfc
import summaries

CSV_FILES = [
    "parameter_sets.csv",
    "replicate_summaries.csv",
    "aggregated_summaries.csv",
    "organ_replicate_summaries.csv",
    "organ_aggregated_summaries.csv",
]
SUMMARY_FEATURE_CSV_FILES = [
    "vaf_class_aggregated_summaries.csv",
    "sharing_aggregated_summaries.csv",
    "vaf_count_spectrum_aggregated_summaries.csv",
]
REPLICATE_FEATURE_GZ_FILES = [
    "vaf_class_replicate_summaries.csv.gz",
    "sharing_replicate_summaries.csv.gz",
    "vaf_count_spectrum_replicate_summaries.csv.gz",
]
RAW_VAF_FILE = "raw_vafs.csv.gz"
REALIZED_EVENT_FILE = "realized_event_truth.csv.gz"


def _resolve_path(value: str, *, base_dir: Path) -> Path:
    return rfc._resolve_path(value, base_dir=base_dir)


def _run_dir_from_cfg(cfg: Dict[str, Any], *, config_path: Path) -> Path:
    outdir_root = rfc._resolve_outdir_root(str(cfg["run"]["outdir_root"]), base_dir=config_path.parent)
    return outdir_root / str(cfg["run"]["experiment_name"])


def _partition_indices(total: int, n_splits: int) -> List[Tuple[int, int]]:
    if total < 1:
        raise ValueError("total must be >= 1")
    n_splits = max(1, min(int(n_splits), int(total)))
    base = total // n_splits
    rem = total % n_splits
    out: List[Tuple[int, int]] = []
    start = 1
    for i in range(n_splits):
        width = base + (1 if i < rem else 0)
        end = start + width - 1
        out.append((start, end))
        start = end + 1
    return out


def _child_experiment_name(base_name: str, split_idx: int, n_splits: int, *, split_axis: str) -> str:
    tag = "psplit" if split_axis == "parameter" else "rsplit"
    return f"{base_name}__{tag}_{split_idx:03d}_of_{n_splits:03d}"


def _relative_path_str(target: Path, *, start: Path) -> str:
    try:
        return Path(os.path.relpath(str(target), str(start))).as_posix()
    except Exception:
        return str(target)


def _retarget_path_for_child(value: str, *, source_base_dir: Path, child_base_dir: Path) -> str:
    resolved = _resolve_path(value, base_dir=source_base_dir)
    return _relative_path_str(resolved, start=child_base_dir)


def _write_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, sort_keys=True), encoding="utf-8")


def _verify_merge_outputs(merged_dir: Path) -> None:
    required = [merged_dir / name for name in CSV_FILES]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Merged outputs incomplete; refusing cleanup because these files are missing: " + ", ".join(missing)
        )


def _copy_master_check(master_run_dir: Path, child_run_dir: Path) -> None:
    src = master_run_dir / "topology_check"
    dst = child_run_dir / "topology_check"
    if not src.exists():
        raise FileNotFoundError(f"Master topology_check directory not found: {src}. Run the check step first.")
    if dst.exists():
        shutil.rmtree(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(src, dst)



def _concat_csv_stream(src_paths: Sequence[Path], dst_path: Path) -> bool:
    """Concatenate split CSV/CSV.GZ files without materializing rows in memory.

    Returns False when no source file exists. Header schemas are required to be
    identical, which is expected for parameter-split runs because each split is
    produced from the same master config and non-overlapping parameter sets.
    """
    existing = [p for p in src_paths if p.exists()]
    if not existing:
        return False
    dst_path.parent.mkdir(parents=True, exist_ok=True)
    opener_out = gzip.open if dst_path.suffix == ".gz" else open
    header: Optional[str] = None
    with opener_out(dst_path, "wt", newline="", encoding="utf-8") as out_fh:
        for src in existing:
            opener_in = gzip.open if src.suffix == ".gz" else open
            with opener_in(src, "rt", newline="", encoding="utf-8") as in_fh:
                src_header = in_fh.readline()
                if src_header == "":
                    continue
                if header is None:
                    header = src_header
                    out_fh.write(src_header)
                elif src_header != header:
                    raise ValueError(f"CSV header mismatch while merging {src} into {dst_path}")
                shutil.copyfileobj(in_fh, out_fh)
    return header is not None


def _merge_split_outputs_parameter_streaming(
    *,
    master_cfg: Dict[str, Any],
    master_config_path: Path,
    split_run_dirs: Sequence[Path],
    split_ranges: Sequence[Tuple[int, int]],
    split_axis: str,
) -> Path:
    master_run_dir = _run_dir_from_cfg(master_cfg, config_path=master_config_path)
    master_grid_dir = master_run_dir / "grid_parameter"
    if master_grid_dir.exists():
        raise FileExistsError(
            f"Merged grid output directory already exists: {master_grid_dir}. Remove or rename it before merging."
        )
    master_grid_dir.mkdir(parents=True, exist_ok=False)

    for split_run_dir in split_run_dirs:
        split_grid_dir = split_run_dir / "grid_parameter"
        if not split_grid_dir.exists():
            raise FileNotFoundError(f"Expected split grid output not found: {split_grid_dir}")

    for filename in CSV_FILES:
        _concat_csv_stream([d / "grid_parameter" / filename for d in split_run_dirs], master_grid_dir / filename)
    for filename in SUMMARY_FEATURE_CSV_FILES:
        _concat_csv_stream([d / "grid_parameter" / filename for d in split_run_dirs], master_grid_dir / filename)
    for filename in REPLICATE_FEATURE_GZ_FILES:
        _concat_csv_stream([d / "grid_parameter" / filename for d in split_run_dirs], master_grid_dir / filename)
    _concat_csv_stream([d / "grid_parameter" / RAW_VAF_FILE for d in split_run_dirs], master_grid_dir / RAW_VAF_FILE)
    _concat_csv_stream([d / "grid_parameter" / REALIZED_EVENT_FILE for d in split_run_dirs], master_grid_dir / REALIZED_EVENT_FILE)

    full_results_dst = master_grid_dir / "full_results"
    _merge_full_results_parameter(split_run_dirs=split_run_dirs, full_results_dst=full_results_dst)

    # For parameter splits, the number of replicates per set is already complete
    # in each split, so no statistical recombination is needed during merge.
    n_param = 0
    for start, end in split_ranges:
        n_param += int(end) - int(start) + 1
    _write_merged_manifest(
        master_cfg=master_cfg,
        master_config_path=master_config_path,
        master_grid_dir=master_grid_dir,
        master_run_dir=master_run_dir,
        split_ranges=split_ranges,
        split_run_dirs=split_run_dirs,
        split_axis=split_axis,
        n_parameter_sets=int(n_param),
        n_sim_total=int(master_cfg["simulation"].get("n_sim", 1)),
    )
    return master_grid_dir

def _read_csv_rows(path: Path) -> Tuple[List[str], List[Dict[str, Any]]]:
    if not path.exists():
        return [], []
    if path.suffix == ".gz":
        with gzip.open(path, "rt", newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            return list(reader.fieldnames or []), list(reader)
    with path.open("r", newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        return list(reader.fieldnames or []), list(reader)


def _write_csv(path: Path, rows: List[Dict[str, Any]], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        with gzip.open(path, "wt", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(fieldnames))
            writer.writeheader()
            for row in rows:
                writer.writerow({k: row.get(k, "") for k in fieldnames})
        return
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(fieldnames))
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def _sort_rows(rows: List[Dict[str, Any]], kind: str) -> List[Dict[str, Any]]:
    def set_num(row: Dict[str, Any]) -> int:
        sid = str(row.get("set_id", ""))
        try:
            return int(sid.split("_", 1)[1])
        except Exception:
            return 10**12

    def rep_num(row: Dict[str, Any]) -> int:
        try:
            return int(row.get("rep", 0))
        except Exception:
            return 0

    def mut_num(row: Dict[str, Any]) -> int:
        try:
            return int(row.get("mutation_id", 0))
        except Exception:
            return 0

    def organ_key(row: Dict[str, Any]) -> List[Any]:
        return rfc._natural_key(str(row.get("organ_id", "")))

    if kind == "parameter_sets":
        return sorted(rows, key=lambda r: (set_num(r), str(r.get("set_label", ""))))
    if kind == "aggregated_summaries":
        return sorted(rows, key=lambda r: (set_num(r), str(r.get("set_label", ""))))
    if kind == "replicate_summaries":
        return sorted(rows, key=lambda r: (set_num(r), rep_num(r)))
    if kind == "organ_replicate_summaries":
        return sorted(rows, key=lambda r: (set_num(r), organ_key(r), rep_num(r)))
    if kind == "organ_aggregated_summaries":
        return sorted(rows, key=lambda r: (set_num(r), organ_key(r)))
    if kind == "realized_event_truth":
        return sorted(rows, key=lambda r: (set_num(r), rep_num(r), str(r.get("event_type", "")), str(r.get("event_id", ""))))
    if kind == "vaf_class_aggregated_summaries":
        return sorted(rows, key=lambda r: (set_num(r), str(r.get("level", "")), organ_key(r), str(r.get("variant_class", "")), _as_int(r.get("bin_index"), 0)))
    if kind == "sharing_aggregated_summaries":
        return sorted(rows, key=lambda r: (set_num(r), str(r.get("statistic", "")), str(r.get("variant_class", "")), str(r.get("sharing_degree", "")), organ_key(r), str(r.get("organ_a", "")), str(r.get("organ_b", ""))))
    if kind == "vaf_count_spectrum_aggregated_summaries":
        return sorted(rows, key=lambda r: (set_num(r), str(r.get("level", "")), organ_key(r), str(r.get("variant_class", "")), str(r.get("sharing_degree", "")), _as_int(r.get("n_sampled_cells"), 0), _as_int(r.get("allele_count"), 0)))
    if kind == "vaf_class_replicate_summaries":
        return sorted(rows, key=lambda r: (set_num(r), rep_num(r), str(r.get("level", "")), organ_key(r), str(r.get("variant_class", "")), _as_int(r.get("bin_index"), 0)))
    if kind == "sharing_replicate_summaries":
        return sorted(rows, key=lambda r: (set_num(r), rep_num(r), str(r.get("statistic", "")), str(r.get("variant_class", "")), str(r.get("sharing_degree", "")), organ_key(r), str(r.get("organ_a", "")), str(r.get("organ_b", ""))))
    if kind == "vaf_count_spectrum_replicate_summaries":
        return sorted(rows, key=lambda r: (set_num(r), rep_num(r), str(r.get("level", "")), organ_key(r), str(r.get("variant_class", "")), str(r.get("sharing_degree", "")), _as_int(r.get("n_sampled_cells"), 0), _as_int(r.get("allele_count"), 0)))
    if kind == "raw_vafs":
        return sorted(rows, key=lambda r: (set_num(r), organ_key(r), rep_num(r), mut_num(r)))
    return rows


def _as_int(value: Any, default: int = 0) -> int:
    try:
        if value in (None, ""):
            return default
        return int(float(value))
    except Exception:
        return default


def _as_float(value: Any) -> Optional[float]:
    try:
        if value in (None, ""):
            return None
        return float(value)
    except Exception:
        return None


def _combine_mean_sd_stats(stat_rows: Sequence[Dict[str, Any]]) -> Dict[str, Any]:
    if not stat_rows:
        return {}
    out = dict(stat_rows[0])
    n_values = [_as_int(r.get("n_replicates"), 0) for r in stat_rows]
    total_n = sum(n_values)
    out["n_replicates"] = int(total_n)

    mean_cols = [c for c in out.keys() if c.endswith("_mean")]
    for mean_col in mean_cols:
        base = mean_col[:-5]
        sd_col = f"{base}_sd"
        if sd_col not in out:
            continue
        groups: List[Tuple[int, float, float]] = []
        for row, n in zip(stat_rows, n_values):
            mu = _as_float(row.get(mean_col))
            sd = _as_float(row.get(sd_col))
            if mu is None or sd is None or n <= 0:
                continue
            groups.append((n, mu, sd))
        if not groups:
            out[mean_col] = ""
            out[sd_col] = ""
            continue
        total = sum(n for n, _, _ in groups)
        mean = sum(n * mu for n, mu, _ in groups) / float(total)
        var = sum(n * ((sd ** 2) + (mu - mean) ** 2) for n, mu, sd in groups) / float(total)
        out[mean_col] = mean
        out[sd_col] = var ** 0.5
    return out


def _dedupe_parameter_rows(rows: Sequence[Dict[str, Any]], *, n_sim_total: int) -> List[Dict[str, Any]]:
    by_set: Dict[str, Dict[str, Any]] = {}
    for row in rows:
        sid = str(row.get("set_id", ""))
        if sid not in by_set:
            by_set[sid] = dict(row)
            by_set[sid]["n_sim"] = int(n_sim_total)
    return _sort_rows(list(by_set.values()), "parameter_sets")


def _ensure_unique(rows: Sequence[Dict[str, Any]], *, kind: str, key_fields: Sequence[str]) -> None:
    seen = set()
    for row in rows:
        key = tuple(str(row.get(k, "")) for k in key_fields)
        if key in seen:
            raise ValueError(f"Duplicate {kind} row during merge for key={key}")
        seen.add(key)


def _combine_aggregate_rows(rows: Sequence[Dict[str, Any]], *, kind: str) -> List[Dict[str, Any]]:
    groups: Dict[Tuple[str, ...], List[Dict[str, Any]]] = {}
    if kind == "aggregated_summaries":
        key_fields = ("set_id",)
    elif kind == "organ_aggregated_summaries":
        key_fields = ("set_id", "organ_id")
    elif kind == "vaf_class_aggregated_summaries":
        key_fields = ("set_id", "level", "organ_id", "variant_class", "bin_index")
    elif kind == "sharing_aggregated_summaries":
        key_fields = ("set_id", "statistic", "variant_class", "sharing_degree", "organ_id", "organ_a", "organ_b")
    elif kind == "vaf_count_spectrum_aggregated_summaries":
        key_fields = ("set_id", "level", "organ_id", "variant_class", "sharing_degree", "n_sampled_cells", "allele_count")
    else:
        raise ValueError(f"Unsupported kind: {kind}")
    for row in rows:
        key = tuple(str(row.get(k, "")) for k in key_fields)
        groups.setdefault(key, []).append(row)
    merged = [_combine_mean_sd_stats(group_rows) for _, group_rows in sorted(groups.items())]
    return _sort_rows(merged, kind)


def _merge_full_results_parameter(*, split_run_dirs: Sequence[Path], full_results_dst: Path) -> bool:
    saw_full_results = False
    for split_run_dir in split_run_dirs:
        split_full_results = split_run_dir / "grid_parameter" / "full_results"
        if not split_full_results.exists():
            continue
        saw_full_results = True
        full_results_dst.mkdir(parents=True, exist_ok=True)
        for src in sorted(split_full_results.glob("*.json")):
            dst = full_results_dst / src.name
            if dst.exists():
                raise FileExistsError(f"Duplicate full-result file during merge: {dst}")
            shutil.copy2(src, dst)
    return saw_full_results


def _merge_full_results_replicate(*, split_run_dirs: Sequence[Path], full_results_dst: Path) -> bool:
    grouped: Dict[str, List[Path]] = {}
    for split_run_dir in split_run_dirs:
        split_full_results = split_run_dir / "grid_parameter" / "full_results"
        if not split_full_results.exists():
            continue
        for src in sorted(split_full_results.glob("*.json")):
            grouped.setdefault(src.name, []).append(src)
    if not grouped:
        return False
    full_results_dst.mkdir(parents=True, exist_ok=True)
    for name, paths in grouped.items():
        merged_payload: Optional[Dict[str, Any]] = None
        merged_sims: List[Dict[str, Any]] = []
        for src in paths:
            payload = json.loads(src.read_text(encoding="utf-8"))
            runs = payload.get("runs", [])
            if len(runs) != 1:
                raise RuntimeError(f"Unexpected full-results structure in {src}: expected exactly one run entry.")
            run_entry = runs[0]
            sims = list(run_entry.get("sims", []))
            if merged_payload is None:
                merged_payload = payload
                merged_payload["runs"][0]["sims"] = []
            merged_sims.extend(sims)
        assert merged_payload is not None
        merged_sims.sort(key=lambda sim: _as_int(sim.get("rep"), 0))
        merged_payload["runs"][0]["sims"] = merged_sims
        rep_summaries = [sim.get("summary") for sim in merged_sims if isinstance(sim.get("summary"), dict)]
        if len(rep_summaries) == len(merged_sims) and rep_summaries:
            merged_payload["runs"][0]["summary_aggregate"] = summaries.aggregate_replicate_summaries(rep_summaries)
        merged_payload["runs"][0].setdefault("params", {})["n_sim"] = len(merged_sims)
        (full_results_dst / name).write_text(json.dumps(merged_payload, indent=2, sort_keys=True), encoding="utf-8")
    return True


def _write_merged_manifest(
    *,
    master_cfg: Dict[str, Any],
    master_config_path: Path,
    master_grid_dir: Path,
    master_run_dir: Path,
    split_ranges: Sequence[Tuple[int, int]],
    split_run_dirs: Sequence[Path],
    split_axis: str,
    n_parameter_sets: int,
    n_sim_total: int,
) -> None:
    if split_axis == "parameter":
        split_desc = [
            {
                "split_index": i + 1,
                "start_index": int(start),
                "end_index": int(end),
                "experiment_name": split_run_dirs[i].name,
            }
            for i, (start, end) in enumerate(split_ranges)
        ]
    else:
        split_desc = [
            {
                "split_index": i + 1,
                "start_rep": int(start) - 1,
                "end_rep": int(end) - 1,
                "start_rep_index": int(start),
                "end_rep_index": int(end),
                "experiment_name": split_run_dirs[i].name,
            }
            for i, (start, end) in enumerate(split_ranges)
        ]
    manifest = {
        "mode": "grid_parameter",
        "merge_mode": "split_subconfigs",
        "split_axis": split_axis,
        "n_parameter_sets_full": rfc._grid_total_n_sets(master_cfg["simulation"]["modules"]),
        "n_parameter_sets": int(n_parameter_sets),
        "n_sim_per_set": int(n_sim_total),
        "store_full_results": bool(master_cfg["simulation"].get("store_full_results", False)),
        "export_raw_vafs": bool(master_cfg["simulation"].get("export_raw_vafs", False)),
        "topology_json_input": _relative_path_str(_resolve_path(str(master_cfg["topology"]["topology_json"]), base_dir=master_config_path.parent), start=master_grid_dir),
        "topology_json_checked": _relative_path_str(master_run_dir / "topology_check" / "topology_internal_steps.json", start=master_grid_dir),
        "topology_json_used": _relative_path_str(_resolve_path(str(master_cfg["topology"]["topology_json"]), base_dir=master_config_path.parent), start=master_grid_dir),
        "split_ranges": split_desc,
    }
    _write_json(master_grid_dir / "manifest.json", manifest)


def _merge_split_outputs(
    *,
    master_cfg: Dict[str, Any],
    master_config_path: Path,
    split_run_dirs: Sequence[Path],
    split_ranges: Sequence[Tuple[int, int]],
    split_axis: str,
) -> Path:
    if split_axis == "parameter":
        return _merge_split_outputs_parameter_streaming(
            master_cfg=master_cfg,
            master_config_path=master_config_path,
            split_run_dirs=split_run_dirs,
            split_ranges=split_ranges,
            split_axis=split_axis,
        )

    master_run_dir = _run_dir_from_cfg(master_cfg, config_path=master_config_path)
    master_grid_dir = master_run_dir / "grid_parameter"
    if master_grid_dir.exists():
        raise FileExistsError(
            f"Merged grid output directory already exists: {master_grid_dir}. Remove or rename it before merging."
        )
    master_grid_dir.mkdir(parents=True, exist_ok=False)

    merged_rows: Dict[str, List[Dict[str, Any]]] = {
        "parameter_sets": [],
        "replicate_summaries": [],
        "aggregated_summaries": [],
        "organ_replicate_summaries": [],
        "organ_aggregated_summaries": [],
        "vaf_class_aggregated_summaries": [],
        "sharing_aggregated_summaries": [],
        "vaf_count_spectrum_aggregated_summaries": [],
        "vaf_class_replicate_summaries": [],
        "sharing_replicate_summaries": [],
        "vaf_count_spectrum_replicate_summaries": [],
        "raw_vafs": [],
        "realized_event_truth": [],
    }
    fieldnames: Dict[str, List[str]] = {k: [] for k in merged_rows}
    saw_raw_vafs = False
    saw_realized_event_truth = False
    saw_optional: Dict[str, bool] = {k: False for k in merged_rows}

    for split_run_dir in split_run_dirs:
        split_grid_dir = split_run_dir / "grid_parameter"
        if not split_grid_dir.exists():
            raise FileNotFoundError(f"Expected split grid output not found: {split_grid_dir}")
        for filename in CSV_FILES:
            kind = filename.replace(".csv", "")
            fns, rows = _read_csv_rows(split_grid_dir / filename)
            if rows:
                merged_rows[kind].extend(rows)
            for fn in fns:
                if fn not in fieldnames[kind]:
                    fieldnames[kind].append(fn)
        for filename in SUMMARY_FEATURE_CSV_FILES:
            kind = filename.replace(".csv", "")
            fns, rows = _read_csv_rows(split_grid_dir / filename)
            if rows:
                saw_optional[kind] = True
                merged_rows[kind].extend(rows)
            for fn in fns:
                if fn not in fieldnames[kind]:
                    fieldnames[kind].append(fn)
        for filename in REPLICATE_FEATURE_GZ_FILES:
            kind = filename.replace(".csv.gz", "")
            fns, rows = _read_csv_rows(split_grid_dir / filename)
            if rows:
                saw_optional[kind] = True
                merged_rows[kind].extend(rows)
            for fn in fns:
                if fn not in fieldnames[kind]:
                    fieldnames[kind].append(fn)
        raw_fns, raw_rows = _read_csv_rows(split_grid_dir / RAW_VAF_FILE)
        if raw_rows:
            saw_raw_vafs = True
            merged_rows["raw_vafs"].extend(raw_rows)
        for fn in raw_fns:
            if fn not in fieldnames["raw_vafs"]:
                fieldnames["raw_vafs"].append(fn)
        event_fns, event_rows = _read_csv_rows(split_grid_dir / REALIZED_EVENT_FILE)
        if event_rows:
            saw_realized_event_truth = True
            merged_rows["realized_event_truth"].extend(event_rows)
        for fn in event_fns:
            if fn not in fieldnames["realized_event_truth"]:
                fieldnames["realized_event_truth"].append(fn)

    n_sim_total = int(master_cfg["simulation"].get("n_sim", 1))
    param_rows = _dedupe_parameter_rows(merged_rows["parameter_sets"], n_sim_total=n_sim_total)
    repl_rows = _sort_rows(merged_rows["replicate_summaries"], "replicate_summaries")
    organ_repl_rows = _sort_rows(merged_rows["organ_replicate_summaries"], "organ_replicate_summaries")
    _ensure_unique(repl_rows, kind="replicate_summaries", key_fields=("set_id", "rep"))
    _ensure_unique(organ_repl_rows, kind="organ_replicate_summaries", key_fields=("set_id", "organ_id", "rep"))
    if saw_raw_vafs:
        raw_rows = _sort_rows(merged_rows["raw_vafs"], "raw_vafs")
        _ensure_unique(raw_rows, kind="raw_vafs", key_fields=("set_id", "organ_id", "rep", "mutation_id"))
    else:
        raw_rows = []
    if saw_realized_event_truth:
        realized_event_rows = _sort_rows(merged_rows["realized_event_truth"], "realized_event_truth")
        _ensure_unique(realized_event_rows, kind="realized_event_truth", key_fields=("set_id", "rep", "event_type", "event_id"))
    else:
        realized_event_rows = []

    if split_axis == "parameter":
        agg_rows = _sort_rows(merged_rows["aggregated_summaries"], "aggregated_summaries")
        organ_agg_rows = _sort_rows(merged_rows["organ_aggregated_summaries"], "organ_aggregated_summaries")
        vaf_class_agg_rows = _sort_rows(merged_rows["vaf_class_aggregated_summaries"], "vaf_class_aggregated_summaries")
        sharing_agg_rows = _sort_rows(merged_rows["sharing_aggregated_summaries"], "sharing_aggregated_summaries")
        vaf_count_agg_rows = _sort_rows(merged_rows["vaf_count_spectrum_aggregated_summaries"], "vaf_count_spectrum_aggregated_summaries")
    else:
        agg_rows = _combine_aggregate_rows(merged_rows["aggregated_summaries"], kind="aggregated_summaries")
        organ_agg_rows = _combine_aggregate_rows(merged_rows["organ_aggregated_summaries"], kind="organ_aggregated_summaries")
        vaf_class_agg_rows = _combine_aggregate_rows(merged_rows["vaf_class_aggregated_summaries"], kind="vaf_class_aggregated_summaries") if saw_optional.get("vaf_class_aggregated_summaries") else []
        sharing_agg_rows = _combine_aggregate_rows(merged_rows["sharing_aggregated_summaries"], kind="sharing_aggregated_summaries") if saw_optional.get("sharing_aggregated_summaries") else []
        vaf_count_agg_rows = _combine_aggregate_rows(merged_rows["vaf_count_spectrum_aggregated_summaries"], kind="vaf_count_spectrum_aggregated_summaries") if saw_optional.get("vaf_count_spectrum_aggregated_summaries") else []

    vaf_class_repl_rows = _sort_rows(merged_rows["vaf_class_replicate_summaries"], "vaf_class_replicate_summaries")
    sharing_repl_rows = _sort_rows(merged_rows["sharing_replicate_summaries"], "sharing_replicate_summaries")
    vaf_count_repl_rows = _sort_rows(merged_rows["vaf_count_spectrum_replicate_summaries"], "vaf_count_spectrum_replicate_summaries")

    output_map = {
        "parameter_sets": param_rows,
        "replicate_summaries": repl_rows,
        "aggregated_summaries": agg_rows,
        "organ_replicate_summaries": organ_repl_rows,
        "organ_aggregated_summaries": organ_agg_rows,
    }
    for kind, rows in output_map.items():
        ordered = rfc._ordered_fieldnames(kind, rows) if rows else fieldnames[kind]
        _write_csv(master_grid_dir / f"{kind}.csv", rows, ordered)

    feature_outputs = {
        "vaf_class_aggregated_summaries": (vaf_class_agg_rows, "vaf_class_summaries", "vaf_class_aggregated_summaries.csv"),
        "sharing_aggregated_summaries": (sharing_agg_rows, "sharing_summaries", "sharing_aggregated_summaries.csv"),
        "vaf_count_spectrum_aggregated_summaries": (vaf_count_agg_rows, "vaf_count_spectrum_summaries", "vaf_count_spectrum_aggregated_summaries.csv"),
        "vaf_class_replicate_summaries": (vaf_class_repl_rows, "vaf_class_summaries", "vaf_class_replicate_summaries.csv.gz"),
        "sharing_replicate_summaries": (sharing_repl_rows, "sharing_summaries", "sharing_replicate_summaries.csv.gz"),
        "vaf_count_spectrum_replicate_summaries": (vaf_count_repl_rows, "vaf_count_spectrum_summaries", "vaf_count_spectrum_replicate_summaries.csv.gz"),
    }
    for kind, (rows, rfc_kind, filename) in feature_outputs.items():
        if not saw_optional.get(kind) and not rows:
            continue
        ordered = rfc._ordered_fieldnames(rfc_kind, rows) if rows else fieldnames[kind]
        _write_csv(master_grid_dir / filename, rows, ordered)

    if saw_raw_vafs:
        ordered = rfc._ordered_fieldnames("raw_vafs", raw_rows) if raw_rows else fieldnames["raw_vafs"]
        _write_csv(master_grid_dir / RAW_VAF_FILE, raw_rows, ordered)
    if saw_realized_event_truth:
        ordered = fieldnames["realized_event_truth"]
        _write_csv(master_grid_dir / REALIZED_EVENT_FILE, realized_event_rows, ordered)

    full_results_dst = master_grid_dir / "full_results"
    if split_axis == "parameter":
        _merge_full_results_parameter(split_run_dirs=split_run_dirs, full_results_dst=full_results_dst)
    else:
        _merge_full_results_replicate(split_run_dirs=split_run_dirs, full_results_dst=full_results_dst)

    _write_merged_manifest(
        master_cfg=master_cfg,
        master_config_path=master_config_path,
        master_grid_dir=master_grid_dir,
        master_run_dir=master_run_dir,
        split_ranges=split_ranges,
        split_run_dirs=split_run_dirs,
        split_axis=split_axis,
        n_parameter_sets=len(param_rows),
        n_sim_total=n_sim_total,
    )
    return master_grid_dir


def _cleanup_split_workspace(*, split_root: Path, split_run_dirs: Sequence[Path], merged_dir: Path) -> None:
    _verify_merge_outputs(merged_dir)

    plan_src = split_root / "plan.json"
    if plan_src.exists():
        shutil.copy2(plan_src, merged_dir / "split_plan.json")

    for run_dir in split_run_dirs:
        if run_dir.exists():
            shutil.rmtree(run_dir)

    if split_root.exists():
        shutil.rmtree(split_root)


def _launch_commands(commands: Sequence[List[str]], *, max_parallel: int, cwd: Path) -> None:
    max_parallel = max(1, int(max_parallel))
    running: List[Tuple[subprocess.Popen[str], List[str]]] = []
    pending = list(commands)

    while pending or running:
        while pending and len(running) < max_parallel:
            cmd = pending.pop(0)
            proc = subprocess.Popen(cmd, cwd=str(cwd), text=True)
            running.append((proc, cmd))

        next_running: List[Tuple[subprocess.Popen[str], List[str]]] = []
        for proc, cmd in running:
            ret = proc.poll()
            if ret is None:
                next_running.append((proc, cmd))
                continue
            if ret != 0:
                for other, _ in next_running:
                    if other.poll() is None:
                        other.terminate()
                raise RuntimeError(f"Split run failed with exit code {ret}: {' '.join(cmd)}")
        if len(next_running) == len(running):
            proc, cmd = next_running[0]
            ret = proc.wait()
            if ret != 0:
                for other, _ in next_running[1:]:
                    if other.poll() is None:
                        other.terminate()
                raise RuntimeError(f"Split run failed with exit code {ret}: {' '.join(cmd)}")
            next_running = [(p, c) for p, c in next_running[1:] if p.poll() is None]
        running = next_running


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description="Split, launch, and merge simSOMA grid runs from one master config.")
    parser.add_argument("--master-config", type=Path, required=True, help="Grouped grid_parameter master config JSON.")
    parser.add_argument("--n-splits", type=int, default=1, help="How many subconfigs to generate.")
    parser.add_argument("--jobs", type=int, default=1, help="How many split runs to launch concurrently.")
    parser.add_argument("--split-axis", choices=["parameter", "replicate"], default="parameter", help="Split over parameter sets or over replicate indices.")
    parser.add_argument("--skip-check", action="store_true", help="Do not run the master topology check step.")
    parser.add_argument("--dry-run", action="store_true", help="Write subconfigs and print the planned ranges without launching runs.")
    parser.add_argument("--cleanup-splits", action="store_true", help="Remove grid_splits workspace after a successful merge.")
    args = parser.parse_args(argv)

    master_config_path = args.master_config.expanduser().resolve()
    if not master_config_path.exists():
        raise FileNotFoundError(f"Master config not found: {master_config_path}")

    pipeline_dir = Path(__file__).resolve().parent
    raw_cfg = json.loads(master_config_path.read_text(encoding="utf-8"))
    cfg = rfc._normalize_config(raw_cfg)
    rfc._validate_config(cfg)

    total_sets = rfc._grid_total_n_sets(cfg["modules"])
    total_replicates = int(cfg["n_sim"])
    total_units = total_sets if args.split_axis == "parameter" else total_replicates
    split_ranges = _partition_indices(total_units, args.n_splits)
    n_splits_eff = len(split_ranges)

    master_run_dir = _run_dir_from_cfg(raw_cfg, config_path=master_config_path)
    split_root = master_run_dir / "grid_splits"
    subconfig_dir = split_root / "subconfigs"
    split_root.mkdir(parents=True, exist_ok=True)
    subconfig_dir.mkdir(parents=True, exist_ok=True)

    master_check_dir = master_run_dir / "topology_check"
    if not args.skip_check:
        if master_check_dir.exists():
            print(f"Reusing existing master topology check: {master_check_dir}")
        else:
            subprocess.run(
                [sys.executable, str(pipeline_dir / "run_from_config.py"), "--config", str(master_config_path), "--step", "check"],
                check=True,
                cwd=str(pipeline_dir),
            )

    commands: List[List[str]] = []
    split_run_dirs: List[Path] = []

    for split_idx, (start_idx, end_idx) in enumerate(split_ranges, start=1):
        child_cfg = copy.deepcopy(raw_cfg)
        base_name = str(raw_cfg["run"]["experiment_name"])
        child_name = _child_experiment_name(base_name, split_idx, n_splits_eff, split_axis=args.split_axis)
        child_cfg["run"]["experiment_name"] = child_name
        child_config_path = subconfig_dir / f"{child_name}.json"
        # Keep project-root-style paths unchanged. The child process inherits
        # SIMSOMA_HOME, SIMSOMA_OUTPUT_DIR, and SIMSOMA_INPUT_DIR from 00_pipeline.sh,
        # so paths such as simSOMA_output and simSOMA_inputs/... remain unambiguous
        # even though child configs live inside grid_splits/subconfigs/.
        child_cfg["run"]["outdir_root"] = str(raw_cfg["run"]["outdir_root"])
        child_cfg["topology"]["topology_json"] = str(raw_cfg["topology"]["topology_json"])
        child_cfg.setdefault("simulation", {})
        master_observation = raw_cfg.get("observation_model")
        if isinstance(master_observation, dict):
            # Observation post-processing occurs once after the split outputs are merged.
            child_cfg.pop("observation_model", None)
            if bool(master_observation.get("export_fitsoma", False)):
                child_cfg["simulation"]["export_raw_vafs"] = True
                child_cfg["simulation"]["export_realized_event_truth"] = True
        if args.split_axis == "parameter":
            child_cfg["simulation"]["grid_subset"] = {
                "start_index": int(start_idx),
                "end_index": int(end_idx),
            }
            child_cfg["simulation"].pop("replicate_subset", None)
        else:
            child_cfg["simulation"]["replicate_subset"] = {
                "start_index": int(start_idx),
                "end_index": int(end_idx),
            }
            child_cfg["simulation"].pop("grid_subset", None)
        child_cfg["simulation"]["split_metadata"] = {
            "split_index": int(split_idx),
            "n_splits": int(n_splits_eff),
            "split_axis": args.split_axis,
        }
        _write_json(child_config_path, child_cfg)

        child_run_dir = _run_dir_from_cfg(child_cfg, config_path=child_config_path)
        split_run_dirs.append(child_run_dir)
        _copy_master_check(master_run_dir, child_run_dir)

        commands.append([
            sys.executable,
            str(pipeline_dir / "run_from_config.py"),
            "--config",
            str(child_config_path),
            "--step",
            "run",
        ])

    if args.split_axis == "parameter":
        split_desc = [
            {
                "split_index": i + 1,
                "start_index": int(start),
                "end_index": int(end),
                "n_parameter_sets": int(end - start + 1),
                "subconfig": _relative_path_str(
                    subconfig_dir / f"{_child_experiment_name(str(raw_cfg['run']['experiment_name']), i + 1, n_splits_eff, split_axis=args.split_axis)}.json",
                    start=split_root,
                ),
                "split_run_dir": _relative_path_str(split_run_dirs[i], start=split_root),
            }
            for i, (start, end) in enumerate(split_ranges)
        ]
    else:
        split_desc = [
            {
                "split_index": i + 1,
                "start_rep": int(start) - 1,
                "end_rep": int(end) - 1,
                "start_rep_index": int(start),
                "end_rep_index": int(end),
                "n_replicates": int(end - start + 1),
                "subconfig": _relative_path_str(
                    subconfig_dir / f"{_child_experiment_name(str(raw_cfg['run']['experiment_name']), i + 1, n_splits_eff, split_axis=args.split_axis)}.json",
                    start=split_root,
                ),
                "split_run_dir": _relative_path_str(split_run_dirs[i], start=split_root),
            }
            for i, (start, end) in enumerate(split_ranges)
        ]

    plan = {
        "master_config": _relative_path_str(master_config_path, start=split_root),
        "master_run_dir": ".",
        "subconfig_dir": _relative_path_str(subconfig_dir, start=split_root),
        "split_axis": args.split_axis,
        "n_parameter_sets_full": int(total_sets),
        "n_sim_full": int(total_replicates),
        "n_splits": int(n_splits_eff),
        "jobs": int(max(1, args.jobs)),
        "cleanup_splits": bool(args.cleanup_splits),
        "split_ranges": split_desc,
    }
    _write_json(split_root / "plan.json", plan)

    print(f"Master config: {master_config_path}")
    print(f"Split axis: {args.split_axis}")
    print(f"Total parameter sets: {total_sets}")
    print(f"Total replicates per set: {total_replicates}")
    print(f"Generated {n_splits_eff} split config(s) under: {subconfig_dir}")
    print(f"Cleanup split workspace after merge: {bool(args.cleanup_splits)}")
    for i, (start_idx, end_idx) in enumerate(split_ranges, start=1):
        if args.split_axis == "parameter":
            print(f"  split {i}: set_{start_idx:05d} .. set_{end_idx:05d} ({end_idx - start_idx + 1} sets)")
        else:
            print(f"  split {i}: rep_{start_idx:05d} .. rep_{end_idx:05d} ({end_idx - start_idx + 1} replicates)")

    if args.dry_run:
        print("Dry run only; no split runs launched.")
        return

    _launch_commands(commands, max_parallel=max(1, args.jobs), cwd=pipeline_dir)
    merged_dir = _merge_split_outputs(
        master_cfg=raw_cfg,
        master_config_path=master_config_path,
        split_run_dirs=split_run_dirs,
        split_ranges=split_ranges,
        split_axis=args.split_axis,
    )
    print(f"Merged outputs written to: {merged_dir}")

    if raw_cfg.get("observation_model") is not None:
        from configured_observation_model import run_configured_observation_model
        report = run_configured_observation_model(merged_dir.parent, raw_cfg)
        if report is not None:
            print(f"Observation-model outputs written to: {report['scenario_dir']}")

    if args.cleanup_splits:
        _cleanup_split_workspace(split_root=split_root, split_run_dirs=split_run_dirs, merged_dir=merged_dir)
        print(f"Removed split workspace: {split_root}")


if __name__ == "__main__":
    main()
