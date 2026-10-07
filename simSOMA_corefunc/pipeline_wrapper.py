"""pipeline_wrapper.py

Plumbing-only orchestrator.

Responsibilities
----------------
- Read a topology spec dict (from topology_io).
- Traverse branches using an explicit DFS stack.
- For each branch, simulate self-renewal sequentially between events.
- At an event time t_local (in self-renewal steps), run:
    * pre-branching amplification at fixed t_age=event_age
    * branching founder sampling OR organ formation

Critical time semantics
-----------------------
We track time as (t_age, t_dev) conceptually:
- t_age advances only with self-renewal steps along a branch
- t_dev advances only within fast developmental bursts (pre-branching, organ expansion)
  and does not consume t_age

Outputs
-------
- organ_events retain raw per-organ allele_counts_by_mutation so VAF spectra are derived downstream in summaries.py.
- branch_competition records realized branch-local displacement-bias values and
  favored seed positions used at branch start.
- victim choice uses one bounded locality parameter `victim_locality` in [0,1].
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional
from dataclasses import replace
import argparse
import json
import os
from pathlib import Path
import random
import hashlib
import itertools
import math
import copy

import numpy as np

import self_renewal
import pre_branching
import branching
import organ
import topology_io
import summaries
import branch_bias
import phyllotaxy


def _build_children_map(topo: Dict[str, Any]) -> Dict[str, List[str]]:
    children: Dict[str, List[str]] = {}
    for bid, b in topo["branches"].items():
        pid = b.get("parent_id")
        if pid is None:
            continue
        children.setdefault(pid, []).append(bid)
    return children


def _stable_seed_int(*parts: Any) -> int:
    h = hashlib.blake2b(digest_size=16)
    for p in parts:
        h.update(str(p).encode("utf-8"))
        h.update(b"\0")
    return int.from_bytes(h.digest(), "big") & ((1 << 63) - 1)


def _event_sort_key(e: Dict[str, Any]) -> tuple:
    t = int(e["time"])
    et = str(e["type"]).upper()
    pri = 0 if et == "BRANCH" else 1
    tgt = str(e.get("target_id", ""))
    return (t, pri, tgt)


def validate_topology_contract(topo: Dict[str, Any]) -> None:
    if not isinstance(topo, dict):
        raise TypeError("topo must be a dict")
    if "root_id" not in topo or "branches" not in topo:
        raise KeyError("topo must contain keys: 'root_id', 'branches'")
    branches = topo["branches"]
    if not isinstance(branches, dict) or len(branches) == 0:
        raise ValueError("topo['branches'] must be a non-empty dict")
    root = topo["root_id"]
    if root not in branches:
        raise KeyError(f"root_id {root} not present in topo['branches']")

    for bid, spec in branches.items():
        if "T" not in spec:
            raise KeyError(f"Branch {bid} missing 'T'")
        T = int(spec["T"])
        if T < 0:
            raise ValueError(f"Branch {bid} has negative T={T}")
        pid = spec.get("parent_id")
        if pid is None:
            if bid != root:
                raise ValueError(f"Non-root branch {bid} has parent_id=None")
        else:
            if pid not in branches:
                raise KeyError(f"Branch {bid} parent_id {pid} not present in branches")

        for e in list(spec.get("events", [])):
            if "time" not in e or "type" not in e or "target_id" not in e:
                raise KeyError(f"Branch {bid} has malformed event (needs time/type/target_id): {e}")
            t = int(e["time"])
            if t < 0 or t > T:
                raise ValueError(f"Branch {bid} event time {t} out of bounds [0,{T}]")
            et = str(e["type"]).upper()
            if et not in {"BRANCH", "ORGAN"}:
                raise ValueError(f"Branch {bid} has unknown event type {et}")
            if et == "BRANCH":
                child = str(e["target_id"])
                if child not in branches:
                    raise KeyError(f"Branch {bid} BRANCH targets missing child {child}")
                if str(branches[child].get("parent_id")) != str(bid):
                    raise ValueError(
                        f"BRANCH event mismatch: parent {bid} targets {child}, "
                        f"but child.parent_id={branches[child].get('parent_id')}"
                    )


def summarize_topology(topo: Dict[str, Any]) -> Dict[str, Any]:
    branches = topo.get("branches", {})
    n_branches = len(branches)
    n_events = 0
    n_branch_events = 0
    n_organ_events = 0
    organ_ids = set()

    for bid, spec in branches.items():
        evs = list(spec.get("events", []))
        n_events += len(evs)
        for e in evs:
            et = str(e.get("type", "")).upper()
            if et == "BRANCH":
                n_branch_events += 1
            elif et == "ORGAN":
                n_organ_events += 1
                organ_ids.add(str(e.get("target_id")))

    root = topo.get("root_id")
    children = _build_children_map(topo) if isinstance(branches, dict) else {}
    max_depth = 0
    stack = [(root, 0)] if root is not None else []
    visited = set()
    while stack:
        node, d = stack.pop()
        if node in visited:
            continue
        visited.add(node)
        max_depth = max(max_depth, d)
        for ch in children.get(node, []):
            stack.append((ch, d + 1))

    Ts = []
    for bid, spec in branches.items():
        try:
            Ts.append(int(spec.get("T")))
        except Exception:
            pass
    T_min = min(Ts) if Ts else None
    T_max = max(Ts) if Ts else None
    T_sum = sum(Ts) if Ts else None

    return {
        "n_branches": n_branches,
        "n_events": n_events,
        "n_branch_events": n_branch_events,
        "n_organ_events": n_organ_events,
        "n_unique_organs": len(organ_ids),
        "max_depth": max_depth,
        "T_min": T_min,
        "T_max": T_max,
        "T_sum": T_sum,
    }


def run_pipeline(
    topo: Dict[str, Any],
    sr_params: self_renewal.SelfRenewalParams,
    *,
    sam_boundary_cells: int,
    branch_precursor_number: int,
    organ_precursor_number: int,
    organ_total_cells: int,
    sequenced_cells: Optional[int],
    seed: Optional[int] = None,
    branch_comp_params: Optional[Dict[str, Any]] = None,
    module_paths: Optional[Dict[str, str]] = None,
    phyllotaxy_params: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    validate_topology_contract(topo)

    base_seed = int(seed) if seed is not None else (int.from_bytes(os.urandom(8), "big") & ((1 << 63) - 1))

    phyllo_summary = dict(topo.get("phyllotaxy", {}))
    if phyllotaxy_params is not None and not phyllo_summary.get("annotated", False):
        topo, phyllo_summary = topology_io.annotate_topology_phyllotaxy(
            topo,
            phyllotaxy_config=phyllotaxy_params,
            seed=_stable_seed_int(base_seed, "phyllotaxy"),
        )

    import module_registry
    from plugin_api import ModuleRegistry as _MR
    defaults = module_registry.default_registry()
    reg = _MR.from_paths(defaults, module_paths=module_paths)

    if branch_comp_params is None:
        branch_comp_params = {"bias_mode": "fixed", "fixed_value": 0.0}

    map_values = None
    if branch_comp_params.get("bias_mode") == "map":
        if branch_comp_params.get("map_values") is not None:
            map_values = {str(k): float(v) for k, v in dict(branch_comp_params["map_values"]).items()}
        elif branch_comp_params.get("map_json") is not None:
            map_values = branch_bias.BranchBiasState.load_map_json(branch_comp_params["map_json"])
        else:
            raise ValueError("bias_mode='map' requires map_values or map_json")

    comp_seed = _stable_seed_int(base_seed, "branch_comp")
    comp_params = branch_bias.BranchBiasParams(
        seed=int(comp_seed),
        bias_mode=str(branch_comp_params.get("bias_mode", "fixed")),
        fixed_value=float(branch_comp_params.get("fixed_value", 0.0)),
        draw_mean=float(branch_comp_params.get("draw_mean", 0.0)),
        draw_kappa=float(branch_comp_params.get("draw_kappa", 50.0)),
        map_values=map_values,
    )
    comp_state = branch_bias.BranchBiasState(comp_params)

    results: Dict[str, Any] = {
        "seed": seed,
        "branch_comp_params": {
            k: v for k, v in branch_comp_params.items() if k != "map_values"
        },
        "module_paths": dict(module_paths) if module_paths else None,
        "phyllotaxy": phyllo_summary,
        "root_id": topo["root_id"],
        "organ_events": [],
        "branch_events": [],
        "event_records": [],
        "branch_starts": [],
    }

    def _np_rng(*parts: Any) -> np.random.Generator:
        return np.random.default_rng(_stable_seed_int(base_seed, *parts))

    def _py_rng(*parts: Any) -> random.Random:
        return random.Random(_stable_seed_int(base_seed, *parts))

    def _prepare_branch_start(branch_id: str, raw_state: List[self_renewal.CellState], age_offset: int) -> tuple[List[self_renewal.CellState], int, float]:
        seed_index = int(comp_state.branch_seed(branch_id, niche_size=len(raw_state)))
        branch_w = float(comp_state.branch_bias(branch_id))
        seeded_state = self_renewal.assign_branch_comp_seed(raw_state, seed_index)
        results["branch_starts"].append(
            {
                "branch_id": str(branch_id),
                "age_offset": int(age_offset),
                "favored_seed_index": int(seed_index),
                "branch_comp_bias": float(branch_w),
            }
        )
        return seeded_state, seed_index, branch_w

    # ---------------------------------------------------------------------------
    # Iterative DFS traversal (replaces the previous recursive simulate_branch).
    #
    # Why iterative: Python's default recursion limit (~1000 frames) would be hit
    # by deep or bushy topologies (>~300 branches in a chain, or combinatorially
    # deep trees).  An explicit stack is unlimited and has the same O(branches)
    # memory footprint as the recursive version.
    #
    # Traversal order: DFS pre-order, children pushed in reverse event order so
    # that the first BRANCH event target is processed first — identical to the
    # depth-first left-to-right order the recursive version produced.
    #
    # State threading: each branch's (init_state, age_offset) is written into
    # _pending by its parent's BRANCH event handler and consumed when that branch
    # is popped from the stack.  The root entry is seeded directly.
    # ---------------------------------------------------------------------------

    # Maps branch_id -> (init_state, age_offset) for branches not yet processed.
    # Root has no parent, so its entry is seeded here with init_state=None.
    _pending: Dict[str, tuple] = {str(topo["root_id"]): (None, 0)}

    # DFS stack of branch_ids to process.  Root goes on first.
    _stack: List[str] = [str(topo["root_id"])]

    while _stack:
        branch_id = _stack.pop()
        init_state, age_offset = _pending.pop(branch_id)

        spec = topo["branches"][branch_id]
        T = int(spec["T"])
        events = sorted(list(spec.get("events", [])), key=_event_sort_key)

        if init_state is None:
            raw_state = reg.self_renewal.initialize_founders(sr_params, rng=_np_rng("sr_init", branch_id))
        else:
            raw_state = list(init_state)
            if len(raw_state) != int(sr_params.m):
                raise ValueError(f"init_state length must be m={sr_params.m}")

        state, _, branch_w = _prepare_branch_start(str(branch_id), raw_state, int(age_offset))
        t_prev = 0

        # Collect child branch ids in the order they are encountered so we can
        # push them onto the stack in reverse (preserving left-to-right DFS order).
        _children_this_branch: List[str] = []

        for e_idx, e in enumerate(events):
            t_local = int(e["time"])
            if t_local < t_prev or t_local > T:
                raise ValueError(f"Invalid event time {t_local} on branch {branch_id} with T={T}")

            dt = t_local - t_prev
            if dt > 0:
                sr_local = replace(sr_params, branch_comp_bias=float(branch_w))
                state = reg.self_renewal.simulate_segment(
                    sr_local,
                    branch_id=str(branch_id),
                    T=int(dt),
                    age_offset=int(age_offset + t_prev),
                    init_state=list(state),
                    rng=_np_rng("sr", branch_id, "seg", age_offset + t_prev, dt),
                    snapshot_times=None,
                )

            event_age = int(age_offset + t_local)
            etype = str(e["type"]).upper()
            target = str(e["target_id"])
            phyllo_angle_rad = e.get("phyllo_angle_rad")
            phyllo_angle_deg = e.get("phyllo_angle_deg")
            phyllo_site_id = e.get("phyllo_site_id")
            phyllo_site_rank = e.get("phyllo_site_rank")

            pre_py = _py_rng("pre", branch_id, event_age, e_idx)
            pre_np = _np_rng("pre", branch_id, event_age, e_idx)
            ring, end_dev = reg.pre_branching.amplify(
                sr_params,
                state,
                branch_id=str(branch_id),
                event_age=int(event_age),
                boundary_cells_target=int(sam_boundary_cells),
                start_dev=0,
                seed=int(pre_py.randrange(2**31)),
                rng=pre_np,
            )

            if etype == "BRANCH":
                child_id = target
                focal_index = (
                    phyllotaxy.angle_to_focal_index(float(phyllo_angle_rad), len(ring))
                    if phyllo_angle_rad is not None else None
                )
                sel = reg.branching.sample_child_niche(
                    ring,
                    branch_id=str(branch_id),
                    child_id=str(child_id),
                    event_time=(int(event_age), int(end_dev)),
                    precursor_number=int(branch_precursor_number),
                    niche_size_m=int(sr_params.m),
                    rng=_py_rng("branch", branch_id, child_id, event_age, e_idx),
                    focal_index=focal_index,
                )

                be = {
                    "branch_id": str(sel["branch_id"]),
                    "child_id": str(sel["child_id"]),
                    "event_time": [int(sel["event_time"][0]), int(sel["event_time"][1])],
                    "niche_size_m": int(sel["niche_size_m"]),
                    "precursor_number": int(sel["precursor_number"]),
                    "requested_precursor_number": int(sel.get("requested_precursor_number", sel["precursor_number"])),
                    "precursor_number_capped": bool(sel.get("precursor_number_capped", False)),
                    "focal_index": int(sel["focal_index"]),
                    "sampled_indices": list(sel["sampled_indices"]),
                    "precursor_quotas": list(sel.get("precursor_quotas", [])),
                    "precursor_wraparound": bool(sel.get("precursor_wraparound", False)),
                    "precursor_min_index": sel.get("precursor_min_index"),
                    "precursor_max_index": sel.get("precursor_max_index"),
                    "founder_lineage_counts": dict(sel.get("founder_lineage_counts", {})),
                    "founder_sector_count": int(sel.get("founder_sector_count", 0)),
                    "founder_effective_sectors": float(sel.get("founder_effective_sectors", 0.0)),
                    "founder_diversity": float(sel.get("founder_diversity", 0.0)),
                    "founder_dominant_fraction": float(sel.get("founder_dominant_fraction", 0.0)),
                }
                if phyllo_site_id is not None:
                    be["phyllo_site_id"] = str(phyllo_site_id)
                if phyllo_site_rank is not None:
                    be["phyllo_site_rank"] = int(phyllo_site_rank)
                if phyllo_angle_deg is not None:
                    be["phyllo_angle_deg"] = float(phyllo_angle_deg)
                results["branch_events"].append(be)
                results["event_records"].append(
                    {
                        "module": "branching",
                        "event_type": "BRANCH",
                        "branch_id": str(branch_id),
                        "target_id": str(child_id),
                        "event_time": be["event_time"],
                        "outputs": be,
                    }
                )

                # Register child state for when it is popped from the stack.
                _pending[str(child_id)] = (list(sel["expanded_state"]), int(event_age))
                _children_this_branch.append(str(child_id))

            elif etype == "ORGAN":
                org_py = _py_rng("organ", branch_id, target, event_age, e_idx)
                org_np = _np_rng("organ", branch_id, target, event_age, e_idx)
                focal_index = (
                    phyllotaxy.angle_to_focal_index(float(phyllo_angle_rad), len(ring))
                    if phyllo_angle_rad is not None else None
                )
                oe = dict(reg.organ.form_organ(
                    sr_params,
                    ring,
                    organ_id=str(target),
                    event_age=int(event_age),
                    start_dev=int(end_dev),
                    precursor_number=int(organ_precursor_number),
                    organ_total_cells=int(organ_total_cells),
                    sequenced_cells=sequenced_cells,
                    rng=org_py,
                    np_rng=org_np,
                    focal_index=focal_index,
                ))
                if phyllo_site_id is not None:
                    oe["phyllo_site_id"] = str(phyllo_site_id)
                if phyllo_site_rank is not None:
                    oe["phyllo_site_rank"] = int(phyllo_site_rank)
                if phyllo_angle_deg is not None:
                    oe["phyllo_angle_deg"] = float(phyllo_angle_deg)
                if "event_time" in oe and isinstance(oe["event_time"], tuple):
                    oe["event_time"] = [int(oe["event_time"][0]), int(oe["event_time"][1])]
                if "allele_counts_by_mutation" in oe:
                    oe["allele_counts_by_mutation"] = {int(k): int(v) for k, v in oe["allele_counts_by_mutation"].items()}

                results["organ_events"].append(oe)
                results["event_records"].append(
                    {
                        "module": "organ",
                        "event_type": "ORGAN",
                        "branch_id": str(branch_id),
                        "target_id": str(target),
                        "event_time": oe["event_time"],
                        "outputs": oe,
                    }
                )
            else:
                raise ValueError(f"Unknown event type {etype} on branch {branch_id}")

            t_prev = t_local

        # Tail-segment self-renewal (T - t_last) is intentionally omitted.
        # Per the simSOMA model (Section 2.5 of the paper), self-renewal runs
        # only *between* events.  After the last event on a branch there are no
        # further module calls whose inputs depend on the niche state, so
        # simulating the tail would accumulate mutations that are never observed.

        # Push children in reverse order so the first child encountered is
        # processed first (preserving DFS left-to-right order).
        for child_id in reversed(_children_this_branch):
            _stack.append(child_id)

    # branch_competition summary.
    # Note: comp_state.realize_all() is intentionally not called here.
    # Every branch already called branch_bias() and branch_seed() lazily via
    # _prepare_branch_start() during the traversal above, so the cache is fully
    # populated at this point.  realize_all() would be a no-op.
    branch_ids = list(topo.get("branches", {}).keys())
    results["branch_competition"] = {
        "bias_mode": str(comp_params.bias_mode),
        "bias_by_branch": {str(bid): float(comp_state.branch_bias(str(bid))) for bid in branch_ids},
        "seed_index_by_branch": {str(bid): int(comp_state.branch_seed(str(bid), niche_size=int(sr_params.m))) for bid in branch_ids},
    }

    agg_union = set()
    for oe in results.get("organ_events", []):
        acm = oe.get("allele_counts_by_mutation")
        if isinstance(acm, dict):
            agg_union |= set(map(int, acm.keys()))
    results["aggregate_total_variants"] = int(len(agg_union))

    return results


def _strip_internal_keys(result: Dict[str, Any]) -> None:
    for oe in result.get("organ_events", []):
        for k in list(oe.keys()):
            if str(k).startswith("__"):
                del oe[k]
    for er in result.get("event_records", []):
        outs = er.get("outputs")
        if isinstance(outs, dict):
            for k in list(outs.keys()):
                if str(k).startswith("__"):
                    del outs[k]


def _parse_float_list(s: str) -> List[float]:
    parts = [p.strip() for p in str(s).split(",") if p.strip() != ""]
    return [float(x) for x in parts]


def _parse_int_list(s: str) -> List[int]:
    parts = [p.strip() for p in str(s).split(",") if p.strip() != ""]
    return [int(x) for x in parts]


def _cli() -> None:
    p = argparse.ArgumentParser(description="Run the plant development pipeline on a topology spec (JSON or CSV).")

    p.add_argument("--topology_json", type=Path, required=False, help="Topology JSON (user-friendly or internal schema).")
    p.add_argument("--branches_csv", type=Path, required=False, help="Branches CSV file (manual entry friendly).")
    p.add_argument("--events_csv", type=Path, required=False, help="Events CSV file (manual entry friendly).")

    p.add_argument(
        "--export_topology_json",
        type=str,
        default="auto",
        choices=["auto", "yes", "no"],
        help="If 'yes', write the validated internal topology JSON next to --out_json (or to --export_topology_path). If 'auto', exports only when CSV input is used.",
    )
    p.add_argument(
        "--export_topology_path",
        type=Path,
        default=None,
        help="Optional output path for exported internal topology JSON. If omitted, uses <out_json stem>.topology_<idx>.json",
    )
    p.add_argument(
        "--stop_after_topology",
        type=str,
        default="no",
        choices=["yes", "no"],
        help="If 'yes', run only topology import/conversion/validation/mapping (and optional export), then write a topology report to --out_json and exit.",
    )
    p.add_argument("--out_json", type=Path, required=True)

    p.add_argument(
        "--topology_unit",
        type=str,
        default="auto",
        choices=["auto", "steps", "years", "meters"],
        help="Unit of T and event times in topology JSON. If not 'steps', values will be mapped to SR steps using --kappa_sr. If 'auto', unit is taken from user-friendly topology JSON if available, otherwise 'steps'.",
    )
    p.add_argument(
        "--topology_mapping_mode",
        type=str,
        default="deterministic",
        choices=["deterministic"],
        help="How to map observed branch lengths to SR steps.",
    )
    p.add_argument(
        "--topology_mapping_seed",
        type=int,
        default=None,
        help="Seed for the topology mapping stochasticity (defaults to --seed or per-run derived seed).",
    )

    p.add_argument("--m", type=int, required=True)
    p.add_argument("--rho", type=float, required=True)
    p.add_argument(
        "--mu_unit",
        "--mu_year",
        dest="mu_unit",
        type=str,
        required=True,
        help="Somatic mutation rate per lineage per topology unit (per genome). --mu_year is a deprecated alias. Can be a comma-separated list for sweeps.",
    )
    p.add_argument(
        "--kappa_sr",
        type=str,
        required=True,
        help="Self-renewal divisions per topology unit. This value is also used for topology-unit to SR-step conversion.",
    )
    p.add_argument(
        "--mu_div",
        type=str,
        default=None,
        help="Optional mutation rate per division. If provided together with mu_unit and kappa_sr, the values must be consistent (mu_div ~= mu_unit/kappa_sr). Can be a comma-separated list for sweeps.",
    )
    p.add_argument("--victim_locality", type=float, default=0.0, help="Victim locality in [0,1]: 0 = uniform across the niche, 1 = nearest-neighbor-only.")

    p.add_argument(
        "--phyllotaxy_mode",
        type=str,
        default="off",
        choices=["off", "random", "spiral", "distichous", "tristichous", "decussate"],
        help="Optional site-based phyllotaxy overlay applied to observed events on each axis.",
    )
    p.add_argument(
        "--phyllotaxy_divergence_deg",
        type=float,
        default=None,
        help="Divergence angle in degrees for --phyllotaxy_mode=spiral.",
    )
    p.add_argument(
        "--phyllotaxy_tie_tol",
        type=float,
        default=1e-9,
        help="Absolute tolerance for grouping same-time events into one shared phyllotactic site.",
    )

    p.add_argument("--bias_mode", type=str, default="fixed", choices=["fixed", "draw", "map"], help="How branch-local displacement biases are assigned.")
    p.add_argument("--branch_bias_value", type=float, default=0.0, help="Fixed branch-local displacement bias w in [0,1] when bias_mode=fixed.")
    p.add_argument("--branch_bias_mean", type=float, default=0.0, help="Mean of the branch-local displacement bias distribution when bias_mode=draw.")
    p.add_argument("--branch_bias_kappa", type=float, default=50.0, help="Concentration of the branch-local displacement bias Beta distribution when bias_mode=draw.")
    p.add_argument("--branch_bias_map_json", type=Path, default=None, help="JSON file mapping branch_id -> bias value in [0,1] for bias_mode=map. Optional key 'default' is allowed.")

    p.add_argument("--sam_boundary_cells", type=int, default=64)
    p.add_argument("--branch_precursor_number", type=int, default=1)
    p.add_argument("--organ_precursor_number", type=int, default=1)
    p.add_argument("--organ_total_cells", type=int, default=64)
    p.add_argument("--sequenced_cells", type=int, default=None)
    p.add_argument("--seed", type=int, default=None)

    p.add_argument("--plugin_organ", type=str, default=None, help="Optional organ module plug-in path (import path).")
    p.add_argument("--plugin_branching", type=str, default=None, help="Optional branching module plug-in path (import path).")
    p.add_argument("--plugin_pre_branching", type=str, default=None, help="Optional pre-branching module plug-in path (import path).")
    p.add_argument("--plugin_self_renewal", type=str, default=None, help="Optional self-renewal module plug-in path (import path).")
    p.add_argument(
        "--summaries",
        type=str,
        default="no",
        choices=["yes", "no"],
        help="If 'yes', compute summary statistics via summaries.py for each simulation result.",
    )
    p.add_argument(
        "--vaf_nbins",
        type=int,
        default=20,
        help="Number of bins for binned VAF spectra in summaries.",
    )
    p.add_argument(
        "--summary_private_shared",
        type=str,
        default="yes",
        choices=["yes", "no"],
        help="If 'yes', compute private vs shared mutation summaries across organs (requires internal per-mutation counts during runtime).",
    )
    p.add_argument(
        "--n_sim",
        type=int,
        default=1,
        help="Number of replicate simulations per parameter combination (uses distinct seeds).",
    )

    args = p.parse_args()

    using_csv = args.branches_csv is not None or args.events_csv is not None
    using_json = args.topology_json is not None
    if using_csv and using_json:
        raise ValueError("Provide either --topology_json OR (--branches_csv and --events_csv), not both.")
    if using_csv:
        if args.branches_csv is None or args.events_csv is None:
            raise ValueError("CSV topology requires both --branches_csv and --events_csv.")
    if not using_csv and not using_json:
        raise ValueError("You must provide a topology via --topology_json or via --branches_csv/--events_csv.")

    mu_unit_list = _parse_float_list(args.mu_unit)
    kappa_list = _parse_float_list(args.kappa_sr)
    mu_div_list = _parse_float_list(args.mu_div) if args.mu_div is not None else []
    grid = list(itertools.product(mu_unit_list, kappa_list, (mu_div_list if mu_div_list else [float("nan")]), [None]))

    if int(args.n_sim) < 1:
        raise ValueError("--n_sim must be >= 1")

    runs: List[Dict[str, Any]] = []
    base_seed = args.seed

    for idx, (mu_unit, kappa_sr, mu_div, _unused) in enumerate(grid):
        if math.isnan(mu_div):
            if kappa_sr <= 0:
                raise ValueError("kappa_sr must be > 0")
            mu_div_eff = float(mu_unit) / float(kappa_sr)
        else:
            mu_div_eff = float(mu_div)
            if mu_div_eff <= 0:
                raise ValueError("mu_div must be > 0")
            denom = max(1e-12, abs(mu_div_eff), abs(float(mu_unit) / float(kappa_sr)))
            if abs(mu_div_eff - float(mu_unit) / float(kappa_sr)) / denom > 1e-6:
                raise ValueError("Inconsistent rates: expected mu_div ~= mu_unit/kappa_sr")

        topo_unit = str(args.topology_unit).lower()
        if topo_unit == "auto":
            if using_json:
                try:
                    raw = json.loads(args.topology_json.read_text(encoding="utf-8"))
                    topo_unit = str(raw.get("unit", "steps")).lower() if isinstance(raw, dict) else "steps"
                except Exception:
                    topo_unit = "steps"
            else:
                topo_unit = "steps"
        if topo_unit not in {"steps", "years", "meters"}:
            raise ValueError("--topology_unit must be one of: auto, steps, years, meters")

        mapping = None
        if topo_unit != "steps":
            mapping = {"unit": topo_unit, "rate": float(kappa_sr), "mode": str(args.topology_mapping_mode)}

        if args.topology_mapping_seed is not None:
            map_seed_combo = int(args.topology_mapping_seed)
        elif base_seed is not None:
            map_seed_combo = _stable_seed_int(int(base_seed), "topology_map", idx)
        else:
            map_seed_combo = int.from_bytes(os.urandom(8), "big") & ((1 << 63) - 1)

        if using_json:
            topo_base = topology_io.load_topology_auto(
                args.topology_json,
                mapping=mapping,
                seed=int(map_seed_combo),
                phyllotaxy_config={
                    "mode": str(args.phyllotaxy_mode),
                    "divergence_angle_deg": args.phyllotaxy_divergence_deg,
                    "tie_tol": float(args.phyllotaxy_tie_tol),
                },
                phyllotaxy_seed=int(map_seed_combo),
            )
        else:
            topo_base = topology_io.read_topology_csv(
                args.branches_csv,
                args.events_csv,
                mapping=mapping,
                seed=int(map_seed_combo),
                phyllotaxy_config={
                    "mode": str(args.phyllotaxy_mode),
                    "divergence_angle_deg": args.phyllotaxy_divergence_deg,
                    "tie_tol": float(args.phyllotaxy_tie_tol),
                },
                phyllotaxy_seed=int(map_seed_combo),
            )

        do_export = str(args.export_topology_json).lower()
        if do_export == "auto":
            do_export = "yes" if using_csv else "no"
        if do_export == "yes":
            if args.export_topology_path is not None:
                export_path = Path(args.export_topology_path)
            else:
                export_path = args.out_json.with_suffix("").with_name(f"{args.out_json.stem}.topology_{idx}.json")
            topology_io.save_topology_json(topo_base, export_path)

        if str(args.stop_after_topology).lower() == "yes":
            runs.append(
                {
                    "idx": idx,
                    "parameters": {
                        "mu_unit": float(mu_unit),
                        "mu_year": float(mu_unit),  # deprecated alias
                        "kappa_sr": float(kappa_sr),
                        "mu_div": float(mu_div_eff),
                        "topology_unit": topo_unit,
                        "topology_mapping_mode": None if mapping is None else str(mapping["mode"]),
                        "topology_rate": None if mapping is None else float(mapping["rate"]),
                        "topology_rate_source": None if mapping is None else "kappa_sr",
                        "topology_mapping_seed": int(map_seed_combo),
                    },
                    "topology_summary": summarize_topology(topo_base),
                    "topology": topo_base,
                }
            )
            continue

        sr_params = self_renewal.SelfRenewalParams(
            m=int(args.m),
            rho=float(args.rho),
            mu_unit=float(mu_unit),
            kappa_sr=float(kappa_sr),
            mu_div=float(mu_div_eff),
            victim_locality=float(args.victim_locality),
            branch_comp_bias=0.0,
        )

        branch_comp_params = {
            "bias_mode": str(args.bias_mode),
            "fixed_value": float(args.branch_bias_value),
            "draw_mean": float(args.branch_bias_mean),
            "draw_kappa": float(args.branch_bias_kappa),
            "map_json": str(args.branch_bias_map_json) if args.branch_bias_map_json is not None else None,
        }

        sims: List[Dict[str, Any]] = []
        for rep in range(int(args.n_sim)):
            run_seed = (int(base_seed) + idx * int(args.n_sim) + rep) if base_seed is not None else None
            topo = copy.deepcopy(topo_base)

            module_paths = {}
            if args.plugin_organ:
                module_paths["organ"] = str(args.plugin_organ)
            if args.plugin_branching:
                module_paths["branching"] = str(args.plugin_branching)
            if args.plugin_pre_branching:
                module_paths["pre_branching"] = str(args.plugin_pre_branching)
            if args.plugin_self_renewal:
                module_paths["self_renewal"] = str(args.plugin_self_renewal)
            if len(module_paths) == 0:
                module_paths = None

            res = run_pipeline(
                topo,
                sr_params,
                sam_boundary_cells=int(args.sam_boundary_cells),
                branch_precursor_number=int(args.branch_precursor_number),
                organ_precursor_number=int(args.organ_precursor_number),
                organ_total_cells=int(args.organ_total_cells),
                sequenced_cells=args.sequenced_cells,
                seed=run_seed,
                branch_comp_params=branch_comp_params,
                module_paths=module_paths,
            )

            sim_entry: Dict[str, Any] = {"rep": rep, "seed": run_seed, "topology_mapping_seed": int(map_seed_combo), "result": res}

            if str(args.summaries).lower() == "yes":
                summ = summaries.compute_summary(
                    res,
                    nbins=int(args.vaf_nbins),
                    private_shared=(str(args.summary_private_shared).lower() == "yes"),
                )
                sim_entry["summary"] = summ

            _strip_internal_keys(res)
            sims.append(sim_entry)

        runs.append(
            {
                "params": {
                    "m": int(args.m),
                    "rho": float(args.rho),
                    "mu_unit": float(mu_unit),
                    "mu_year": float(mu_unit),  # deprecated alias
                    "kappa_sr": float(kappa_sr),
                    "mu_div": float(mu_div_eff),
                    "victim_locality": float(args.victim_locality),
                    "bias_mode": str(args.bias_mode),
                    "branch_bias_value": float(args.branch_bias_value),
                    "branch_bias_mean": float(args.branch_bias_mean),
                    "branch_bias_kappa": float(args.branch_bias_kappa),
                    "branch_bias_map_json": None if args.branch_bias_map_json is None else str(args.branch_bias_map_json),
                    "sam_boundary_cells": int(args.sam_boundary_cells),
                    "branch_precursor_number": int(args.branch_precursor_number),
                    "organ_precursor_number": int(args.organ_precursor_number),
                    "organ_total_cells": int(args.organ_total_cells),
                    "sequenced_cells": args.sequenced_cells,
                    "seed": base_seed,
                    "n_sim": int(args.n_sim),
                    "topology_unit": str(args.topology_unit),
                    "topology_mapping_mode": str(args.topology_mapping_mode),
                    "topology_rate": None if mapping is None else float(mapping["rate"]),
                    "topology_rate_source": None if mapping is None else "kappa_sr",
                    "topology_mapping_seed": None if args.topology_mapping_seed is None else int(args.topology_mapping_seed),
                    "summaries": str(args.summaries).lower(),
                    "vaf_nbins": int(args.vaf_nbins),
                    "summary_private_shared": str(args.summary_private_shared).lower(),
                },
                "sims": sims,
                "summary_aggregate": summaries.aggregate_replicate_summaries([s.get("summary") for s in sims if "summary" in s]) if (str(args.summaries).lower() == "yes" and len(sims) > 1) else None,
            }
        )

    if str(args.stop_after_topology).lower() == "yes":
        out = {"runs": runs}
    elif str(args.summaries).lower() == "yes" or len(runs) != 1 or int(args.n_sim) != 1:
        out = {"runs": runs}
    else:
        out = runs[0]["sims"][0]["result"]

    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")


if __name__ == "__main__":
    _cli()
