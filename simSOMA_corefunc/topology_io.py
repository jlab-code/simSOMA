
"""
topology_io.py

Lightweight IO utilities to convert the pipeline's CSV topology format to JSON
and load/save topology specs.

CSV format (as used previously)
------------------------------
branches.csv columns:
  - branch_id
  - parent_id   (empty for root)
  - T           (either int self-renewal steps, or observed age/length depending on mapping)

events.csv columns:
  - branch_id
  - time        (either int self-renewal steps, or observed time/length depending on mapping)
  - type        ("BRANCH" or "ORGAN")
  - target_id   (child_branch_id for BRANCH, organ_id for ORGAN)

JSON format
-----------
{
  "root_id": "B0",
  "branches": {
    "B0": { "parent_id": null, "T": 400, "events": [ {"time":80,"type":"BRANCH","target_id":"B1"}, ... ] },
    ...
  }
}

Optional mapping
----------------
If your topology is given in *observed units* (e.g. years, meters), you can map it to
integer self-renewal steps via a user-specified rate.

We store times in the JSON under the same keys (T, time) and overwrite them with integer
self-renewal steps after mapping. The original observed branch length is retained as T_obs.

Mapping methods:
- deterministic: T_steps = round(rate * T_obs)
- poisson:       T_steps ~ Poisson(rate * T_obs)

Event times are mapped by preserving within-branch position:
time_frac = time_obs / T_obs, then time_steps = round(time_frac * T_steps).
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
import math
import random

import phyllotaxy


Topology = Dict[str, Any]


def annotate_topology_phyllotaxy(
    topo: Topology,
    *,
    phyllotaxy_config: Optional[Dict[str, Any]] = None,
    seed: Optional[int] = None,
) -> Tuple[Topology, Dict[str, Any]]:
    """Attach site-based phyllotaxy annotations to a topology.

    This is kept in topology_io so topology preparation stays modular: the
    topology layer owns event/site annotation, while the wrapper only consumes
    already-annotated event metadata.
    """
    topo_out, summary = phyllotaxy.annotate_topology_sites(
        topo,
        config=phyllotaxy_config,
        seed=seed,
    )
    topo_out["phyllotaxy"] = dict(summary)
    return topo_out, summary

def validate_topology(topo: Topology, *, unit: str = "steps") -> None:
    """Validate internal topology JSON structure.

    Raises ValueError with concrete messages on invalid input.
    """
    if not isinstance(topo, dict):
        raise ValueError("Topology must be a dict")
    if "root_id" not in topo or "branches" not in topo:
        raise ValueError("Invalid topology: must contain keys root_id and branches")
    branches = topo["branches"]
    if not isinstance(branches, dict) or not branches:
        raise ValueError("Invalid topology: branches must be a non-empty dict")

    # root existence
    root_id = str(topo["root_id"])
    if root_id not in branches:
        raise ValueError(f"root_id {root_id} not found in branches")

    # parent consistency and unique root
    roots = []
    for bid, b in branches.items():
        pid = b.get("parent_id")
        if pid is None:
            roots.append(str(bid))
        else:
            if str(pid) not in branches:
                raise ValueError(f"Branch {bid} references missing parent_id {pid}")
    if len(roots) != 1:
        raise ValueError(f"Expected exactly one root (parent_id null). Found roots={roots}")
    if roots[0] != root_id:
        raise ValueError(f"root_id={root_id} but detected root branch is {roots[0]}")

    # cycle check via DFS
    children: Dict[str, List[str]] = {}
    for bid, b in branches.items():
        pid = b.get("parent_id")
        if pid is not None:
            children.setdefault(str(pid), []).append(str(bid))

    seen: set[str] = set()
    stack: set[str] = set()

    def dfs(u: str) -> None:
        if u in stack:
            raise ValueError(f"Cycle detected in topology at branch {u}")
        if u in seen:
            return
        stack.add(u)
        for v in children.get(u, []):
            dfs(v)
        stack.remove(u)
        seen.add(u)

    dfs(root_id)
    if len(seen) != len(branches):
        missing = sorted(set(map(str, branches.keys())) - seen)
        raise ValueError(f"Topology is disconnected; unreachable branches: {missing}")

    # event checks
    for bid, b in branches.items():
        if "T" not in b:
            raise ValueError(f"Branch {bid} missing T")
        T = b["T"]
        if unit == "steps":
            try:
                T_int = int(T)
            except Exception:
                raise ValueError(f"Branch {bid} T must be int for unit=steps")
            if T_int < 0:
                raise ValueError(f"Branch {bid} T must be >= 0")
        else:
            try:
                T_float = float(T)
            except Exception:
                raise ValueError(f"Branch {bid} T must be numeric for unit={unit}")
            if T_float < 0:
                raise ValueError(f"Branch {bid} T must be >= 0")

        events = list(b.get("events", []))
        last_t = -float("inf")
        for e in events:
            if "time" not in e or "type" not in e or "target_id" not in e:
                raise ValueError(f"Event on branch {bid} missing required keys (time,type,target_id)")
            et = str(e["type"]).upper()
            if et not in {"BRANCH", "ORGAN"}:
                raise ValueError(f"Unknown event type {et} on branch {bid}")
            t = e["time"]
            try:
                t_num = int(t) if unit == "steps" else float(t)
            except Exception:
                raise ValueError(f"Invalid event time {t} on branch {bid} for unit={unit}")
            if t_num < last_t:
                raise ValueError(f"Events on branch {bid} are not non-decreasing by time")
            last_t = t_num
            # bounds check
            if unit == "steps":
                if t_num < 0 or t_num > int(b["T"]):
                    raise ValueError(f"Event time {t_num} out of bounds on branch {bid} (T={b['T']})")
            else:
                if float(t_num) < 0.0 or float(t_num) > float(b["T"]):
                    raise ValueError(f"Event time {t_num} out of bounds on branch {bid} (T={b['T']})")
            # target checks
            target = str(e["target_id"])
            if et == "BRANCH" and target not in branches:
                raise ValueError(f"BRANCH event on {bid} targets missing child branch {target}")

    # organ target uniqueness (helpful)
    organ_targets = []
    for bid, b in branches.items():
        for e in b.get("events", []):
            if str(e["type"]).upper() == "ORGAN":
                organ_targets.append(str(e["target_id"]))
    if len(organ_targets) != len(set(organ_targets)):
        # don't hard fail, but warn via ValueError? Keep as error for now to avoid silent overwrites downstream.
        dup = sorted({x for x in organ_targets if organ_targets.count(x) > 1})
        raise ValueError(f"Duplicate organ target_id(s) found: {dup}")


def read_topology_user_json(path: Path) -> Tuple[Topology, str]:
    """Read a user-friendly topology JSON and convert to internal format.

    User-friendly schema (recommended)
    --------------------------------
    {
      "unit": "years" | "meters" | "steps",
      "branches": [{"id":"B0","parent":null,"length":12.3}, ...],
      "events": [
        {"branch":"B0","pos":0.35,"type":"BRANCH","target":"B1"},
        {"branch":"B1","time":3.2,"type":"ORGAN","target":"O1"}
      ]
    }

    Semantics
    ---------
    - "length" is the observed branch length in the given unit.
    - Each event provides either:
        * "pos" in [0,1] (fraction along the branch), or
        * "time" in the same unit as "length"
    - Output internal topology uses keys root_id, branches{...}, where each branch has:
        * T (observed), and events with numeric time (observed) before any SR-step mapping.

    Returns
    -------
    topo_internal, unit
    """
    raw = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise ValueError("Topology JSON must be an object")

    unit = str(raw.get("unit", "steps")).lower()
    if unit not in {"steps", "years", "meters"}:
        raise ValueError("user topology 'unit' must be one of: steps, years, meters")

    branches_raw = raw.get("branches")
    events_raw = raw.get("events")
    if not isinstance(branches_raw, list) or not branches_raw:
        raise ValueError("user topology must contain non-empty 'branches' list")
    if not isinstance(events_raw, list):
        raise ValueError("user topology must contain 'events' list")

    branches: Dict[str, Dict[str, Any]] = {}
    for b in branches_raw:
        if not isinstance(b, dict):
            raise ValueError("Each branch entry must be an object")
        bid = str(b.get("id", "")).strip()
        if not bid:
            raise ValueError("Branch entry missing 'id'")
        if bid in branches:
            raise ValueError(f"Duplicate branch id {bid} in user topology")
        pid = b.get("parent", None)
        pid = (str(pid).strip() if pid is not None and str(pid).strip() != "" else None)
        if "length" not in b and "T" not in b:
            raise ValueError(f"Branch {bid} missing 'length' (or 'T')")
        length = b.get("length", b.get("T"))
        T_val = float(length) if unit != "steps" else int(length)
        branches[bid] = {"parent_id": pid, "T": T_val, "events": []}

    # determine root
    roots = [bid for bid, b in branches.items() if b["parent_id"] is None]
    if len(roots) != 1:
        raise ValueError(f"Expected exactly one root branch (parent null). Found: {roots}")
    root_id = roots[0]

    # attach events
    for e in events_raw:
        if not isinstance(e, dict):
            raise ValueError("Each event entry must be an object")
        bid = str(e.get("branch", "")).strip()
        if bid not in branches:
            raise ValueError(f"Event references unknown branch '{bid}'")
        etype = str(e.get("type", "")).strip().upper()
        if etype not in {"BRANCH", "ORGAN"}:
            raise ValueError(f"Unknown event type '{etype}' on branch {bid}")
        target = str(e.get("target", e.get("target_id", ""))).strip()
        if not target:
            raise ValueError(f"Event on branch {bid} missing 'target'")

        T_obs = float(branches[bid]["T"]) if unit != "steps" else float(int(branches[bid]["T"]))
        if "pos" in e:
            pos = float(e["pos"])
            if pos < 0.0 or pos > 1.0:
                raise ValueError(f"Event pos must be in [0,1] on branch {bid}")
            t_obs = pos * T_obs
        elif "time" in e:
            t_obs = float(e["time"]) if unit != "steps" else int(e["time"])
        else:
            raise ValueError(f"Event on branch {bid} must provide 'pos' or 'time'")

        branches[bid]["events"].append({"time": t_obs, "type": etype, "target_id": target})

    # sort events by observed time
    for bid in branches:
        branches[bid]["events"] = sorted(branches[bid]["events"], key=lambda x: float(x["time"]))

    topo_internal: Topology = {"root_id": root_id, "branches": branches}
    # validate in observed unit before mapping
    validate_topology(topo_internal, unit=unit)
    return topo_internal, unit


def load_topology_auto(
    path: Path,
    *,
    mapping: Optional[Dict[str, Any]] = None,
    seed: Optional[int] = None,
    phyllotaxy_config: Optional[Dict[str, Any]] = None,
    phyllotaxy_seed: Optional[int] = None,
) -> Topology:
    """Load either internal topology JSON or user-friendly topology JSON.

    Auto-detection:
    - Internal format: has keys root_id and branches (dict).
    - User format: has keys branches (list) and events (list), optionally unit.

    mapping:
      If provided, overrides/defines mapping from observed units to SR steps.
      This should be used when the input is in years/meters and you want SR steps.
    """
    raw = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(raw, dict) and "root_id" in raw and isinstance(raw.get("branches"), dict):
        topo = raw
        # validate in given unit before mapping
        validate_topology(topo, unit=(mapping.get("unit","steps") if mapping else "steps"))
        if mapping is not None:
            topo = map_topology_to_sr_steps(
                topo,
                unit=str(mapping.get("unit", "steps")),
                rate=float(mapping.get("rate")),
                mode=str(mapping.get("mode", "poisson")),
                seed=seed,
                keep_observed=True,
            )
        validate_topology(topo, unit="steps" if (mapping and mapping.get("unit","steps")!="steps") else (mapping.get("unit","steps") if mapping else "steps"))
        topo, _ = annotate_topology_phyllotaxy(topo, phyllotaxy_config=phyllotaxy_config, seed=phyllotaxy_seed)
        return topo

    if isinstance(raw, dict) and isinstance(raw.get("branches"), list):
        # user-friendly
        # Re-read using helper for full validation
        topo_obs, unit = read_topology_user_json(path)
        # Decide mapping: if mapping provided, use it; else if unit != steps, do not map unless mapping provided.
        if mapping is not None and str(mapping.get("unit","steps")).lower() != "steps":
            topo = map_topology_to_sr_steps(
                topo_obs,
                unit=str(mapping.get("unit", unit)),
                rate=float(mapping.get("rate")),
                mode=str(mapping.get("mode","poisson")),
                seed=seed,
                keep_observed=True,
            )
            validate_topology(topo, unit="steps")
            topo, _ = annotate_topology_phyllotaxy(topo, phyllotaxy_config=phyllotaxy_config, seed=phyllotaxy_seed)
            return topo
        # no mapping requested
        validate_topology(topo_obs, unit=unit)
        topo_obs, _ = annotate_topology_phyllotaxy(topo_obs, phyllotaxy_config=phyllotaxy_config, seed=phyllotaxy_seed)
        return topo_obs

    raise ValueError("Unrecognized topology JSON format (expected internal or user-friendly schema)")



def _clamp_int(x: int, lo: int, hi: int) -> int:
    return max(lo, min(hi, int(x)))


def map_topology_to_sr_steps(
    topo: Topology,
    *,
    unit: str,
    rate: float,
    mode: str = "poisson",
    seed: Optional[int] = None,
    keep_observed: bool = True,
) -> Topology:
    """Map observed branch age/length to self-renewal steps.

    Parameters
    ----------
    unit:
        "steps" (no mapping), "years", or "meters".
    rate:
        Turnover rate: events per year (if unit="years") or events per meter (if unit="meters").
    mode:
        "poisson" (default) or "deterministic".

    Semantics
    ---------
    - Input topology must have per-branch numeric `T` and per-event numeric `time` in the *observed* unit.
    - Output topology has integer `T` and integer event `time` measured in self-renewal steps.
    - Event times are mapped by *fractional position* along the branch to preserve ordering under
      stochastic branch-length sampling.
    """
    unit = str(unit).lower()
    mode = str(mode).lower()

    if unit == "steps":
        return topo
    if unit not in {"years", "meters"}:
        raise ValueError("unit must be one of: steps, years, meters")
    if rate <= 0:
        raise ValueError("rate must be > 0")
    if mode not in {"poisson", "deterministic"}:
        raise ValueError("mode must be one of: poisson, deterministic")

    rng = random.Random(seed)

    out: Topology = {"root_id": topo["root_id"], "branches": {}}
    out["mapping"] = {"unit": unit, "rate": float(rate), "mode": mode, "seed": seed}

    for bid, b in topo["branches"].items():
        T_obs = float(b["T"])
        if T_obs < 0:
            raise ValueError(f"Branch {bid} has negative T")

        mean_steps = float(rate) * T_obs
        if mode == "deterministic":
            T_steps = int(round(mean_steps))
        else:
            # Poisson via Knuth's algorithm (numpy-free).
            # This is fine for moderate means typical of branch mapping.
            L = math.exp(-mean_steps)
            k = 0
            p = 1.0
            while p > L:
                k += 1
                p *= rng.random()
            T_steps = k - 1

        if T_steps < 0:
            T_steps = 0

        b_out: Dict[str, Any] = {"parent_id": b.get("parent_id"), "T": int(T_steps), "events": []}
        if keep_observed:
            b_out["T_obs"] = float(T_obs)

        # map events by fractional position along branch
        events = list(b.get("events", []))
        mapped_events: List[Dict[str, Any]] = []
        for e in events:
            t_obs = float(e["time"])
            frac = 0.0 if T_obs == 0 else t_obs / T_obs
            frac = max(0.0, min(1.0, frac))
            t_step = int(round(frac * T_steps))
            t_step = _clamp_int(t_step, 0, T_steps)
            ee: Dict[str, Any] = {"time": int(t_step), "type": e["type"], "target_id": e["target_id"]}
            if keep_observed:
                ee["time_obs"] = float(t_obs)
            mapped_events.append(ee)

        # sort + drop None keys
        mapped_events = sorted(mapped_events, key=lambda ee: int(ee["time"]))
        b_out["events"].extend(mapped_events)

        out["branches"][str(bid)] = b_out

    return out


def read_topology_csv(
    branches_csv: Path,
    events_csv: Path,
    *,
    mapping: Optional[Dict[str, Any]] = None,
    seed: Optional[int] = None,
    phyllotaxy_config: Optional[Dict[str, Any]] = None,
    phyllotaxy_seed: Optional[int] = None,
) -> Topology:
    branches: Dict[str, Dict[str, Any]] = {}
    root_id: Optional[str] = None

    with branches_csv.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        required = {"branch_id", "parent_id", "T"}
        if set(reader.fieldnames or []) < required:
            raise ValueError(f"branches_csv missing columns. Need at least {required}.")
        for row in reader:
            bid = row["branch_id"].strip()
            pid = row["parent_id"].strip()
            T_raw = row["T"]
            T_val = float(T_raw) if mapping is not None else int(T_raw)
            branches[bid] = {"parent_id": (pid if pid else None), "T": T_val, "events": []}

    # determine root
    roots = [bid for bid, b in branches.items() if b["parent_id"] is None]
    if len(roots) != 1:
        raise ValueError(f"Expected exactly one root branch (parent_id empty). Found: {roots}")
    root_id = roots[0]

    with events_csv.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        required = {"branch_id", "time", "type", "target_id"}
        if set(reader.fieldnames or []) < required:
            raise ValueError(f"events_csv missing columns. Need at least {required}.")
        for row in reader:
            bid = row["branch_id"].strip()
            if bid not in branches:
                raise ValueError(f"Event references unknown branch_id {bid}.")
            t_raw = row["time"]
            t = float(t_raw) if mapping is not None else int(t_raw)
            etype = row["type"].strip().upper()
            target = row["target_id"].strip()
            if etype not in {"BRANCH", "ORGAN"}:
                raise ValueError(f"Unknown event type {etype}. Must be BRANCH or ORGAN.")
            branches[bid]["events"].append({"time": t, "type": etype, "target_id": target})

    # sort events by time per branch (float-safe)
    for bid in branches:
        branches[bid]["events"] = sorted(branches[bid]["events"], key=lambda e: float(e["time"]))

    topo: Topology = {"root_id": root_id, "branches": branches}

    # Validate pre-mapping in the appropriate unit
    unit0 = str(mapping.get("unit", "steps")).lower() if mapping is not None else "steps"
    validate_topology(topo, unit=unit0)

    if mapping is not None and unit0 != "steps":
        topo = map_topology_to_sr_steps(
            topo,
            unit=unit0,
            rate=float(mapping.get("rate")),
            mode=str(mapping.get("mode", "poisson")),
            seed=seed,
            keep_observed=True,
        )
        validate_topology(topo, unit="steps")
    else:
        # either no mapping or already in steps
        validate_topology(topo, unit="steps" if unit0 == "steps" else unit0)

    topo, _ = annotate_topology_phyllotaxy(topo, phyllotaxy_config=phyllotaxy_config, seed=phyllotaxy_seed)
    return topo


def save_topology_json(topology: Topology, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(topology, indent=2, sort_keys=True), encoding="utf-8")


def load_topology_json(
    path: Path,
    *,
    mapping: Optional[Dict[str, Any]] = None,
    seed: Optional[int] = None,
    phyllotaxy_config: Optional[Dict[str, Any]] = None,
    phyllotaxy_seed: Optional[int] = None,
) -> Topology:
    topo: Topology = json.loads(path.read_text(encoding="utf-8"))
    if "root_id" not in topo or "branches" not in topo:
        raise ValueError("Invalid topology JSON: must contain root_id and branches.")
    # Validate before mapping in the appropriate unit
    unit0 = str(mapping.get("unit", "steps")).lower() if mapping is not None else "steps"
    validate_topology(topo, unit=unit0)

    if mapping is not None:
        topo = map_topology_to_sr_steps(
            topo,
            unit=str(mapping.get("unit", "steps")),
            rate=float(mapping.get("rate")),
            mode=str(mapping.get("mode", "poisson")),
            seed=seed,
            keep_observed=True,
        )
    # Validate mapped topology
    validate_topology(topo, unit='steps' if unit0 != 'steps' else unit0)
    topo, _ = annotate_topology_phyllotaxy(topo, phyllotaxy_config=phyllotaxy_config, seed=phyllotaxy_seed)
    return topo


def convert_csv_to_json(
    branches_csv: Path,
    events_csv: Path,
    out_json: Path,
    *,
    mapping: Optional[Dict[str, Any]] = None,
    seed: Optional[int] = None,
    phyllotaxy_config: Optional[Dict[str, Any]] = None,
    phyllotaxy_seed: Optional[int] = None,
) -> Topology:
    topo = read_topology_csv(
        branches_csv,
        events_csv,
        mapping=mapping,
        seed=seed,
        phyllotaxy_config=phyllotaxy_config,
        phyllotaxy_seed=phyllotaxy_seed,
    )
    save_topology_json(topo, out_json)
    return topo


def _cli() -> None:
    p = argparse.ArgumentParser(description="Convert pipeline topology CSV files to JSON.")
    p.add_argument("--branches_csv", type=Path, required=True)
    p.add_argument("--events_csv", type=Path, required=True)
    p.add_argument("--out_json", type=Path, required=True)
    p.add_argument("--topology_unit", type=str, default="steps", choices=["steps", "years", "meters"], help="Unit of T and event times in CSV.")
    p.add_argument("--rate", type=float, default=None, help="Turnover rate (events/year or events/meter). Required if topology_unit != steps.")
    p.add_argument("--mapping_mode", type=str, default="poisson", choices=["poisson", "deterministic"], help="How to map observed branch lengths to SR steps.")
    p.add_argument("--seed", type=int, default=None)
    args = p.parse_args()

    mapping = None
    if args.topology_unit != "steps":
        if args.rate is None:
            raise ValueError("--rate is required when --topology_unit is years or meters")
        mapping = {"unit": args.topology_unit, "rate": float(args.rate), "mode": args.mapping_mode}

    topo = convert_csv_to_json(args.branches_csv, args.events_csv, args.out_json, mapping=mapping, seed=args.seed)
    print(f"Wrote topology JSON with root_id={topo['root_id']} to {args.out_json}")


if __name__ == "__main__":
    _cli()
