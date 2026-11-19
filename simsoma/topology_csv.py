#!/usr/bin/env python3
"""Convert branch/organ CSV topology tables to simSOMA user-friendly topology JSON.

This is a preprocessing helper. It does not change the simSOMA config schema:
users convert two CSV files to the standard topology JSON, then point
`topology.topology_json` in the usual config to the generated JSON.

Official CSV headers
--------------------
branches CSV:
    branch_id,parent_id,start,end

organs CSV:
    organ_id,branch_id,position

Accepted aliases are provided for convenience and backward compatibility:
    branches: id -> branch_id; parent -> parent_id; start_age -> start; end_age -> end
    organs:   id -> organ_id; branch -> branch_id; age -> position

The coordinates are unit-neutral. They can represent years, meters, or another
branch-axis coordinate. The generated topology JSON still needs a simSOMA unit
(`years`, `meters`, or `steps`) so the existing config mapping logic can be used.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from collections import defaultdict, deque
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

NULL_STRINGS = {"", "null", "none", "na", "nan", "."}

BRANCH_ALIASES = {
    "branch_id": ("branch_id", "id", "branch", "name"),
    "parent_id": ("parent_id", "parent", "parent_branch", "parent_branch_id"),
    "start": ("start", "start_age", "start_time", "start_position", "start_coordinate"),
    "end": ("end", "end_age", "end_time", "end_position", "end_coordinate"),
    "branch_order": ("branch_order", "order"),
}

ORGAN_ALIASES = {
    "organ_id": ("organ_id", "id", "organ", "name"),
    "branch_id": ("branch_id", "branch", "parent_branch", "parent_branch_id"),
    "position": ("position", "age", "time", "coordinate", "pos_abs"),
}

EVENT_TYPE_ORDER = {"branch": 0, "organ": 1}


def _norm_header(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(s).strip().lower()).strip("_")


def _is_null(value: Any) -> bool:
    if value is None:
        return True
    return str(value).strip().lower() in NULL_STRINGS


def _clean_str(value: Any, *, field: str, row_number: int) -> str:
    if value is None:
        raise ValueError(f"Missing value for {field} in row {row_number}.")
    out = str(value).strip()
    if not out:
        raise ValueError(f"Empty value for {field} in row {row_number}.")
    return out


def _clean_optional_str(value: Any) -> Optional[str]:
    if _is_null(value):
        return None
    return str(value).strip()


def _clean_float(value: Any, *, field: str, row_number: int) -> float:
    if value is None or str(value).strip() == "":
        raise ValueError(f"Missing numeric value for {field} in row {row_number}.")
    try:
        out = float(str(value).strip())
    except Exception as exc:
        raise ValueError(f"Could not parse {field}={value!r} as a number in row {row_number}.") from exc
    if not math.isfinite(out):
        raise ValueError(f"Non-finite value for {field} in row {row_number}: {value!r}.")
    return out


def _choose_columns(fieldnames: Sequence[str], aliases: Dict[str, Sequence[str]], *, table_name: str) -> Dict[str, Optional[str]]:
    normalized_to_original: Dict[str, str] = {}
    for name in fieldnames:
        key = _norm_header(name)
        if key in normalized_to_original:
            raise ValueError(
                f"{table_name} has ambiguous duplicate-like column headers after normalization: "
                f"{normalized_to_original[key]!r} and {name!r}."
            )
        normalized_to_original[key] = name

    out: Dict[str, Optional[str]] = {}
    missing: List[str] = []
    for canonical, choices in aliases.items():
        match = None
        for choice in choices:
            norm = _norm_header(choice)
            if norm in normalized_to_original:
                match = normalized_to_original[norm]
                break
        out[canonical] = match
        if match is None and canonical != "branch_order":
            missing.append(canonical)
    if missing:
        expected = {
            canonical: list(choices)
            for canonical, choices in aliases.items()
            if canonical in missing
        }
        raise ValueError(
            f"{table_name} is missing required columns {missing}. Accepted header aliases: {expected}."
        )
    return out


def _read_branch_rows(path: Path) -> List[Dict[str, Any]]:
    if not path.exists():
        raise FileNotFoundError(f"Branch CSV not found: {path}")
    rows: List[Dict[str, Any]] = []
    with path.open("r", newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        if not reader.fieldnames:
            raise ValueError(f"Branch CSV has no header row: {path}")
        cols = _choose_columns(reader.fieldnames, BRANCH_ALIASES, table_name="branch CSV")
        for row_number, row in enumerate(reader, start=2):
            if not any((v or "").strip() for v in row.values()):
                continue
            branch_id = _clean_str(row.get(cols["branch_id"]), field="branch_id", row_number=row_number)
            parent_id = _clean_optional_str(row.get(cols["parent_id"]))
            start = _clean_float(row.get(cols["start"]), field="start", row_number=row_number)
            end = _clean_float(row.get(cols["end"]), field="end", row_number=row_number)
            branch_order: Optional[int] = None
            order_col = cols.get("branch_order")
            if order_col is not None and not _is_null(row.get(order_col)):
                order_float = _clean_float(row.get(order_col), field="branch_order", row_number=row_number)
                order_int = int(round(order_float))
                if abs(order_float - order_int) > 1e-9:
                    raise ValueError(f"branch_order must be an integer-like value in row {row_number}.")
                branch_order = order_int
            rows.append(
                {
                    "branch_id": branch_id,
                    "parent_id": parent_id,
                    "start": start,
                    "end": end,
                    "branch_order": branch_order,
                    "row_number": row_number,
                }
            )
    if not rows:
        raise ValueError(f"Branch CSV contains no data rows: {path}")
    return rows


def _read_organ_rows(path: Path) -> List[Dict[str, Any]]:
    if not path.exists():
        raise FileNotFoundError(f"Organ CSV not found: {path}")
    rows: List[Dict[str, Any]] = []
    with path.open("r", newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        if not reader.fieldnames:
            raise ValueError(f"Organ CSV has no header row: {path}")
        cols = _choose_columns(reader.fieldnames, ORGAN_ALIASES, table_name="organ CSV")
        for row_number, row in enumerate(reader, start=2):
            if not any((v or "").strip() for v in row.values()):
                continue
            organ_id = _clean_str(row.get(cols["organ_id"]), field="organ_id", row_number=row_number)
            branch_id = _clean_str(row.get(cols["branch_id"]), field="branch_id", row_number=row_number)
            position = _clean_float(row.get(cols["position"]), field="position", row_number=row_number)
            rows.append(
                {
                    "organ_id": organ_id,
                    "branch_id": branch_id,
                    "position": position,
                    "row_number": row_number,
                }
            )
    if not rows:
        raise ValueError(f"Organ CSV contains no data rows: {path}")
    return rows


def _fmt_number(x: float) -> float | int:
    if abs(x - round(x)) <= 1e-12:
        return int(round(x))
    return float(x)


def _check_unique(values: Iterable[str], *, what: str) -> None:
    seen: set[str] = set()
    dup: List[str] = []
    for value in values:
        if value in seen:
            dup.append(value)
        seen.add(value)
    if dup:
        raise ValueError(f"Duplicate {what} identifier(s): {sorted(set(dup))}")


def _compute_branch_orders(branches: Dict[str, Dict[str, Any]], root_id: str) -> None:
    children: Dict[str, List[str]] = defaultdict(list)
    for bid, br in branches.items():
        parent = br["parent_id"]
        if parent is not None:
            children[parent].append(bid)
    for parent in children:
        children[parent].sort(key=lambda b: (branches[b]["start"], b))

    q: deque[Tuple[str, int]] = deque([(root_id, 0)])
    while q:
        bid, order = q.popleft()
        if branches[bid].get("branch_order") is None:
            branches[bid]["branch_order"] = order
        for child in children.get(bid, []):
            q.append((child, int(branches[bid]["branch_order"]) + 1))


def _validate_branch_graph(branches: Dict[str, Dict[str, Any]], *, tol: float) -> str:
    _check_unique(branches.keys(), what="branch")

    roots = [bid for bid, br in branches.items() if br["parent_id"] is None]
    if len(roots) != 1:
        raise ValueError(f"Expected exactly one root branch with empty/null parent_id. Found roots: {roots}")
    root_id = roots[0]

    for bid, br in branches.items():
        if br["end"] <= br["start"]:
            raise ValueError(
                f"Branch {bid!r} must have end > start. Observed start={br['start']}, end={br['end']}."
            )
        parent = br["parent_id"]
        if parent is not None and parent not in branches:
            raise ValueError(f"Branch {bid!r} references missing parent_id {parent!r}.")
        if parent is not None:
            parent_br = branches[parent]
            if br["start"] < parent_br["start"] - tol or br["start"] > parent_br["end"] + tol:
                raise ValueError(
                    f"Branch {bid!r} starts outside its parent interval: child start={br['start']}, "
                    f"parent {parent!r} interval=[{parent_br['start']}, {parent_br['end']}]."
                )

    children: Dict[str, List[str]] = defaultdict(list)
    for bid, br in branches.items():
        parent = br["parent_id"]
        if parent is not None:
            children[parent].append(bid)

    seen: set[str] = set()
    active: set[str] = set()

    def dfs(bid: str) -> None:
        if bid in active:
            raise ValueError(f"Cycle detected at branch {bid!r}.")
        if bid in seen:
            return
        active.add(bid)
        for child in children.get(bid, []):
            dfs(child)
        active.remove(bid)
        seen.add(bid)

    dfs(root_id)
    if len(seen) != len(branches):
        missing = sorted(set(branches) - seen)
        raise ValueError(f"Topology is disconnected; branches not reachable from root {root_id!r}: {missing}")

    _compute_branch_orders(branches, root_id)
    return root_id


def _validate_organs(organs: List[Dict[str, Any]], branches: Dict[str, Dict[str, Any]], *, tol: float) -> None:
    _check_unique((o["organ_id"] for o in organs), what="organ")
    for org in organs:
        bid = org["branch_id"]
        if bid not in branches:
            raise ValueError(f"Organ {org['organ_id']!r} references missing branch_id {bid!r}.")
        br = branches[bid]
        pos = org["position"]
        if pos < br["start"] - tol or pos > br["end"] + tol:
            raise ValueError(
                f"Organ {org['organ_id']!r} position is outside branch {bid!r}: "
                f"position={pos}, branch interval=[{br['start']}, {br['end']}]."
            )


def _fraction_along_branch(value: float, branch: Dict[str, Any]) -> float:
    denom = branch["end"] - branch["start"]
    frac = (value - branch["start"]) / denom
    # Clamp tiny floating-point excursions only.
    if frac < 0 and frac > -1e-12:
        frac = 0.0
    if frac > 1 and frac < 1 + 1e-12:
        frac = 1.0
    return float(frac)


def convert_branch_organ_csv_to_user_json(
    branches_csv: Path,
    organs_csv: Path,
    *,
    unit: str,
    description: Optional[str] = None,
    terminal_tol: float = 1e-9,
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """Convert branch/organ CSV files to the simSOMA user-friendly topology JSON schema."""
    unit = unit.lower()
    if unit not in {"steps", "years", "meters"}:
        raise ValueError("unit must be one of: steps, years, meters")

    branch_rows = _read_branch_rows(branches_csv)
    organ_rows = _read_organ_rows(organs_csv)

    branches: Dict[str, Dict[str, Any]] = {}
    for row in branch_rows:
        bid = row["branch_id"]
        if bid in branches:
            raise ValueError(f"Duplicate branch identifier {bid!r}.")
        branches[bid] = {
            "branch_id": bid,
            "parent_id": row["parent_id"],
            "start": row["start"],
            "end": row["end"],
            "branch_order": row["branch_order"],
            "row_number": row["row_number"],
        }

    root_id = _validate_branch_graph(branches, tol=terminal_tol)
    _validate_organs(organ_rows, branches, tol=terminal_tol)

    tree_max = max(br["end"] for br in branches.values())

    branch_entries: List[Dict[str, Any]] = []
    for bid in sorted(branches, key=lambda b: (branches[b]["start"], branches[b]["branch_order"], b)):
        br = branches[bid]
        length = br["end"] - br["start"]
        branch_entries.append(
            {
                "id": bid,
                "parent": br["parent_id"],
                "length": _fmt_number(length),
                # General, unit-neutral coordinate fields.
                "start": _fmt_number(br["start"]),
                "end": _fmt_number(br["end"]),
                # Backward-compatible aliases used by the bundled topology plotter.
                "start_age": _fmt_number(br["start"]),
                "end_age": _fmt_number(br["end"]),
                "branch_order": int(br["branch_order"]),
            }
        )

    events: List[Dict[str, Any]] = []

    # A child branch becomes a branch-initiation event on its parent.
    for bid, br in branches.items():
        parent = br["parent_id"]
        if parent is None:
            continue
        parent_br = branches[parent]
        pos = _fraction_along_branch(br["start"], parent_br)
        events.append(
            {
                "branch": parent,
                "pos": pos,
                "type": "branch",
                "target": bid,
                "parent": parent,
                "position": _fmt_number(br["start"]),
                "age": _fmt_number(br["start"]),
                "child": bid,
                "branch_order": int(br["branch_order"]),
            }
        )

    for org in organ_rows:
        bid = org["branch_id"]
        br = branches[bid]
        pos = _fraction_along_branch(org["position"], br)
        is_terminal = abs(org["position"] - br["end"]) <= terminal_tol
        events.append(
            {
                "branch": bid,
                "pos": pos,
                "type": "organ",
                "target": org["organ_id"],
                "organ": org["organ_id"],
                "position": _fmt_number(org["position"]),
                "age": _fmt_number(org["position"]),
                "terminal_tip": bool(is_terminal),
            }
        )

    events.sort(key=lambda e: (str(e["branch"]), float(e["pos"]), EVENT_TYPE_ORDER[str(e["type"])], str(e["target"])))

    topology = {
        "unit": unit,
        "coordinate_fields": {
            "branch_start": "start",
            "branch_end": "end",
            "organ_position": "position",
            "note": "Coordinates are unit-neutral; interpret them according to the topology unit and the simulation config mapping_unit/mapping_rate.",
        },
        "tree_coordinate_max": _fmt_number(tree_max),
        # Backward-compatible plotting alias. This remains useful even when unit="meters".
        "tree_age": _fmt_number(tree_max),
        "description": description or "Converted from branch and organ CSV topology files.",
        "source_files": {
            "branches_csv": str(branches_csv),
            "organs_csv": str(organs_csv),
        },
        "branches": branch_entries,
        "events": events,
    }

    report = {
        "root_id": root_id,
        "unit": unit,
        "n_branches": len(branch_entries),
        "n_organs": len(organ_rows),
        "n_events": len(events),
        "tree_coordinate_max": _fmt_number(tree_max),
        "branch_ids": [b["id"] for b in branch_entries],
        "organ_ids": [o["organ_id"] for o in organ_rows],
        "validation": "passed",
    }
    return topology, report


def write_json(path: Path, obj: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding="utf-8")


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Convert unit-neutral branch/organ CSV topology tables to simSOMA user-friendly topology JSON."
    )
    p.add_argument("--branches", "--branches-csv", dest="branches_csv", type=Path, required=True, help="Branch CSV with branch_id,parent_id,start,end columns.")
    p.add_argument("--organs", "--organs-csv", dest="organs_csv", type=Path, required=True, help="Organ CSV with organ_id,branch_id,position columns.")
    p.add_argument("--out", "--out-json", dest="out_json", type=Path, required=True, help="Output user-friendly topology JSON path.")
    p.add_argument("--unit", choices=["steps", "years", "meters"], default="years", help="Topology unit written into the JSON. Coordinates themselves are unit-neutral.")
    p.add_argument("--description", default=None, help="Optional description written into the output JSON.")
    p.add_argument("--terminal-tol", type=float, default=1e-9, help="Tolerance for flagging organ positions at branch end as terminal tips.")
    p.add_argument("--report", type=Path, default=None, help="Optional validation report JSON path.")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        topology, report = convert_branch_organ_csv_to_user_json(
            args.branches_csv,
            args.organs_csv,
            unit=args.unit,
            description=args.description,
            terminal_tol=args.terminal_tol,
        )
        write_json(args.out_json, topology)
        if args.report is not None:
            write_json(args.report, report)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    print(f"Wrote topology JSON: {args.out_json}")
    if args.report is not None:
        print(f"Wrote validation report: {args.report}")
    print(f"Branches: {report['n_branches']} | Organs: {report['n_organs']} | Events: {report['n_events']} | Root: {report['root_id']}")
    print("Next step: use the output JSON in topology.topology_json in the usual simSOMA config.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
