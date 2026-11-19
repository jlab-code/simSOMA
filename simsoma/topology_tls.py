"""Convert a terrestrial-laser-scan (TLS) segment table into a simSOMA topology.

Input: a tab-separated segment table as written by quantitative structure models such as
TreeQSM, one row per segment, with at least the columns

    segment, parent_segment (0 = root), branch_order, length_m, base_distance_m

(other columns, e.g. diameters, height_m, azimuth_deg, are optional and partly kept as
metadata). A segment is a piece of a branch between two branching points.

Conversion (based on the conversion script written with ChatGPT for earlier TLS runs;
reviewed and extended for simSOMA 0.2.0):

1. Axes. Chains of segments of the same branch order are merged into one branch axis (a
   simSOMA branch). A new axis starts at the root and at every segment whose branch order is
   one higher than its parent's. If a segment has several children of its own order (not
   expected for TreeQSM output, but possible), the child with the largest base diameter
   (then the longest) continues the axis and the others start new axes of the same order.
2. Branch events. A child axis is attached to its parent axis at
       pos = (base_distance(child start) - base_distance(parent axis start)) / length(parent axis),
   i.e. by the path distance from the tree base, which is consistent along axes in TreeQSM
   tables. (The named parent segment is not used for the position: in real scans many child
   branches have base distances before or beyond the end of their named parent segment.) Positions are clipped to [0, 1]; the number of clipped events is
   reported.
3. Organs. One candidate organ at the tip (pos = 1) of every axis. Which candidates are kept is
   set by `organs` ("all", "random:<n>", "min_order:<k>", "orders:<k1,k2,...>", or an explicit
   list of source segment IDs of axis tips), optionally restricted to axes of at least
   `min_axis_length` m.
4. Pruning (default). Axes with no kept organ downstream are removed together with their
   branch events. This does not change the distribution of any simulated organ: a branch event
   leaves the parent niche unchanged, and removed sub-trees contain no sampled organ. Random
   draws differ from an unpruned run because self-renewal is split at fewer event times.

Lengths stay in meters; the configuration maps meters to self-renewal rounds
(topology.mapping_unit = "meters", mapping_rate = rounds per meter).
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import random
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Optional, Sequence

TLS_CONVERTER_VERSION = "1.0.0"
REQUIRED_COLUMNS = ("segment", "parent_segment", "branch_order", "length_m", "base_distance_m")


# ------------------------------------------------------------------------------------- reading
def read_segment_table(path: Path) -> dict[int, dict[str, Any]]:
    rows: dict[int, dict[str, Any]] = {}
    with Path(path).open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path}: no header row")
        reader.fieldnames = [f.strip() for f in reader.fieldnames]
        missing = [c for c in REQUIRED_COLUMNS if c not in reader.fieldnames]
        if missing:
            raise ValueError(f"{path}: missing required columns {missing} (need {list(REQUIRED_COLUMNS)})")
        for line, row in enumerate(reader, start=2):
            if not row.get("segment"):
                continue
            try:
                seg = int(float(row["segment"]))
                r = {"parent": int(float(row["parent_segment"])), "order": int(float(row["branch_order"])),
                     "length": float(row["length_m"]), "base": float(row["base_distance_m"]),
                     "diameter": _num(row.get("diameter_base_cm")), "height": _num(row.get("height_m")),
                     "azimuth": _num(row.get("azimuth_deg"))}
            except (TypeError, ValueError) as exc:
                raise ValueError(f"{path}, line {line}: cannot parse required fields ({exc})") from exc
            if seg in rows:
                raise ValueError(f"{path}, line {line}: duplicate segment {seg}")
            if r["length"] < 0 or not math.isfinite(r["length"]):
                raise ValueError(f"{path}, line {line}: invalid length_m {row['length_m']}")
            rows[seg] = r
    if not rows:
        raise ValueError(f"{path}: no segments")
    return rows


def _num(x: Any) -> Optional[float]:
    try:
        y = float(x)
        return y if math.isfinite(y) else None
    except (TypeError, ValueError):
        return None


# ------------------------------------------------------------------------------------ structure
def check_tree(rows: dict[int, dict[str, Any]]) -> tuple[int, dict[int, list[int]], list[str]]:
    children: dict[int, list[int]] = defaultdict(list)
    roots = []
    for seg, r in rows.items():
        if r["parent"] == 0:
            roots.append(seg)
        elif r["parent"] not in rows:
            raise ValueError(f"segment {seg}: parent_segment {r['parent']} not in table")
        else:
            children[r["parent"]].append(seg)
    if len(roots) != 1:
        raise ValueError(f"expected exactly one root segment (parent_segment = 0), found {len(roots)}: {sorted(roots)[:10]}")
    root = roots[0]
    seen, stack = set(), [root]
    while stack:
        s = stack.pop()
        if s in seen:
            raise ValueError(f"cycle at segment {s}")
        seen.add(s)
        stack.extend(children.get(s, []))
    if len(seen) != len(rows):
        raise ValueError(f"{len(rows) - len(seen)} segments not reachable from the root, e.g. {sorted(set(rows) - seen)[:10]}")
    warnings = []
    bad = [s for s, r in rows.items() if r["parent"] and r["order"] not in (rows[r["parent"]]["order"], rows[r["parent"]]["order"] + 1)]
    if bad:
        warnings.append(f"{len(bad)} segments with branch_order not equal to or one above their parent's (e.g. {bad[0]}); "
                        "they start new axes")
    return root, children, warnings


def build_axes(rows, children, root) -> tuple[list[dict[str, Any]], dict[int, int], int]:
    """Merge same-order chains into axes. Returns axes (start-ordered), segment -> axis index,
    and the number of segments with several same-order children."""
    def rank(s):
        r = rows[s]
        return (-(r["diameter"] or 0.0), -r["length"], s)

    starts, n_multi = [root], 0
    continuation: dict[int, int] = {}
    for seg in rows:
        same = [c for c in children.get(seg, []) if rows[c]["order"] == rows[seg]["order"]]
        if same:
            same.sort(key=rank)
            continuation[seg] = same[0]
            starts.extend(same[1:])
            n_multi += len(same) > 1
        for c in children.get(seg, []):
            if rows[c]["order"] != rows[seg]["order"]:
                starts.append(c)
    starts = sorted(set(starts), key=lambda s: (rows[s]["base"], s))
    seg_axis: dict[int, int] = {}
    axes = []
    for i, st in enumerate(starts):
        segs, cur = [], st
        while True:
            seg_axis[cur] = i
            segs.append(cur)
            if cur not in continuation:
                break
            cur = continuation[cur]
        axes.append({"segments": segs, "order": rows[st]["order"], "start": st, "tip": segs[-1],
                     "length": sum(rows[s]["length"] for s in segs), "base": rows[st]["base"]})
    assert len(seg_axis) == len(rows)
    return axes, seg_axis, n_multi


# ---------------------------------------------------------------------------------- conversion
def _parse_organ_policy(organs: Any, axes, rows, min_axis_length: float, seed: int) -> list[int]:
    cand = [i for i, a in enumerate(axes) if a["length"] >= min_axis_length]
    if isinstance(organs, (list, tuple, set)):
        tips = {a["tip"]: i for i, a in enumerate(axes)}
        unknown = [s for s in organs if int(s) not in tips]
        if unknown:
            raise ValueError(f"organ segments are not axis tips: {unknown[:10]}")
        return sorted(tips[int(s)] for s in organs)
    spec = str(organs)
    if spec == "all":
        return cand
    kind, _, arg = spec.partition(":")
    if kind == "random":
        n = int(arg)
        if n > len(cand):
            raise ValueError(f"random:{n} requested but only {len(cand)} candidate tips")
        return sorted(random.Random(seed).sample(cand, n))
    if kind == "min_order":
        return [i for i in cand if axes[i]["order"] >= int(arg)]
    if kind == "orders":
        keep = {int(x) for x in arg.split(",") if x}
        return [i for i in cand if axes[i]["order"] in keep]
    raise ValueError("organs must be 'all', 'random:<n>', 'min_order:<k>', 'orders:<k1,k2>' or a list of tip segment IDs")


def convert(rows: dict[int, dict[str, Any]], *, organs: Any = "all", min_axis_length: float = 0.0,
            prune: bool = True, seed: int = 0, source_name: str = "") -> tuple[dict[str, Any], dict[str, Any]]:
    """Return (topology JSON in the simSOMA user format, report)."""
    root, children, warnings = check_tree(rows)
    axes, seg_axis, n_multi = build_axes(rows, children, root)
    if n_multi:
        warnings.append(f"{n_multi} segments with several same-order children; extra children start new axes")
    parent = [None] * len(axes)
    pos = [None] * len(axes)
    n_clip = 0
    for i, a in enumerate(axes):
        p = rows[a["start"]]["parent"]
        if p == 0:
            continue
        j = seg_axis[p]
        P = axes[j]
        x = (a["base"] - P["base"]) / P["length"] if P["length"] > 0 else 0.0
        if x < 0 or x > 1:
            n_clip += 1
        parent[i], pos[i] = j, min(1.0, max(0.0, x))
    if n_clip:
        warnings.append(f"{n_clip} branch positions outside [0, 1] were clipped")

    kept_organs = _parse_organ_policy(organs, axes, rows, min_axis_length, seed)
    if not kept_organs:
        raise ValueError("no organs selected")
    keep = set(range(len(axes)))
    if prune:
        keep = set()
        for i in kept_organs:
            while i is not None and i not in keep:
                keep.add(i)
                i = parent[i]
    order = [i for i in range(len(axes)) if i in keep]
    bid = {i: f"B{k:05d}" for k, i in enumerate(order)}

    # start/end coordinates along the path from the tree base (m), consistent with lengths and
    # positions (start = parent start + pos * parent length); used by the topology plotter.
    start: list[Optional[float]] = [None] * len(axes)
    for i0 in range(len(axes)):
        path, i = [], i0
        while i is not None and start[i] is None:
            path.append(i); i = parent[i]
        base = 0.0 if i is None else start[i]
        for j in reversed(path):
            base = 0.0 if parent[j] is None else base + pos[j] * axes[parent[j]]["length"]
            start[j] = base
    branches, events = [], []
    organ_set = set(kept_organs)
    for i in order:
        a = axes[i]
        branches.append({"id": bid[i], "parent": bid[parent[i]] if parent[i] is not None else None,
                         "length": round(a["length"], 6), "start_age": round(start[i], 6),
                         "end_age": round(start[i] + a["length"], 6), "branch_order": a["order"],
                         "source_start_segment": a["start"], "source_tip_segment": a["tip"],
                         "source_n_segments": len(a["segments"])})
        if parent[i] is not None:
            events.append({"branch": bid[parent[i]], "pos": round(pos[i], 8), "type": "branch", "target": bid[i],
                           "age": round(start[i], 6)})
        if i in organ_set:
            events.append({"branch": bid[i], "pos": 1.0, "type": "organ", "target": f"O_{a['tip']}",
                           "age": round(start[i] + a["length"], 6)})
    events.sort(key=lambda e: (e["branch"], e["pos"], e["type"] != "branch", e["target"]))

    lengths = [b["length"] for b in branches]
    report = {
        "converter": "simSOMA topology_tls", "converter_version": TLS_CONVERTER_VERSION, "source_file": source_name,
        "n_segments": len(rows), "n_axes_total": len(axes), "n_candidate_organs": len(axes),
        "organ_policy": organs if not isinstance(organs, (list, tuple, set)) else f"list of {len(organs)} tip segments",
        "min_axis_length_m": min_axis_length, "pruned": bool(prune), "seed": seed,
        "n_branches": len(branches), "n_branch_events": sum(e["type"] == "branch" for e in events),
        "n_organs": sum(e["type"] == "organ" for e in events),
        "axis_order_counts": {str(k): v for k, v in sorted(_count(axes[i]["order"] for i in order).items())},
        "organ_order_counts": {str(k): v for k, v in sorted(_count(axes[i]["order"] for i in kept_organs).items())},
        "total_length_m": round(sum(lengths), 4), "max_axis_length_m": round(max(lengths), 4),
        "max_tip_distance_m": round(max(axes[i]["base"] + axes[i]["length"] for i in kept_organs), 4),
        "warnings": warnings,
    }
    topo = {"unit": "meters",
            "description": f"Converted from TLS segment table {source_name} (simSOMA topology_tls {TLS_CONVERTER_VERSION}); "
                           f"{report['n_organs']} organs at axis tips, policy {report['organ_policy']}.",
            "branches": branches, "events": events, "conversion": {k: v for k, v in report.items() if k != "warnings"}}
    return topo, report


def _count(xs: Iterable[int]) -> dict[int, int]:
    out: dict[int, int] = defaultdict(int)
    for x in xs:
        out[x] += 1
    return out


def write_axis_table(rows, path: Path) -> None:
    """Segment -> axis membership (for mapping results back to the scan)."""
    root, children, _ = check_tree(rows)
    axes, _, _ = build_axes(rows, children, root)
    with Path(path).open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["axis_index", "segment", "branch_order", "axis_start_segment", "axis_tip_segment"])
        for i, a in enumerate(axes):
            for s in a["segments"]:
                w.writerow([i, s, a["order"], a["start"], a["tip"]])


# ------------------------------------------------------------------------------------------ CLI
def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(prog="simsoma topology-from-tls",
                                 description="Convert a TLS (TreeQSM-style) segment table into a simSOMA topology (meters).")
    ap.add_argument("segments", type=Path, help="tab-separated segment table")
    ap.add_argument("out_json", type=Path)
    ap.add_argument("--organs", default="all",
                    help="all | random:<n> | min_order:<k> | orders:<k1,k2> | @file (one tip segment ID per line)")
    ap.add_argument("--min-axis-length", type=float, default=0.0, help="only axes at least this long (m) carry organs")
    ap.add_argument("--no-prune", action="store_true", help="keep axes without sampled organs downstream")
    ap.add_argument("--seed", type=int, default=0, help="seed for random organ selection")
    ap.add_argument("--report", type=Path, default=None, help="write the conversion report (JSON)")
    ap.add_argument("--axis-table", type=Path, default=None, help="write segment -> axis membership (CSV)")
    a = ap.parse_args(argv)
    rows = read_segment_table(a.segments)
    organs: Any = a.organs
    if organs.startswith("@"):
        organs = [int(x) for x in Path(organs[1:]).read_text().split() if x.strip()]
    topo, report = convert(rows, organs=organs, min_axis_length=a.min_axis_length, prune=not a.no_prune,
                           seed=a.seed, source_name=a.segments.name)
    a.out_json.parent.mkdir(parents=True, exist_ok=True)
    a.out_json.write_text(json.dumps(topo, indent=1) + "\n")
    if a.report:
        a.report.write_text(json.dumps(report, indent=2) + "\n")
    if a.axis_table:
        write_axis_table(rows, a.axis_table)
    print(json.dumps({k: report[k] for k in ("n_segments", "n_axes_total", "n_branches", "n_organs", "total_length_m",
                                              "warnings")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
