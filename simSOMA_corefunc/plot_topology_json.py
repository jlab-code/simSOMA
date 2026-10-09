#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Set

try:
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
except ModuleNotFoundError as exc:
    raise SystemExit("matplotlib is required. Install it with: python -m pip install matplotlib") from exc


def fnum(x, default=None):
    try:
        if x is None:
            return default
        return float(x)
    except Exception:
        return default


def parse_ticks(s):
    if s is None:
        return None
    text = str(s).strip()
    if not text or text.lower() in {"auto", "none", "null"}:
        return None
    return [float(x.strip()) for x in text.split(",") if x.strip()]


def topology_axis_label(topo, override: str = "auto") -> str:
    if override and str(override).strip().lower() not in {"auto", ""}:
        return str(override)
    unit = str(topo.get("unit", "")).strip().lower()
    if unit in {"year", "years", "yr", "yrs"}:
        return "Age (years)"
    if unit in {"meter", "meters", "metre", "metres", "m"}:
        return "Position (meters)"
    if unit in {"step", "steps", "division", "divisions", "self-renewal steps", "self_renewal_steps"}:
        return "Self-renewal steps"
    if unit:
        return f"Topology coordinate ({topo.get('unit')})"
    return "Topology coordinate"


def apply_y_ticks(ax, tick_text, ymin: float, ymax: float):
    """Apply manual y ticks only if they fit the plotted coordinate range.

    Matplotlib expands the visible axis limits when set_yticks includes values
    outside the current limits. That made small tutorial topologies with fixed
    ticks such as 0,100,200 appear compressed near the bottom of the plot.
    Here, unsuitable manual ticks are ignored and Matplotlib chooses compact
    automatic ticks instead.
    """
    yt = parse_ticks(tick_text)
    if yt is not None:
        eps = max(1e-9, 1e-9 * max(abs(ymin), abs(ymax), 1.0))
        kept = [t for t in yt if (ymin - eps) <= t <= (ymax + eps)]
        if len(kept) >= 2:
            ax.set_yticks(kept)
            return {
                "requested_y_ticks": yt,
                "applied_y_ticks": kept,
                "mode": "manual_filtered" if len(kept) != len(yt) else "manual",
            }
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5, min_n_ticks=3))
    return {
        "requested_y_ticks": yt,
        "applied_y_ticks": None,
        "mode": "auto",
    }


def fill_missing_coordinates(topo) -> None:
    """Derive start_age / end_age (branches) and age (events) where a topology gives only
    lengths and relative positions (pos, or branch-local time). Existing values are kept;
    start / end / position are accepted as aliases."""
    branches = {str(b["id"]): b for b in topo.get("branches", []) if isinstance(b, dict) and "id" in b}
    events = [e for e in topo.get("events", []) if isinstance(e, dict)]
    for b in branches.values():
        if b.get("start_age") is None and b.get("start") is not None:
            b["start_age"] = b["start"]
        if b.get("end_age") is None and b.get("end") is not None:
            b["end_age"] = b["end"]
    for e in events:
        if e.get("age") is None and e.get("position") is not None:
            e["age"] = e["position"]

    def length(b):
        if b.get("length") is not None:
            return fnum(b.get("length"), 0.0)
        if b.get("start_age") is not None and b.get("end_age") is not None:
            return fnum(b["end_age"], 0.0) - fnum(b["start_age"], 0.0)
        return 0.0

    def rel_pos(e, b):
        if e.get("pos") is not None:
            return fnum(e["pos"], 0.0)
        if e.get("time") is not None and length(b) > 0:
            return fnum(e["time"], 0.0) / length(b)
        return None

    attach = {str(e.get("target")): e for e in events if e.get("type") == "branch"}
    resolving = set()

    def start(bid):
        b = branches[bid]
        if b.get("start_age") is not None:
            return fnum(b["start_age"], 0.0)
        if bid in resolving:            # cycle: leave unresolved
            return 0.0
        resolving.add(bid)
        parent = b.get("parent")
        ev = attach.get(bid)
        if parent is None or str(parent) not in branches:
            s = 0.0
        else:
            pb = branches[str(parent)]
            pos = rel_pos(ev, pb) if ev is not None else None
            s = start(str(parent)) + (pos if pos is not None else 1.0) * length(pb)
        b["start_age"] = s
        resolving.discard(bid)
        return s

    for bid, b in branches.items():
        s = start(bid)
        if b.get("end_age") is None:
            b["end_age"] = s + length(b)
    for e in events:
        if e.get("age") is None:
            bid = str(e.get("branch") or e.get("parent") or "")
            if bid in branches:
                pos = rel_pos(e, branches[bid])
                if pos is not None:
                    e["age"] = fnum(branches[bid]["start_age"], 0.0) + pos * length(branches[bid])


def load_topology(path: Path):
    if path.is_dir():
        raise SystemExit(
            f"Topology input is a directory, not a JSON file: {path}\n"
            "Pass the complete topology .json path."
        )
    if not path.exists():
        raise SystemExit(f"Topology JSON not found: {path}")
    with path.open() as fh:
        topo = json.load(fh)
    fill_missing_coordinates(topo)
    branches = {str(b["id"]): b for b in topo.get("branches", []) if isinstance(b, dict) and "id" in b}
    children = defaultdict(list)
    roots = []
    for bid, br in branches.items():
        parent = br.get("parent")
        if parent in branches:
            children[str(parent)].append(bid)
        else:
            roots.append(bid)
    for parent in children:
        children[parent].sort(key=lambda c: (fnum(branches[c].get("start_age"), 0.0), c))
    roots.sort(key=lambda c: (fnum(branches[c].get("start_age"), 0.0), c))
    return topo, branches, children, roots


def organ_events_by_branch(topo, branches, terminal_tol: float):
    terminal = defaultdict(list)
    nonterminal = defaultdict(list)
    ignored = []
    for ev in topo.get("events", []):
        if not isinstance(ev, dict) or ev.get("type") != "organ":
            continue
        bid = str(ev.get("branch") or ev.get("parent") or "")
        if bid not in branches:
            ignored.append({"reason": "organ event references missing branch", "event": ev})
            continue
        br = branches[bid]
        age = fnum(ev.get("age", ev.get("attach_age")), fnum(br.get("end_age"), 0.0))
        end_age = fnum(br.get("end_age"), age)
        organ = str(ev.get("organ") or ev.get("target") or "")
        if not organ:
            ignored.append({"reason": "organ event lacks organ/target id", "event": ev})
            continue
        is_terminal = bool(ev.get("terminal_tip", False)) or abs(age - end_age) <= terminal_tol
        rec = {
            "organ_id": organ,
            "branch_id": bid,
            "age": age,
            "branch_order": ev.get("branch_order", br.get("branch_order")),
            "terminal_tip": is_terminal,
        }
        if is_terminal:
            terminal[bid].append(rec)
        else:
            nonterminal[bid].append(rec)
    for d in (terminal, nonterminal):
        for bid in d:
            d[bid].sort(key=lambda r: (r["age"], r["organ_id"]))
    return terminal, nonterminal, ignored


def count_sampled_descendants(branches, children, organs_by_branch):
    memo = {}

    def rec(bid):
        if bid in memo:
            return memo[bid]
        n = len(organs_by_branch.get(bid, []))
        for child in children.get(bid, []):
            n += rec(child)
        memo[bid] = n
        return n

    for bid in branches:
        rec(bid)
    return memo


def allocate_x(branches, children, roots, keep, organs_by_branch):
    branch_x = {}
    organ_x = {}
    leaf_order = []
    next_x = 0.0

    def rec(bid):
        nonlocal next_x
        kept_children = [c for c in children.get(bid, []) if c in keep]
        child_roots = [rec(child) for child in kept_children]
        own_organs = organs_by_branch.get(bid, [])
        if own_organs:
            xs = []
            for org in own_organs:
                x = next_x
                next_x += 1.0
                organ_x[org["organ_id"]] = x
                leaf_order.append(org["organ_id"])
                xs.append(x)
            branch_x[bid] = sum(xs) / len(xs)
        elif child_roots:
            branch_x[bid] = sum(child_roots) / len(child_roots)
        else:
            branch_x[bid] = next_x
            next_x += 1.0
        return branch_x[bid]

    for root in roots:
        if root in keep:
            rec(root)
    return branch_x, organ_x, leaf_order


def descendants(branch_id: str, children) -> Set[str]:
    out = {branch_id}
    stack = list(children.get(branch_id, []))
    while stack:
        x = stack.pop()
        out.add(x)
        stack.extend(children.get(x, []))
    return out


def shift_subtree(branch_id: str, delta: float, children, branch_x: Dict[str, float], organ_x: Dict[str, float], organs_by_branch):
    for bid in descendants(branch_id, children):
        if bid in branch_x:
            branch_x[bid] += delta
        for org in organs_by_branch.get(bid, []):
            oid = org.get("organ_id")
            if oid in organ_x:
                organ_x[oid] += delta


def avoid_exact_parent_child_overlap(branches, children, keep, branch_x, organ_x, organs_by_branch, min_offset: float):
    """Move visible child subtrees slightly when child and parent would overlap exactly.

    This prevents a lateral branch from being drawn on top of the parent axis, which
    otherwise creates a misleading short tick that looks like a branch that goes nowhere.
    """
    if min_offset <= 0:
        return []
    changes = []
    for parent in sorted(keep, key=lambda b: (fnum(branches[b].get("start_age"), 0.0), b)):
        xp = branch_x.get(parent)
        if xp is None:
            continue
        visible_children = [c for c in children.get(parent, []) if c in keep]
        for i, child in enumerate(visible_children):
            xc = branch_x.get(child)
            if xc is None:
                continue
            if abs(xc - xp) < min_offset:
                # Alternate directions; choose right for the first exact-overlap child.
                side = 1.0 if (i % 2 == 0) else -1.0
                delta = side * (min_offset - (xc - xp))
                shift_subtree(child, delta, children, branch_x, organ_x, organs_by_branch)
                changes.append({"child": child, "parent": parent, "delta": delta})
    return changes


def max_relevant_age(bid, branches, children, keep, terminal_organs, nonterminal_organs, show_nonterminal_organs):
    ages = []
    ages.extend(o["age"] for o in terminal_organs.get(bid, []))
    if show_nonterminal_organs:
        ages.extend(o["age"] for o in nonterminal_organs.get(bid, []))
    ages.extend(
        fnum(branches[c].get("start_age"), fnum(branches[bid].get("start_age"), 0.0))
        for c in children.get(bid, [])
        if c in keep
    )
    if ages:
        return max(ages)
    return fnum(branches[bid].get("end_age"), fnum(branches[bid].get("start_age"), 0.0))


def make_order_colors(branches, keep, mode: str):
    if mode == "black":
        return {}, []
    orders = sorted({branches[b].get("branch_order") for b in keep if isinstance(branches[b].get("branch_order"), int)})
    cmap = plt.get_cmap("tab10")
    return {order: cmap(i % 10) for i, order in enumerate(orders)}, orders


def main(argv=None):
    ap = argparse.ArgumentParser(description="Axis-preserving pruned developmental topology plot for simSOMA JSON.")
    ap.add_argument("--topology_json", "--input", dest="topology_json", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--out_name", default="topology_developmental_pruned")
    ap.add_argument("--fig_width", type=float, default=2.6)
    ap.add_argument("--fig_height", type=float, default=1.8)
    ap.add_argument("--dpi", type=int, default=1000)
    ap.add_argument("--tick_size", type=float, default=6)
    ap.add_argument("--label_size", type=float, default=7)
    ap.add_argument("--tick_length", type=float, default=1.5)
    ap.add_argument("--tick_width", type=float, default=0.45)
    ap.add_argument("--tick_pad", type=float, default=2.0)
    ap.add_argument("--axis_linewidth", type=float, default=0.45)
    ap.add_argument("--terminal_tol", type=float, default=1e-5)
    ap.add_argument("--y_ticks", default="auto", help="Comma-separated y-axis ticks, or auto. Out-of-range tick sets are ignored to avoid compressing small topologies.")
    ap.add_argument("--y_label", default="auto", help="Y-axis label. Use auto for a unit-aware label from topology unit.")
    ap.add_argument("--branch_linewidth", type=float, default=0.85)
    ap.add_argument("--connector_linewidth", type=float, default=0.65)
    ap.add_argument("--branchpoint_size", type=float, default=6)
    ap.add_argument("--tip_size", type=float, default=5.0)
    ap.add_argument("--organ_size", type=float, default=10)
    ap.add_argument("--show_tip_labels", action="store_true")
    ap.add_argument("--show_nonterminal_organs", action="store_true", help="Draw sampled organ events that occur before the branch end; useful for organs along a trunk.")
    ap.add_argument("--organ_stub_length", type=float, default=0.28, help="Horizontal length of non-terminal organ stubs.")
    ap.add_argument("--nonterminal_organ_size", type=float, default=None, help="Marker size for non-terminal organs; defaults to --organ_size.")
    ap.add_argument("--nonterminal_organ_style", choices=["dotted", "dashed", "solid"], default="dotted")
    ap.add_argument("--nonterminal_organ_color", default="#595959")
    ap.add_argument("--min_branch_offset", type=float, default=0.35, help="Minimum x-offset between parent axis and visible child axis; prevents hidden/overlapping lateral branches.")
    ap.add_argument("--branch_color_mode", choices=["order", "black"], default="order")
    ap.add_argument("--title", default="")
    args = ap.parse_args(argv)

    topo, branches, children, roots = load_topology(args.topology_json)
    terminal_organs, nonterminal_organs, ignored_organs = organ_events_by_branch(topo, branches, args.terminal_tol)

    # Pruning is based on terminal sampled organs, because that is the original
    # publication-panel intent. If requested, non-terminal organs on kept branches
    # are overlaid as organ emergence events rather than branch-like axes.
    sampled_count = count_sampled_descendants(branches, children, terminal_organs)
    keep = {bid for bid, n in sampled_count.items() if n > 0}
    if not keep and args.show_nonterminal_organs:
        # For pure trunk-with-nonterminal-organs cases, keep branches that carry
        # nonterminal organ events even if no terminal organ was flagged.
        keep = {bid for bid, orgs in nonterminal_organs.items() if orgs}
    if not keep:
        raise SystemExit("No sampled terminal organs found. Use --show_nonterminal_organs if this topology samples organs along a non-terminal axis.")

    branch_x, organ_x, leaf_order = allocate_x(branches, children, roots, keep, terminal_organs)
    overlap_changes = avoid_exact_parent_child_overlap(
        branches, children, keep, branch_x, organ_x, terminal_organs, args.min_branch_offset
    )
    colors, orders = make_order_colors(branches, keep, args.branch_color_mode)

    def color_of(bid):
        if args.branch_color_mode == "black":
            return "black"
        order = branches[bid].get("branch_order")
        return colors.get(order, "black")

    fig, ax = plt.subplots(figsize=(args.fig_width, args.fig_height))
    rows = []

    # Draw vertical branch axes.
    for bid in sorted(keep, key=lambda b: (branches[b].get("branch_order", 999), fnum(branches[b].get("start_age"), 0.0), b)):
        br = branches[bid]
        x = branch_x[bid]
        y0 = fnum(br.get("start_age"), 0.0)
        y1 = max_relevant_age(bid, branches, children, keep, terminal_organs, nonterminal_organs, args.show_nonterminal_organs)
        ax.plot([x, x], [y0, y1], color=color_of(bid), linewidth=args.branch_linewidth, solid_capstyle="butt", zorder=2)
        if args.tip_size > 0:
            ax.scatter([x], [y1], marker="o", s=args.tip_size, facecolor=color_of(bid), edgecolor=color_of(bid), linewidth=0.2, zorder=3)
        rows.append({
            "type": "branch",
            "id": bid,
            "parent": br.get("parent"),
            "branch_order": br.get("branch_order"),
            "start_age": y0,
            "draw_end_age": y1,
            "true_end_age": br.get("end_age"),
            "x": x,
            "n_terminal_organs_in_subtree": sampled_count.get(bid, 0),
        })

    # Draw true branch connections, one per sampled-bearing lateral branch.
    for parent in sorted(keep, key=lambda b: (fnum(branches[b].get("start_age"), 0.0), b)):
        xp = branch_x[parent]
        for child in children.get(parent, []):
            if child not in keep:
                continue
            xc = branch_x[child]
            y = fnum(branches[child].get("start_age"), fnum(branches[parent].get("start_age"), 0.0))
            child_color = color_of(child)
            ax.plot([xp, xc], [y, y], color=child_color, linewidth=args.connector_linewidth, solid_capstyle="butt", zorder=2)
            ax.scatter([xp], [y], marker="o", s=args.branchpoint_size, facecolor="#4d4d4d", edgecolor="#4d4d4d", linewidth=0.25, zorder=4)
            rows.append({
                "type": "branch_event",
                "id": child,
                "parent": parent,
                "branch_order": branches[child].get("branch_order"),
                "start_age": y,
                "draw_end_age": "",
                "true_end_age": branches[child].get("end_age"),
                "x": xc,
                "n_terminal_organs_in_subtree": sampled_count.get(child, 0),
            })

    # Draw sampled terminal organs.
    for bid, orgs in terminal_organs.items():
        if bid not in keep:
            continue
        for org in orgs:
            x = branch_x[bid]
            y = org["age"]
            ax.scatter([x], [y], marker="s", s=args.organ_size, facecolor=color_of(bid), edgecolor="black", linewidth=0.35, zorder=5)
            if args.show_tip_labels:
                ax.text(x, y + 0.015 * max(1.0, y), org["organ_id"], fontsize=args.tick_size, ha="center", va="bottom", rotation=90)
            rows.append({
                "type": "terminal_organ",
                "id": org["organ_id"],
                "parent": bid,
                "branch_order": org.get("branch_order"),
                "start_age": y,
                "draw_end_age": "",
                "true_end_age": y,
                "x": x,
                "n_terminal_organs_in_subtree": "",
            })

    # Draw non-terminal organ events as organ stubs, not branches.
    linestyle = {"dotted": (0, (1.0, 1.4)), "dashed": (0, (2.4, 1.6)), "solid": "solid"}[args.nonterminal_organ_style]
    nonterm_drawn = 0
    if args.show_nonterminal_organs:
        nonterm_size = args.nonterminal_organ_size if args.nonterminal_organ_size is not None else args.organ_size
        for bid, orgs in nonterminal_organs.items():
            if bid not in keep:
                continue
            x0 = branch_x[bid]
            for i, org in enumerate(orgs):
                # Alternate sides to avoid giving the impression of a single long branch.
                side = 1.0 if i % 2 == 0 else -1.0
                x1 = x0 + side * args.organ_stub_length
                y = org["age"]
                ax.plot(
                    [x0, x1], [y, y],
                    color=args.nonterminal_organ_color,
                    linewidth=max(0.35, args.connector_linewidth * 0.75),
                    linestyle=linestyle,
                    solid_capstyle="butt",
                    dash_capstyle="round",
                    zorder=3,
                )
                ax.scatter([x1], [y], marker="o", s=nonterm_size, facecolor="white", edgecolor=args.nonterminal_organ_color, linewidth=0.55, zorder=5)
                if args.show_tip_labels:
                    ha = "left" if side > 0 else "right"
                    ax.text(x1 + side * 0.035, y, org["organ_id"], fontsize=args.tick_size, ha=ha, va="center")
                rows.append({
                    "type": "nonterminal_organ",
                    "id": org["organ_id"],
                    "parent": bid,
                    "branch_order": org.get("branch_order"),
                    "start_age": y,
                    "draw_end_age": "",
                    "true_end_age": y,
                    "x": x1,
                    "n_terminal_organs_in_subtree": "",
                })
                nonterm_drawn += 1

    all_y = [fnum(topo.get("tree_age"), 1.0), 1.0]
    all_y.extend(o["age"] for orgs in terminal_organs.values() for o in orgs)
    if args.show_nonterminal_organs:
        all_y.extend(o["age"] for orgs in nonterminal_organs.values() for o in orgs)
    max_age = max(all_y)

    all_x = list(branch_x.values())
    all_x.extend([float(r["x"]) for r in rows if r.get("x") not in (None, "")])
    if not all_x:
        all_x = [0.0]
    xmin, xmax = min(all_x), max(all_x)
    pad = max(0.45, 0.08 * max(1.0, xmax - xmin))

    ymin = 0.0
    ymax = max(max_age * 1.03, max_age + 0.5 if max_age < 10 else max_age * 1.03)
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(xmin - pad, xmax + pad)
    ax.set_xlabel("")
    ax.set_ylabel(topology_axis_label(topo, args.y_label), fontsize=args.label_size)
    ax.set_xticks([])
    y_tick_report = apply_y_ticks(ax, args.y_ticks, ymin, ymax)
    ax.set_ylim(ymin, ymax)
    ax.tick_params(axis="y", labelsize=args.tick_size, length=args.tick_length, width=args.tick_width, pad=args.tick_pad)
    ax.grid(True, axis="y", alpha=0.25)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    for spine in ["left", "bottom"]:
        ax.spines[spine].set_linewidth(args.axis_linewidth)
    if args.title:
        ax.set_title(args.title, fontsize=args.label_size)
    fig.subplots_adjust(left=0.18, right=0.98, bottom=0.10, top=0.98)

    args.outdir.mkdir(parents=True, exist_ok=True)
    out_png = args.outdir / f"{args.out_name}.png"
    out_pdf = args.outdir / f"{args.out_name}.pdf"
    out_csv = args.outdir / f"{args.out_name}_layout.csv"
    out_report = args.outdir / f"{args.out_name}_report.json"
    fig.savefig(out_png, dpi=args.dpi)
    fig.savefig(out_pdf, dpi=args.dpi)
    plt.close(fig)

    fieldnames = ["type", "id", "parent", "branch_order", "start_age", "draw_end_age", "true_end_age", "x", "n_terminal_organs_in_subtree"]
    with out_csv.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    report = {
        "plot_type": "axis_preserving_pruned_developmental_topology",
        "n_branches_total": len(branches),
        "n_branches_drawn": len(keep),
        "n_terminal_organs_drawn": sum(len(v) for v in terminal_organs.values()),
        "n_nonterminal_organs_available": sum(len(v) for v in nonterminal_organs.values()),
        "n_nonterminal_organs_drawn": nonterm_drawn,
        "branch_orders_colored": orders,
        "pruned_branches": sorted([bid for bid in branches if bid not in keep]),
        "overlap_adjustments": overlap_changes,
        "topology_unit": topo.get("unit"),
        "y_axis_label": topology_axis_label(topo, args.y_label),
        "y_axis_limits": [ymin, ymax],
        "y_ticks": y_tick_report,
        "note": "True branches are drawn as solid axes/connectors. Non-terminal organ events can be overlaid as short dotted/dashed organ stubs with circular markers, so organ emergence is visually distinct from branch formation. Visible child branches that would exactly overlap their parent axis are offset slightly to avoid hidden branches."
    }
    out_report.write_text(json.dumps(report, indent=2))

    print(f"Wrote: {out_png}")
    print(f"Wrote: {out_pdf}")
    print(f"Wrote: {out_csv}")
    print(f"Wrote: {out_report}")
    print(f"Branches drawn: {len(keep)} / {len(branches)}")
    print(f"Terminal organs drawn: {sum(len(v) for v in terminal_organs.values())}")
    print(f"Non-terminal organs drawn: {nonterm_drawn}")
    print(f"Pruned branches: {', '.join(report['pruned_branches']) if report['pruned_branches'] else 'none'}")


if __name__ == "__main__":
    main()
