"""inspect_topology.py

Standalone topology inspection utility.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, Optional

import topology_io
from pipeline_wrapper import validate_topology_contract, summarize_topology


def _mapping_from_args(args: argparse.Namespace) -> Optional[Dict[str, Any]]:
    unit = str(args.topology_unit).lower()
    if unit == "steps":
        return None
    if args.kappa_sr is None:
        raise SystemExit("Error: --kappa_sr is required when --topology_unit is years or meters")
    return {
        "unit": unit,
        "rate": float(args.kappa_sr),
        "mode": str(args.topology_mapping_mode).lower(),
    }


def _load_topology(args: argparse.Namespace) -> Dict[str, Any]:
    mapping = _mapping_from_args(args)
    seed = args.mapping_seed

    if args.topology_json is not None:
        return topology_io.load_topology_auto(Path(args.topology_json), mapping=mapping, seed=seed)

    if args.branches_csv is None or args.events_csv is None:
        raise SystemExit("Error: provide either --topology_json OR both --branches_csv and --events_csv")

    return topology_io.read_topology_csv(
        Path(args.branches_csv),
        Path(args.events_csv),
        mapping=mapping,
        seed=seed,
    )


def _write_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, sort_keys=True), encoding="utf-8")


def main() -> None:
    p = argparse.ArgumentParser(description="Inspect (load/convert/map/validate) a topology without running the simulation.")
    inp = p.add_argument_group("Input")
    inp.add_argument("--branches_csv", type=str, default=None, help="Path to branches.csv")
    inp.add_argument("--events_csv", type=str, default=None, help="Path to events.csv")
    inp.add_argument("--topology_json", type=str, default=None, help="Path to topology JSON (internal or user-friendly)")

    mapg = p.add_argument_group("Mapping (observed units -> SR steps)")
    mapg.add_argument("--topology_unit", type=str, default="steps", choices=["steps", "years", "meters"], help="Unit in the input topology. If years/meters, mapping is applied using --kappa_sr.")
    mapg.add_argument("--kappa_sr", type=float, default=None, help="SR divisions per topology unit used for mapping.")
    mapg.add_argument("--topology_mapping_mode", type=str, default="deterministic", choices=["deterministic"], help="Mapping mode for branch lengths (deterministic only)")
    mapg.add_argument("--mapping_seed", type=int, default=123, help="Seed for stochastic mapping (Poisson).")

    out = p.add_argument_group("Output")
    out.add_argument("--export_topology_json", type=str, default=None, help="If set, write the validated internal topology JSON to this path.")
    out.add_argument("--report_json", type=str, default=None, help="If set, write a report JSON (summary + topology) to this path.")
    out.add_argument("--print_topology", action="store_true", help="Print the full internal topology JSON to stdout")

    args = p.parse_args()

    topo = _load_topology(args)
    validate_topology_contract(topo)
    summary = summarize_topology(topo)

    print(json.dumps({"topology_summary": summary}, indent=2, sort_keys=True))
    if args.print_topology:
        print(json.dumps({"topology": topo}, indent=2, sort_keys=True))
    if args.export_topology_json is not None:
        topology_io.save_topology_json(topo, Path(args.export_topology_json))
    if args.report_json is not None:
        report = {
            "inputs": {
                "branches_csv": args.branches_csv,
                "events_csv": args.events_csv,
                "topology_json": args.topology_json,
                "topology_unit": args.topology_unit,
                "kappa_sr": args.kappa_sr,
                "topology_mapping_mode": args.topology_mapping_mode,
                "mapping_seed": args.mapping_seed,
            },
            "topology_summary": summary,
            "topology": topo,
        }
        _write_json(Path(args.report_json), report)


if __name__ == "__main__":
    main()
