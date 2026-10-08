"""simSOMA command-line interface.

    simsoma run CONFIG [--splits N --jobs J] [--skip-check]   check topology, then simulate (grid workflow)
    simsoma check CONFIG                                       topology check / plots only
    simsoma layers CONFIG                                      per-layer simulation + read-level observation
    simsoma topology-from-csv BRANCHES ORGANS OUT [--unit years|meters|steps]
    simsoma template layered                                   print a layered-config template
    simsoma version
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import simsoma


def _versions() -> dict:
    out = {"simsoma": simsoma.__version__}
    try:
        import plantsoma_obs
        out["plantsoma_obs"] = plantsoma_obs.__version__
    except ImportError:
        out["plantsoma_obs"] = None
    out["core_path"] = str(simsoma.core_path())
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="simsoma", description="simSOMA: somatic VAF spectra from plant cell lineages")
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="topology check + grid simulation from a JSON config")
    r.add_argument("config", type=Path)
    r.add_argument("--splits", type=int, default=1, help="split the parameter grid into N sub-runs")
    r.add_argument("--jobs", type=int, default=1, help="concurrent sub-runs when --splits > 1")
    r.add_argument("--split-axis", choices=["parameter", "replicate"], default="parameter")
    r.add_argument("--skip-check", action="store_true", help="do not run the topology check first")
    c = sub.add_parser("check", help="topology check / plots only")
    c.add_argument("config", type=Path)
    l = sub.add_parser("layers", help="per-layer simulation and read-level observation (layered organs)")
    l.add_argument("config", type=Path)
    t = sub.add_parser("topology-from-csv", help="convert branch/organ CSV tables to topology JSON")
    t.add_argument("branches_csv", type=Path); t.add_argument("organs_csv", type=Path); t.add_argument("out_json", type=Path)
    t.add_argument("--unit", default="years", choices=["years", "meters", "steps"])
    t.add_argument("--report", type=Path, default=None)
    tp = sub.add_parser("template", help="print a config template")
    tp.add_argument("kind", choices=["layered"])
    sub.add_parser("version", help="print versions")
    a = ap.parse_args(argv)

    simsoma.use_core()
    os.environ.setdefault("SIMSOMA_SKIP_TOPOLOGY_CONFIRM", "1")
    try:
        if a.cmd == "version":
            print(json.dumps(_versions(), indent=2)); return 0
        if a.cmd == "template":
            import layered
            print(json.dumps(layered.LAYERED_CONFIG_TEMPLATE, indent=2)); return 0
        if a.cmd == "topology-from-csv":
            import topology_csv
            args = [str(a.branches_csv), str(a.organs_csv), str(a.out_json), "--unit", a.unit]
            if a.report:
                args += ["--report", str(a.report)]
            return int(topology_csv.main(args))
        if not a.config.exists():
            print(f"simsoma: config not found: {a.config}", file=sys.stderr); return 2
        if a.cmd == "layers":
            import layered
            res = layered.run_layered_config(a.config)
            print(f"wrote {len(res['files'])} files to {res['output_dir']}"); return 0
        import run_from_config
        if a.cmd == "check":
            run_from_config.main(["--config", str(a.config), "--step", "check"]); return 0
        if a.splits > 1:
            import launch_grid_splits
            args = ["--master-config", str(a.config), "--n-splits", str(a.splits), "--jobs", str(a.jobs),
                    "--split-axis", a.split_axis]
            if a.skip_check:
                args.append("--skip-check")
            launch_grid_splits.main(args); return 0
        if not a.skip_check:
            run_from_config.main(["--config", str(a.config), "--step", "check"])
        run_from_config.main(["--config", str(a.config), "--step", "run"]); return 0
    except (ValueError, FileNotFoundError, KeyError) as exc:   # configuration errors: short message
        print(f"simsoma: error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
