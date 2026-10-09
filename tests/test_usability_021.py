"""Usability fixes for 0.2.1: CSV CLI, plotter coordinates, missing topology file."""
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "simSOMA_corefunc"))
sys.path.insert(0, str(ROOT))
import plot_topology_json  # noqa: E402

TOPO = {"unit": "years",
        "branches": [{"id": "trunk", "parent": None, "length": 30},
                     {"id": "A", "parent": "trunk", "length": 20},
                     {"id": "A1", "parent": "A", "length": 8}],
        "events": [{"branch": "trunk", "pos": 1 / 3, "type": "branch", "target": "A"},
                   {"branch": "A", "pos": 0.6, "type": "branch", "target": "A1"},
                   {"branch": "A1", "pos": 1.0, "type": "organ", "target": "leaf"},
                   {"branch": "trunk", "time": 15.0, "type": "organ", "target": "mid"}]}


class UsabilityTests(unittest.TestCase):
    def test_fill_missing_coordinates(self):
        t = json.loads(json.dumps(TOPO))
        plot_topology_json.fill_missing_coordinates(t)
        b = {x["id"]: x for x in t["branches"]}
        self.assertAlmostEqual(b["A"]["start_age"], 10.0)
        self.assertAlmostEqual(b["A1"]["start_age"], 22.0)
        self.assertAlmostEqual(b["A1"]["end_age"], 30.0)
        ev = {e["target"]: e for e in t["events"]}
        self.assertAlmostEqual(ev["leaf"]["age"], 30.0)
        self.assertAlmostEqual(ev["mid"]["age"], 15.0)

    def test_existing_coordinates_kept(self):
        t = json.loads(json.dumps(TOPO))
        t["branches"][1]["start_age"] = 11.0
        plot_topology_json.fill_missing_coordinates(t)
        self.assertEqual(t["branches"][1]["start_age"], 11.0)

    def test_cli_topology_from_csv(self):
        d = Path(tempfile.mkdtemp())
        (d / "b.csv").write_text("branch_id,parent_id,start,end\ntrunk,,0,30\nA,trunk,10,30\n")
        (d / "o.csv").write_text("organ_id,branch_id,position\nleaf_top,trunk,30\nleaf_A,A,30\n")
        from simsoma import cli
        rc = cli.main(["topology-from-csv", str(d / "b.csv"), str(d / "o.csv"), str(d / "t.json"),
                       "--unit", "years"])
        self.assertEqual(rc, 0)
        self.assertEqual(len(json.loads((d / "t.json").read_text())["branches"]), 2)

    def test_missing_topology_file_message(self):
        d = Path(tempfile.mkdtemp())
        cfg = {"run": {"experiment_name": "e", "outdir_root": str(d / "out")},
               "topology": {"topology_json": "missing.json", "mapping_unit": "years", "mapping_rate": 5},
               "simulation": {"n_sim": 1, "modules": {
                   "self_renewal": {"m": 2, "rho": 0, "mu_unit": 1, "victim_locality": 0, "bias_mode": "fixed",
                                    "branch_bias_value": 0, "branch_bias_mean": 0, "branch_bias_kappa": 1},
                   "pre_branching": {"sam_boundary_cells": 8}, "branching": {"branch_precursor_number": 2},
                   "organ": {"organ_precursor_number": 2, "organ_total_cells": 16, "seq_fraction": 1}}}}
        p = d / "c.json"; p.write_text(json.dumps(cfg))
        import run_from_config
        with self.assertRaisesRegex(FileNotFoundError, "relative to the folder of the config file"):
            run_from_config.main(["--config", str(p), "--step", "check"])



class SplitPathTests(unittest.TestCase):
    def test_config_relative_topology_with_splits(self):
        import shutil
        d = Path(tempfile.mkdtemp())
        (d / "b.csv").write_text("branch_id,parent_id,start,end\ntrunk,,0,4\nA,trunk,2,4\n")
        (d / "o.csv").write_text("organ_id,branch_id,position\nleaf_top,trunk,4\nleaf_A,A,4\n")
        from simsoma import cli
        self.assertEqual(cli.main(["topology-from-csv", str(d / "b.csv"), str(d / "o.csv"),
                                   str(d / "t.json"), "--unit", "years"]), 0)
        cfg = {"run": {"experiment_name": "split", "outdir_root": "out", "seed": 1},
               "topology": {"topology_json": "t.json", "mapping_unit": "years", "mapping_rate": 2},
               "simulation": {"n_sim": 1, "modules": {
                   "self_renewal": {"m": {"values": [2, 3]}, "rho": 0, "mu_unit": 1, "victim_locality": 0,
                                    "bias_mode": "fixed", "branch_bias_value": 0, "branch_bias_mean": 0,
                                    "branch_bias_kappa": 1},
                   "pre_branching": {"sam_boundary_cells": 8}, "branching": {"branch_precursor_number": 2},
                   "organ": {"organ_precursor_number": 2, "organ_total_cells": 16, "seq_fraction": 1}}}}
        (d / "c.json").write_text(json.dumps(cfg))
        import os
        os.environ.setdefault("SIMSOMA_SKIP_TOPOLOGY_CONFIRM", "1")
        rc = cli.main(["run", str(d / "c.json"), "--splits", "2", "--jobs", "1"])
        self.assertEqual(rc, 0)
        self.assertTrue((d / "out" / "split" / "grid_parameter" / "parameter_sets.csv").exists())
        shutil.rmtree(d)


if __name__ == "__main__":
    unittest.main()
