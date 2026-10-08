"""Tests for the TLS segment-table converter (topology_tls)."""
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "simSOMA_corefunc"))
import topology_io, topology_tls  # noqa: E402

# trunk: 1 -> 2 -> 3 (order 0, 2 m each, base 0/2/4); branch 4 (order 1) on segment 2 at base 3.0,
# continued by 5; branch 6 (order 1) on segment 3 at base 5.0; 7 (order 2) on 5.
HEADER = "segment\tparent_segment\tbranch_order\tlength_m\tbase_distance_m\tdiameter_base_cm\n"
ROWS = [(1, 0, 0, 2.0, 0.0, 30), (2, 1, 0, 2.0, 2.0, 25), (3, 2, 0, 2.0, 4.0, 20),
        (4, 2, 1, 1.0, 3.0, 8), (5, 4, 1, 1.0, 4.0, 6), (6, 3, 1, 0.5, 5.0, 5), (7, 5, 2, 0.4, 4.5, 2)]


def write(rows, extra=""):
    d = Path(tempfile.mkdtemp())
    p = d / "seg.txt"
    p.write_text(HEADER + "".join("\t".join(map(str, r)) + "\n" for r in rows) + extra)
    return p


class TLSTests(unittest.TestCase):
    def test_axes_positions_and_load(self):
        rows = topology_tls.read_segment_table(write(ROWS))
        topo, rep = topology_tls.convert(rows, prune=False)
        self.assertEqual(rep["n_axes_total"], 4)
        br = {b["source_start_segment"]: b for b in topo["branches"]}
        self.assertAlmostEqual(br[1]["length"], 6.0); self.assertAlmostEqual(br[4]["length"], 2.0)
        ev = {e["target"]: e for e in topo["events"] if e["type"] == "branch"}
        self.assertAlmostEqual(ev[br[4]["id"]]["pos"], 3.0 / 6.0)     # base 3 on a 6 m trunk
        self.assertAlmostEqual(ev[br[6]["id"]]["pos"], 5.0 / 6.0)
        self.assertAlmostEqual(ev[br[7]["id"]]["pos"], (4.5 - 3.0) / 2.0)
        self.assertEqual(rep["n_organs"], 4)
        self.assertAlmostEqual(br[7]["start_age"], 3.0 + 1.5)              # path coordinate
        self.assertAlmostEqual(br[7]["end_age"] - br[7]["start_age"], br[7]["length"])
        p = Path(tempfile.mkdtemp()) / "t.json"; p.write_text(json.dumps(topo))
        t = topology_io.load_topology_auto(p, mapping={"unit": "meters", "rate": 10.0, "mode": "deterministic"})
        self.assertEqual(len(t["branches"]), 4)

    def test_pruning_keeps_ancestors_only(self):
        rows = topology_tls.read_segment_table(write(ROWS))
        topo, rep = topology_tls.convert(rows, organs=[7])           # tip of the order-2 axis
        self.assertEqual(rep["n_branches"], 3)                         # trunk, axis 4-5, axis 7
        self.assertEqual(rep["n_organs"], 1)
        self.assertEqual({b["source_start_segment"] for b in topo["branches"]}, {1, 4, 7})

    def test_organ_policies(self):
        rows = topology_tls.read_segment_table(write(ROWS))
        self.assertEqual(topology_tls.convert(rows, organs="min_order:1")[1]["n_organs"], 3)
        self.assertEqual(topology_tls.convert(rows, organs="orders:2")[1]["n_organs"], 1)
        a = topology_tls.convert(rows, organs="random:2", seed=3)[0]
        b = topology_tls.convert(rows, organs="random:2", seed=3)[0]
        self.assertEqual(a, b)
        self.assertEqual(topology_tls.convert(rows, min_axis_length=1.0)[1]["n_organs"], 2)
        with self.assertRaises(ValueError):
            topology_tls.convert(rows, organs=[2])                     # not an axis tip

    def test_several_same_order_children(self):
        rows = topology_tls.read_segment_table(write(ROWS + [(8, 2, 0, 0.3, 4.0, 4)]))
        topo, rep = topology_tls.convert(rows, prune=False)
        self.assertEqual(rep["n_axes_total"], 5)
        self.assertTrue(any("same-order" in w for w in rep["warnings"]))
        trunk = [b for b in topo["branches"] if b["parent"] is None][0]
        self.assertEqual(trunk["source_tip_segment"], 3)                # thicker child continues the trunk

    def test_errors(self):
        with self.assertRaises(ValueError):                             # missing parent
            topology_tls.convert(topology_tls.read_segment_table(write([(1, 0, 0, 1.0, 0.0, 1), (2, 9, 0, 1.0, 1.0, 1)])))
        with self.assertRaises(ValueError):
            topology_tls.convert(topology_tls.read_segment_table(write([(1, 0, 0, 1.0, 0.0, 1), (2, 0, 0, 1.0, 0.0, 1)])))  # two roots
        p = Path(tempfile.mkdtemp()) / "bad.txt"; p.write_text("segment\tlength_m\n1\t2\n")
        with self.assertRaises(ValueError):
            topology_tls.read_segment_table(p)


if __name__ == "__main__":
    unittest.main()
