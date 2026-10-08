"""Founder clonal composition, definition local_sector_v1 (see founder_diversity.py).

Checks: closed form; agreement of realized polyclonal fraction with the exact expectation;
independence of turnover rho and branch order (the old root-label statistic failed both).
"""
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "simSOMA_corefunc"))

import founder_diversity as fd  # noqa: E402
import pipeline_wrapper  # noqa: E402
import self_renewal  # noqa: E402
import topology_io  # noqa: E402

TOPO = ROOT / "simSOMA_inputs" / "examples" / "example_03_nested_branches.json"


def _events(rho, P_b, P_o, reps, m=4, C=64):
    topo = topology_io.load_topology_auto(TOPO, mapping={"unit": "years", "rate": 5.0, "mode": "deterministic"})
    order = {}
    for b, spec in topo["branches"].items():
        d, x = 0, b
        while topo["branches"][x]["parent_id"] is not None:
            x = topo["branches"][x]["parent_id"]; d += 1
        order[b] = d
    br, org = [], []
    for s in range(reps):
        sr = self_renewal.SelfRenewalParams(m=m, rho=rho, mu_unit=1.0, kappa_sr=5.0)
        r = pipeline_wrapper.run_pipeline(topo, sr, sam_boundary_cells=C, branch_precursor_number=P_b,
                                          organ_precursor_number=P_o, organ_total_cells=64, sequenced_cells=64, seed=s)
        br += [(order[e["child_id"]], e) for e in r["branch_events"]]
        org += [e for e in r["organ_events"]]
    return br, org


class ClosedFormTests(unittest.TestCase):
    def test_closed_form_divisible(self):
        for m, C in [(2, 16), (4, 64), (8, 64), (3, 48)]:
            for P in range(1, C // m + 2):
                e = fd.polyclonal_founding_expectation(P, C, m)
                self.assertAlmostEqual(e["probability"], min(1.0, (P - 1) * m / C))
                self.assertAlmostEqual(e["expected_sector_count"], 1.0 + e["probability"])
                self.assertAlmostEqual(e["phi"], (P - 1) * m / C)

    def test_single_precursor_is_monoclonal(self):
        self.assertEqual(fd.polyclonal_founding_expectation(1, 50, 7)["probability"], 0.0)

    def test_summary_fields(self):
        s = fd.summarize_founder_sectors([15, 16], 64, 4)   # straddles the sector boundary at 16
        self.assertEqual(s["founder_sector_count"], 2)
        self.assertEqual(s["founder_polyclonal"], 1)
        self.assertAlmostEqual(s["founder_effective_sectors"], 2.0)
        self.assertEqual(s["founder_definition"], fd.FOUNDER_DEFINITION_VERSION)


class SimulationTests(unittest.TestCase):
    def _check(self, values, expected, n):
        se = (expected * (1 - expected) / n) ** 0.5
        self.assertLess(abs(sum(values) / n - expected), 4 * se + 1e-9)

    def test_organs_match_expectation_and_ignore_turnover(self):
        exp = fd.polyclonal_founding_expectation(8, 64, 4)["probability"]   # 7/16
        for rho in (0.0, 0.5):
            _, org = _events(rho, P_b=4, P_o=8, reps=150)
            self._check([e["founder_polyclonal"] for e in org], exp, len(org))
        # the deprecated root-label statistic collapses under turnover
        _, org = _events(0.5, P_b=4, P_o=8, reps=30)
        self.assertEqual(sum(e["root_lineage_sector_count"] > 1 for e in org), 0)

    def test_branches_match_expectation_at_every_order(self):
        exp = fd.polyclonal_founding_expectation(4, 64, 4)["probability"]   # 3/16
        br, _ = _events(0.2, P_b=4, P_o=8, reps=300)
        for k in sorted({o for o, _ in br}):
            vals = [e["founder_polyclonal"] for o, e in br if o == k]
            self._check(vals, exp, len(vals))


if __name__ == "__main__":
    unittest.main()
