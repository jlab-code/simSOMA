"""Regression tests for the 2026-10 code fixes (see CHANGELOG, 0.2.0-dev).

1/2  Stochastic (Poisson) topology mapping removed; deterministic mapping exact for long branches.
4    Displacement: the victim receives a second, independently mutated daughter of the displacer.
6    mu_unit is the canonical name; mu_year remains a deprecated alias (config + params).
"""
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "simSOMA_corefunc"))

import topology_io  # noqa: E402
import self_renewal  # noqa: E402
import run_from_config  # noqa: E402


def _topo(T_obs):
    return {"root_id": "B0", "branches": {"B0": {"parent_id": None, "T": T_obs,
            "events": [{"time": T_obs / 2, "type": "ORGAN", "target_id": "o1"}]}}}


class MappingTests(unittest.TestCase):
    def test_poisson_mode_rejected_with_clear_error(self):
        with self.assertRaisesRegex(ValueError, "removed"):
            topology_io.map_topology_to_sr_steps(_topo(10.0), unit="years", rate=5.0, mode="poisson")

    def test_default_is_deterministic(self):
        out = topology_io.map_topology_to_sr_steps(_topo(10.0), unit="years", rate=5.0)
        self.assertEqual(out["mapping"]["mode"], "deterministic")

    def test_long_branch_not_truncated(self):
        # the removed Poisson sampler saturated at ~745 steps
        out = topology_io.map_topology_to_sr_steps(_topo(1000.0), unit="years", rate=5.0)
        self.assertEqual(out["branches"]["B0"]["T"], 5000)
        self.assertEqual(out["branches"]["B0"]["events"][0]["time"], 2500)


class DisplacementTests(unittest.TestCase):
    def test_victim_gets_independent_daughter(self):
        params = self_renewal.SelfRenewalParams(m=2, rho=1.0, mu_div=20.0)
        sim = self_renewal.SelfRenewalSimulator(rng=__import__("numpy").random.default_rng(1))
        out = sim.simulate_segment(params, self_renewal.SelfRenewalInputs(branch_id="B0", T=1))
        a, b = (c.genotype for c in out.final_state)
        self.assertEqual(out.num_displacements, 1)
        # both cells descend from the same (mutation-free) parent; new mutations must be independent
        self.assertTrue(len(a) > 0 and len(b) > 0)
        self.assertEqual(len(a & b), 0)
        self.assertEqual(out.final_state[0].lineage_id, out.final_state[1].lineage_id)


class NamingTests(unittest.TestCase):
    def test_params_alias(self):
        p = self_renewal.SelfRenewalParams(m=3, rho=0.0, mu_year=1.5, kappa_sr=5.0)
        self.assertEqual(p.mu_unit, 1.5)
        self.assertEqual(p.mu_year, 1.5)
        self.assertAlmostEqual(self_renewal.SelfRenewalSimulator().mu_div(p), 0.3)

    def test_params_alias_conflict(self):
        with self.assertRaises(ValueError):
            self_renewal.SelfRenewalParams(m=3, rho=0.0, mu_unit=1.0, mu_year=2.0, kappa_sr=5.0)

    def test_config_alias_mapped(self):
        mods = {"self_renewal": {"mu_year": {"values": [1.0]}}, "branching": {}}
        out = run_from_config._normalize_module_param_locations(mods)
        self.assertIn("mu_unit", out["self_renewal"])
        self.assertNotIn("mu_year", out["self_renewal"])

    def test_config_alias_conflict(self):
        mods = {"self_renewal": {"mu_year": 1.0, "mu_unit": 2.0}, "branching": {}}
        with self.assertRaises(ValueError):
            run_from_config._normalize_module_param_locations(mods)


if __name__ == "__main__":
    unittest.main()
