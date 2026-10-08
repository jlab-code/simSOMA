"""Shared observation model (plantsoma_obs), per-layer runs (layered.py), organ mutation
multiplier, and CLI smoke tests."""
import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "simSOMA_corefunc"))

import plantsoma_obs  # noqa: E402
import layered  # noqa: E402
import pipeline_wrapper, self_renewal, topology_io  # noqa: E402,E401

TOPO = ROOT / "simSOMA_inputs" / "examples" / "example_03_nested_branches.json"
MAP = {"unit": "years", "rate": 5.0, "mode": "deterministic"}


class ObservationModelTests(unittest.TestCase):
    def test_validation_errors_name_the_key(self):
        for bad, key in [({"depth": {"mode": "magic"}}, "depth.mode"), ({"reads": {"type": "x"}}, "reads.type"),
                         ({"background": {"distribution": "weird"}}, "background.distribution"),
                         ({"depht": {}}, "sections"), ({"depth": {"meen": 3}}, "depth")]:
            with self.assertRaisesRegex(ValueError, key):
                plantsoma_obs.normalize_config(bad)

    def test_hash_stable_and_sensitive(self):
        a = plantsoma_obs.config_sha256({"depth": {"mean": 60}})
        self.assertEqual(a, plantsoma_obs.config_sha256({"depth": {"mean": 60.0, "mode": "normal"}}))
        self.assertNotEqual(a, plantsoma_obs.config_sha256({"depth": {"mean": 61}}))

    def test_normal_depth_matches_direct_draws(self):
        cfg = plantsoma_obs.normalize_config({"depth": {"mode": "normal", "mean": 80, "sd": 10}})
        d = plantsoma_obs.draw_depths(50, 3, cfg["depth"], np.random.default_rng(3))
        rng = np.random.default_rng(3)
        ref = np.column_stack([np.maximum(np.rint(rng.normal(80, 10, 50)).astype(int), 1) for _ in range(3)])
        self.assertTrue((d == ref).all())

    def test_lognormal_depth_v4(self):
        cfg = plantsoma_obs.normalize_config({"depth": {"mode": "lognormal_site_sample", "mean": 60}})
        self.assertEqual((cfg["depth"]["site_sdlog"], cfg["depth"]["sample_sdlog"], cfg["depth"]["sample_factor_sdlog"]),
                         (0.58, 0.17, 0.27))
        d = plantsoma_obs.draw_depths(20000, 14, cfg["depth"], np.random.default_rng(1))
        # site factor truncated at 3: E[g | g <= 3] = Phi(z - s) / Phi(z), z = (ln 3 + s^2/2) / s
        from statistics import NormalDist
        s = 0.58; z = (np.log(3) + s * s / 2) / s
        # sample factor h (sdlog 0.27, 14 samples) adds noise to the overall mean
        self.assertAlmostEqual(d.mean() / 60, NormalDist().cdf(z - s) / NormalDist().cdf(z), delta=0.12)
        # organ factor: log of per-sample mean depths has sd ~ 0.27
        self.assertAlmostEqual(np.log(d.mean(axis=0)).std(ddof=1), 0.27, delta=0.09)
        # matches the v4 generator draw order g, e, h exactly
        rng = np.random.default_rng(7); S, O = 500, 6
        g = np.exp(rng.normal(-0.58 ** 2 / 2, 0.58, 4 * S)); g = g[g <= 3][:S]
        e = np.exp(rng.normal(-0.17 ** 2 / 2, 0.17, (S, O)))
        h = np.exp(rng.normal(-0.27 ** 2 / 2, 0.27, O))
        ref = rng.poisson(60 * g[:, None] * e * h[None, :])
        cfg0 = plantsoma_obs.normalize_config({"depth": {"mode": "lognormal_site_sample", "mean": 60, "minimum": 0}})
        self.assertTrue((plantsoma_obs.draw_depths(S, O, cfg0["depth"], np.random.default_rng(7)) == ref).all())

    def test_background_and_ascertainment(self):
        cfg = {"background": {"n_sites": 500, "distribution": "gamma:2", "mean": 0.023},
               "reads": {"type": "binomial", "sequencing_error": 0.0},
               "caller": {"min_depth": 0, "min_alt_reads": 2}}
        obs = plantsoma_obs.observe(np.zeros((0, 5)), cfg, np.random.default_rng(0))
        self.assertEqual(obs["is_background"].sum(), 500)
        self.assertLessEqual(obs["assay_vaf"].max(), 0.45)
        self.assertTrue(np.all(obs["assay_vaf"] == obs["assay_vaf"][:, :1]))   # same q in all samples
        self.assertTrue(0 < obs["keep"].mean() < 1)
        self.assertEqual(obs["version"], plantsoma_obs.OBSERVATION_MODEL_VERSION)

    def test_fixed_vaf_read_expectation(self):
        cfg = {"depth": {"mode": "fixed", "mean": 200}, "reads": {"type": "binomial", "sequencing_error": 0.0},
               "caller": {"retain_called_any": False}}
        obs = plantsoma_obs.observe(np.full((4000, 1), 0.25), cfg, np.random.default_rng(0))
        self.assertAlmostEqual(obs["observed_vaf"].mean(), 0.25, delta=0.003)


class LayeredTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.topo = topology_io.load_topology_auto(TOPO, mapping=MAP)
        cls.car = layered.simulate_layers(
            cls.topo, {"L1": {"mu_unit": 2.0}, "L2": {"mu_unit": 1.0}}, m=3, rho=0.0, kappa_sr=5.0,
            sam_boundary_cells=30, branch_precursor_number=2, organ_precursor_number=6, organ_total_cells=200,
            sequenced_cells=200, seed=1)

    def test_layers_independent_and_rate_scaled(self):
        ids = {k: set(g.mutation_id) for k, g in self.car.groupby("layer")}
        self.assertFalse(ids["L1"] & ids["L2"])
        n = {"L1": 0, "L2": 0}
        for rep in range(8):   # single trees vary strongly (clonal fixation); pool replicates
            c = layered.simulate_layers(self.topo, {"L1": {"mu_unit": 2.0}, "L2": {"mu_unit": 1.0}}, m=3, rho=0.0,
                                        kappa_sr=5.0, sam_boundary_cells=30, branch_precursor_number=2,
                                        organ_precursor_number=6, organ_total_cells=200, sequenced_cells=200,
                                        seed=1, replicate=rep)
            for k in n:
                n[k] += c[c.layer == k].mutation_id.nunique()
        self.assertAlmostEqual(n["L1"] / n["L2"], 2.0, delta=0.25)

    def test_assay_vaf_formula(self):
        organs = sorted(self.car.organ_id.unique())
        cm = layered.resolve_contributions({"L1": 0.3, "L2": 0.7}, organs, ["L1", "L2"])
        sites, organs, V, F = layered.assay_vaf_matrix(self.car, cm, "unphased")
        c = np.where(sites.layer.to_numpy() == "L1", 0.3, 0.7)[:, None]
        self.assertTrue(np.allclose(V, 0.5 * c * F))

    def test_contribution_validation(self):
        with self.assertRaisesRegex(ValueError, "sum to 1"):
            layered.resolve_contributions({"L1": 0.3, "L2": 0.3}, ["o"], ["L1", "L2"])
        with self.assertRaisesRegex(ValueError, "unknown layers"):
            layered.resolve_contributions({"L1": 0.5, "LX": 0.5}, ["o"], ["L1", "L2"])

    def test_per_organ_contributions(self):
        organs = sorted(self.car.organ_id.unique())
        spec = {"default": {"L1": 0.5, "L2": 0.5}, organs[0]: {"L1": 1.0}}
        cm = layered.resolve_contributions(spec, organs, ["L1", "L2"])
        self.assertEqual(cm.loc[organs[0], "L1"], 1.0)
        self.assertEqual(cm.loc[organs[-1], "L2"], 0.5)


class OrganMultiplierTests(unittest.TestCase):
    def test_multiplier_scales_private_low_vaf(self):
        topo = topology_io.load_topology_auto(TOPO, mapping=MAP)
        low = {}
        for mult in (1.0, 4.0):
            n = 0
            for s in range(6):
                sr = self_renewal.SelfRenewalParams(m=3, rho=0.0, mu_unit=1.0, kappa_sr=5.0)
                r = pipeline_wrapper.run_pipeline(topo, sr, sam_boundary_cells=30, branch_precursor_number=2,
                                                  organ_precursor_number=4, organ_total_cells=256, sequenced_cells=256,
                                                  seed=s, organ_mu_multiplier=mult)
                n += sum(1 for oe in r["organ_events"] for a in oe["allele_counts_by_mutation"].values() if a / 256 < 0.05)
            low[mult] = n
        self.assertAlmostEqual(low[4.0] / low[1.0], 4.0, delta=0.6)


class CliTests(unittest.TestCase):
    def test_version_and_layers(self):
        from simsoma.cli import main
        self.assertEqual(main(["version"]), 0)
        with tempfile.TemporaryDirectory() as d:
            cfg = json.loads(json.dumps(layered.LAYERED_CONFIG_TEMPLATE))
            cfg["topology"]["topology_json"] = str(TOPO)
            cfg["simulation"]["parameters"].update(organ_total_cells=200, sequenced_cells=200)
            cfg["observation"]["model"]["background"]["n_sites"] = 100
            cfg["output"]["dir"] = str(Path(d) / "out")
            p = Path(d) / "c.json"; p.write_text(json.dumps(cfg))
            self.assertEqual(main(["layers", str(p)]), 0)
            out = Path(d) / "out" / "replicate_0000"
            self.assertTrue((out / "layer_carriers.csv.gz").exists() and (out / "read_evidence.csv.gz").exists())
            self.assertTrue(list(out.glob("vafsoma_dp_*_vaf.csv")))
            self.assertEqual(main(["run", str(Path(d) / "missing.json")]), 2)


if __name__ == "__main__":
    unittest.main()
