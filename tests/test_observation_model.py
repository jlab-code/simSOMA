#!/usr/bin/env python3
"""Focused tests for config-driven observation output and fitSOMA handoff helpers."""
from __future__ import annotations

import json
import math
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CORE = ROOT / "simSOMA_corefunc"
if str(CORE) not in sys.path:
    sys.path.insert(0, str(CORE))

from configured_observation_model import (  # noqa: E402
    _canonical_truth_aliases,
    _contract_sha256,
    _deterministic_evidence,
    _observation_contract,
    _read_count_evidence,
    _sample_info,
    normalize_observation_model_config,
)


class ObservationConfigTests(unittest.TestCase):
    def test_deterministic_defaults_are_identity(self) -> None:
        cfg = normalize_observation_model_config({})
        self.assertEqual(cfg["mode"], "deterministic")
        self.assertEqual(cfg["layers"], ["layer_equivalent"])
        self.assertEqual(cfg["layer_weights"], [1.0])
        self.assertEqual(cfg["sampling"], "layer_specific")
        self.assertEqual(cfg["target_layer"], "layer_equivalent")
        self.assertEqual(cfg["phase"], "phased")


    def test_top_level_retain_called_any_alias_is_honoured(self) -> None:
        cfg = normalize_observation_model_config({
            "mode": "read_counts",
            "layers": ["L2"],
            "layer_weights": [1.0],
            "retain_called_any": False,
        })
        self.assertFalse(cfg["read_counts"]["retain_called_any"])
        contract = _observation_contract(cfg)
        self.assertFalse(contract["layers"][0]["ascertainment"]["retain_called_any"])
        self.assertEqual(len(_contract_sha256(contract)), 64)

    def test_nested_retain_called_any_precedence_is_backward_compatible(self) -> None:
        cfg = normalize_observation_model_config({
            "mode": "read_counts",
            "layers": ["L2"],
            "layer_weights": [1.0],
            "retain_called_any": False,
            "ascertainment": {"retain_called_any": False},
            "read_counts": {"retain_called_any": True},
        })
        self.assertTrue(cfg["read_counts"]["retain_called_any"])

    def test_invalid_weights_are_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "sum to 1.0"):
            normalize_observation_model_config({
                "layers": ["L1", "L2"],
                "layer_weights": [0.2, 0.2],
            })

    def test_deterministic_identity_preserves_exact_vaf_and_matrix(self) -> None:
        raw = pd.DataFrame([
            {"mutation_id": "m1", "organ_id": "O1", "allele_count": 1, "sequenced_cells": 4},
            {"mutation_id": "m1", "organ_id": "O2", "allele_count": 2, "sequenced_cells": 4},
            {"mutation_id": "m2", "organ_id": "O2", "allele_count": 1, "sequenced_cells": 2},
        ])
        cfg = normalize_observation_model_config({
            "mode": "deterministic",
            "layers": ["layer_equivalent"],
            "layer_weights": [1.0],
            "sampling": "layer_specific",
            "target_layer": "layer_equivalent",
            "phase": "phased",
            "deterministic_depth": 1000,
        })
        info = _sample_info(["O1", "O2"], cfg)
        ev = _deterministic_evidence(raw, info, cfg)
        self.assertEqual(len(ev), 4)
        self.assertEqual(ev["mutation_id"].nunique(), 2)
        self.assertTrue((ev.groupby("mutation_id")["sample_id"].nunique() == 2).all())
        lookup = ev.set_index(["mutation_id", "sample_id"])["exact_vaf"].to_dict()
        self.assertEqual(lookup[("layer_equivalent::m1", "O1")], 0.25)
        self.assertEqual(lookup[("layer_equivalent::m1", "O2")], 0.50)
        self.assertEqual(lookup[("layer_equivalent::m2", "O1")], 0.00)
        self.assertEqual(lookup[("layer_equivalent::m2", "O2")], 0.50)

    def test_deterministic_bulk_unphased_scaling(self) -> None:
        raw = pd.DataFrame([
            {"mutation_id": "m1", "organ_id": "O1", "allele_count": 4, "sequenced_cells": 4},
        ])
        cfg = normalize_observation_model_config({
            "mode": "deterministic",
            "layers": ["L1", "L2", "L3"],
            "layer_weights": [0.1, 0.7, 0.2],
            "sampling": "bulk",
            "phase": "unphased",
            "deterministic_depth": 10000,
        })
        ev = _deterministic_evidence(raw, _sample_info(["O1"], cfg), cfg)
        by_layer = ev.set_index("source_layer")["exact_vaf"].to_dict()
        self.assertAlmostEqual(by_layer["L1"], 0.05)
        self.assertAlmostEqual(by_layer["L2"], 0.35)
        self.assertAlmostEqual(by_layer["L3"], 0.10)

    def test_read_counts_are_reproducible_and_valid(self) -> None:
        raw = pd.DataFrame([
            {"mutation_id": "m1", "organ_id": "O1", "allele_count": 1, "sequenced_cells": 4},
            {"mutation_id": "m1", "organ_id": "O2", "allele_count": 2, "sequenced_cells": 4},
            {"mutation_id": "m2", "organ_id": "O2", "allele_count": 1, "sequenced_cells": 2},
        ])
        cfg = normalize_observation_model_config({
            "mode": "read_counts",
            "layers": ["L2"],
            "layer_weights": [1.0],
            "sampling": "layer_specific",
            "target_layer": "L2",
            "phase": "phased",
            "depth": {"mode": "fixed", "value": 100},
            "read_sampling": {"distribution": "beta_binomial", "concentration": 200},
            "sequencing_error": {"reference_to_alternate": 0.001, "alternate_to_reference": 0.001},
            "caller": {"minimum_depth": 20, "minimum_alt_reads": 3, "minimum_vaf": 0.02},
        })
        info = _sample_info(["O1", "O2"], cfg)
        a, _ = _read_count_evidence(raw, info, cfg, 12345)
        b, _ = _read_count_evidence(raw, info, cfg, 12345)
        pd.testing.assert_frame_equal(a.reset_index(drop=True), b.reset_index(drop=True))
        self.assertTrue((a["alt_count"] >= 0).all())
        self.assertTrue((a["alt_count"] <= a["depth"]).all())
        self.assertTrue((a["depth"] == 100).all())
        self.assertTrue((a.groupby("mutation_id")["sample_id"].nunique() == 2).all())

    def test_read_counts_preserve_assay_scaling_before_noise(self) -> None:
        raw = pd.DataFrame([
            {"mutation_id": "m1", "organ_id": "O1", "allele_count": 4, "sequenced_cells": 4},
        ])
        cfg = normalize_observation_model_config({
            "mode": "read_counts",
            "layers": ["L1", "L2", "L3"],
            "layer_weights": [0.1, 0.7, 0.2],
            "sampling": "bulk",
            "phase": "unphased",
            "depth": {"mode": "fixed", "value": 100},
            "read_sampling": {"distribution": "binomial", "concentration": 200},
            "sequencing_error": {"reference_to_alternate": 0.0, "alternate_to_reference": 0.0},
            "caller": {"minimum_depth": 1, "minimum_alt_reads": 0, "minimum_vaf": 0.0},
            "read_counts": {"retain_called_any": False},
        })
        ev, _ = _read_count_evidence(raw, _sample_info(["O1"], cfg), cfg, 456)
        by_layer = ev.set_index("source_layer")["assay_vaf"].to_dict()
        self.assertAlmostEqual(by_layer["L1"], 0.05)
        self.assertAlmostEqual(by_layer["L2"], 0.35)
        self.assertAlmostEqual(by_layer["L3"], 0.10)

    def test_asymmetric_error_is_rejected(self) -> None:
        raw = pd.DataFrame([
            {"mutation_id": "m1", "organ_id": "O1", "allele_count": 1, "sequenced_cells": 2},
        ])
        cfg = normalize_observation_model_config({
            "mode": "read_counts",
            "layers": ["L2"],
            "layer_weights": [1.0],
            "sampling": "layer_specific",
            "target_layer": "L2",
            "phase": "phased",
            "sequencing_error": {"reference_to_alternate": 0.001, "alternate_to_reference": 0.002},
        })
        with self.assertRaisesRegex(ValueError, "symmetric sequencing-error"):
            _read_count_evidence(raw, _sample_info(["O1"], cfg), cfg, 789)

    def test_truth_aliases_use_correct_turnover_mapping(self) -> None:
        configured = {"mu_year": 0.3, "effective_kappa_sr": 5.0, "rho": 0.2, "m": 4, "sam_boundary_cells": 12}
        aliases = _canonical_truth_aliases(configured, {"topology": {"phyllotaxy": {"mode": "off"}}}, {"mode": "deterministic"})
        expected_lambda = -5.0 * math.log(0.8)
        self.assertAlmostEqual(aliases["lambda_turn"], expected_lambda)
        self.assertAlmostEqual(aliases["rho_turn"], 0.2)
        self.assertAlmostEqual(aliases["p_turn"], 1.0 - math.exp(-expected_lambda))
        self.assertEqual(aliases["mu_unit"], 0.3)
        self.assertEqual(aliases["K_sr"], 5.0)


if __name__ == "__main__":
    unittest.main()
