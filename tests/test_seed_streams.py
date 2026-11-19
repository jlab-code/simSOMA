"""Regression test for the seed collisions of v0.2.3 (pre-branching amplification re-seeded with 31-bit integers).

Up to v0.2.3, every pre-branching amplification was re-seeded with ``random.randrange(2**31)``, so two events could
receive the same seed and draw identical mutation identifiers in unrelated lineages (probability about N**2 / 2**32
per genealogy for N branch and organ events). From v0.3.0 each amplification draws from its event's own stream.
"""
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

from simsoma import pipeline_wrapper, pre_branching, self_renewal, topology_io

ROOT = Path(__file__).resolve().parents[1]
TOPO = ROOT / "simSOMA_inputs" / "examples" / "example_03_nested_branches.json"


def _run(seed=1, patch_randrange=False):
    topo = topology_io.load_topology_auto(TOPO, mapping={"unit": "years", "rate": 5.0, "mode": "deterministic"})
    sr = self_renewal.SelfRenewalParams(m=4, rho=0.0, mu_unit=1.0, kappa_sr=5.0)
    calls = []
    real = pre_branching.pre_branching_amplify

    def spy(params, source_state, inputs, rng):
        out = real(params, source_state, inputs, rng)
        new_ids = set()
        for cell in out.ring_state:
            new_ids.update(cell.genotype)
        for cell in source_state:
            new_ids.difference_update(cell.genotype)
        calls.append({"seed": inputs.seed, "new_ids": new_ids})
        return out

    patches = [mock.patch.object(pre_branching, "pre_branching_amplify", spy)]
    if patch_randrange:   # the v0.2.3 failure mode: every event receives the same 31-bit seed
        patches.append(mock.patch("random.Random.randrange", lambda self, *a, **k: 7))
    for p in patches:
        p.start()
    try:
        res = pipeline_wrapper.run_pipeline(topo, sr, sam_boundary_cells=32, branch_precursor_number=2,
                                            organ_precursor_number=4, organ_total_cells=32, sequenced_cells=32,
                                            seed=seed)
    finally:
        for p in patches:
            p.stop()
    return res, calls


class SeedStreamTests(unittest.TestCase):
    def test_amplification_uses_event_stream(self):
        _, calls = _run()
        self.assertGreater(len(calls), 3)
        self.assertTrue(all(c["seed"] is None for c in calls))

    def test_no_shared_new_mutations_between_events(self):
        # even if every 31-bit seed were identical, the amplifications must draw distinct identifiers
        _, calls = _run(patch_randrange=True)
        seen = set()
        for c in calls:
            self.assertFalse(seen & c["new_ids"], "two pre-branching amplifications drew the same mutation id")
            seen |= c["new_ids"]

    def test_reproducible_for_a_given_seed(self):
        a, _ = _run(seed=11)
        b, _ = _run(seed=11)
        ka = [sorted(o["allele_counts_by_mutation"].items()) for o in a["organ_events"]]
        kb = [sorted(o["allele_counts_by_mutation"].items()) for o in b["organ_events"]]
        self.assertEqual(ka, kb)

    def test_explicit_seed_still_supported(self):
        # direct callers (and plug-ins) may still pass a seed; it then defines the stream
        sr = self_renewal.SelfRenewalParams(m=2, rho=0.0, mu_unit=1.0, kappa_sr=1.0)
        state = [self_renewal.CellState(genotype=frozenset(), lineage_id=i, comp_label=0) for i in range(2)]
        inp = pre_branching.PreBranchingInputs(branch_id="b", event_age=0, boundary_cells_target=8, seed=5)
        r1 = pre_branching.pre_branching_amplify(sr, state, inp, rng=np.random.default_rng(1))
        r2 = pre_branching.pre_branching_amplify(sr, state, inp, rng=np.random.default_rng(2))
        self.assertEqual([c.genotype for c in r1.ring_state], [c.genotype for c in r2.ring_state])


if __name__ == "__main__":
    unittest.main()
