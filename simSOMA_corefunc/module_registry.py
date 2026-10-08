"""module_registry.py

Default adapters around built-in simulator modules.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple
import random
import numpy as np

import self_renewal as _sr
import pre_branching as _pre
import branching as _br
import organ as _org
from founder_diversity import FOUNDER_EVENT_FIELDS, summarize_founder_lineages, summarize_founder_sectors


class DefaultSelfRenewalAdapter:
    def initialize_founders(self, sr_params: _sr.SelfRenewalParams, *, rng: np.random.Generator) -> List[_sr.CellState]:
        sim = _sr.SelfRenewalSimulator(rng=rng)
        return list(sim.initialize_founders(sr_params))

    def simulate_segment(
        self,
        sr_params: _sr.SelfRenewalParams,
        *,
        branch_id: str,
        T: int,
        age_offset: int,
        init_state: Sequence[_sr.CellState],
        rng: np.random.Generator,
        snapshot_times: Optional[List[int]] = None,
    ) -> List[_sr.CellState]:
        sim = _sr.SelfRenewalSimulator(rng=rng)
        out = sim.simulate_segment(
            sr_params,
            _sr.SelfRenewalInputs(
                branch_id=str(branch_id),
                T=int(T),
                age_offset=int(age_offset),
                init_state=list(init_state),
                snapshot_times=snapshot_times,
            ),
        )
        return list(out.final_state)


class DefaultPreBranchingAdapter:
    def amplify(
        self,
        sr_params: _sr.SelfRenewalParams,
        niche_state: Sequence[_sr.CellState],
        *,
        branch_id: str,
        event_age: int,
        boundary_cells_target: int,
        start_dev: int,
        seed: int,
        rng: np.random.Generator,
    ) -> Tuple[List[_sr.CellState], int]:
        pre_in = _pre.PreBranchingInputs(
            branch_id=str(branch_id),
            event_age=int(event_age),
            boundary_cells_target=int(boundary_cells_target),
            start_dev=int(start_dev),
            seed=int(seed),
        )
        out = _pre.pre_branching_amplify(sr_params, list(niche_state), pre_in, rng=rng)
        return list(out.ring_state), int(out.end_time[1])


class DefaultBranchingAdapter:
    def sample_child_niche(
        self,
        ring_state: Sequence[_sr.CellState],
        *,
        branch_id: str,
        child_id: str,
        event_time: Tuple[int, int],
        precursor_number: int,
        niche_size_m: int,
        rng: random.Random,
        focal_index: Optional[int] = None,
    ) -> Dict[str, Any]:
        sel = _br.sample_branch_niche(
            ring_state=list(ring_state),
            branch_id=str(branch_id),
            child_id=str(child_id),
            event_time=(int(event_time[0]), int(event_time[1])),
            precursor_number=int(precursor_number),
            niche_size_m=int(niche_size_m),
            rng=rng,
            focal_index=(int(focal_index) if focal_index is not None else None),
        )
        founder = summarize_founder_sectors(sel.sampled_indices, len(ring_state), int(niche_size_m))
        founder.update(summarize_founder_lineages(
            [ring_state[i] for i in sel.sampled_indices], int(niche_size_m)))
        return {
            'branch_id': sel.branch_id,
            'child_id': sel.child_id,
            'event_time': (int(sel.event_time[0]), int(sel.event_time[1])),
            'niche_size_m': int(sel.niche_size_m),
            'precursor_number': int(sel.precursor_number),
            'requested_precursor_number': int(sel.requested_precursor_number),
            'precursor_number_capped': bool(sel.precursor_number_capped),
            'focal_index': int(sel.focal_index),
            'sampled_indices': list(sel.sampled_indices),
            'precursor_quotas': list(sel.precursor_quotas),
            'expanded_state': list(sel.expanded_state),
            'precursor_wraparound': bool(sel.precursor_wraparound),
            'precursor_min_index': sel.precursor_min_index,
            'precursor_max_index': sel.precursor_max_index,
            **founder,
        }


class DefaultOrganAdapter:
    def form_organ(
        self,
        sr_params: _sr.SelfRenewalParams,
        ring_state: Sequence[_sr.CellState],
        *,
        organ_id: str,
        event_age: int,
        start_dev: int,
        precursor_number: int,
        organ_total_cells: int,
        sequenced_cells: Optional[int],
        rng: random.Random,
        np_rng: np.random.Generator,
        focal_index: Optional[int] = None,
        mu_multiplier: float = 1.0,
    ) -> Dict[str, Any]:
        inputs = _org.OrganInputs(
            organ_id=str(organ_id),
            event_age=int(event_age),
            start_dev=int(start_dev),
            precursor_number=int(precursor_number),
            organ_total_cells=int(organ_total_cells),
            sequenced_cells=sequenced_cells,
            focal_index=(int(focal_index) if focal_index is not None else None),
            mu_multiplier=float(mu_multiplier),
        )
        org_ev = _org.run_organ_event(sr_params, list(ring_state), inputs, rng=rng, np_rng=np_rng)
        return {
            'organ_id': org_ev.organ_id,
            'event_time': (int(org_ev.event_time[0]), int(org_ev.event_time[1])),
            'precursor_number': int(org_ev.precursor_number),
            'precursor_number_requested': int(org_ev.precursor_number_requested),
            'precursor_number_capped': bool(org_ev.precursor_number_capped),
            'boundary_cells': int(org_ev.boundary_cells),
            'organ_total_cells': int(org_ev.organ_total_cells),
            'focal_index': int(org_ev.focal_index),
            'precursor_indices': list(org_ev.precursor_indices),
            'precursor_quotas': list(org_ev.precursor_quotas),
            'precursor_mutations': list(getattr(org_ev, 'precursor_mutations', ())),
            'sequenced_cells': int(org_ev.sequenced_cells),
            'allele_counts_by_mutation': {int(k): int(v) for k, v in org_ev.allele_counts_by_mutation.items()},
            'lineage_counts': {int(k): int(v) for k, v in org_ev.lineage_counts.items()},
            'dominant_lineage_id': int(org_ev.dominant_lineage_id),
            'dominant_fraction': float(org_ev.dominant_fraction),
            'lineage_entropy': float(org_ev.lineage_entropy),
            'precursor_wraparound': bool(org_ev.precursor_wraparound),
            'precursor_min_index': org_ev.precursor_min_index,
            'precursor_max_index': org_ev.precursor_max_index,
            'max_dev_depth': int(org_ev.max_dev_depth),
            **{k: getattr(org_ev, k) for k in FOUNDER_EVENT_FIELDS},
        }


def default_registry() -> "plugin_api.ModuleRegistry":
    from plugin_api import ModuleRegistry
    return ModuleRegistry(
        self_renewal=DefaultSelfRenewalAdapter(),
        pre_branching=DefaultPreBranchingAdapter(),
        branching=DefaultBranchingAdapter(),
        organ=DefaultOrganAdapter(),
    )
