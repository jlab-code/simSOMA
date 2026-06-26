"""organ.py

Organ formation from a fixed-size boundary ring into a fixed-size mature organ.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import random
import numpy as np

from self_renewal import CellState, SelfRenewalParams, Time2D
from recruitment_utils import allocate_equalized_quotas, expand_founders_to_exact_leaves, select_contiguous_precursors


@dataclass(frozen=True)
class OrganInputs:
    organ_id: str
    event_age: int
    start_dev: int
    precursor_number: int
    organ_total_cells: int
    sequenced_cells: Optional[int] = None
    focal_index: Optional[int] = None


@dataclass(frozen=True)
class OrganEvent:
    organ_id: str
    event_time: Time2D
    precursor_number: int
    precursor_number_requested: int
    precursor_number_capped: bool
    boundary_cells: int
    organ_total_cells: int
    focal_index: int
    precursor_indices: Tuple[int, ...]
    precursor_quotas: Tuple[int, ...]
    precursor_mutations: Tuple[int, ...]
    sequenced_cells: int
    allele_counts_by_mutation: Dict[int, int]
    lineage_counts: Dict[int, int]
    dominant_lineage_id: int
    dominant_fraction: float
    lineage_entropy: float
    precursor_wraparound: bool
    precursor_min_index: int | None
    precursor_max_index: int | None
    max_dev_depth: int



def _compute_allele_counts_from_genotypes(genotypes: Sequence[frozenset[int]]) -> Dict[int, int]:
    mut_counts: Dict[int, int] = {}
    for g in genotypes:
        for mid in g:
            mut_counts[int(mid)] = mut_counts.get(int(mid), 0) + 1
    return mut_counts



def run_organ_event(
    params: SelfRenewalParams,
    ring_state: Sequence[CellState],
    inputs: OrganInputs,
    *,
    rng: random.Random,
    np_rng: Optional[np.random.Generator] = None,
) -> OrganEvent:
    N = len(ring_state)
    if N <= 0:
        raise ValueError('ring_state must be non-empty')
    O = int(inputs.organ_total_cells)
    if O < 1:
        raise ValueError('organ_total_cells must be >= 1')
    recruitment = select_contiguous_precursors(
        ring_state,
        requested_cells=int(inputs.precursor_number),
        max_realized_cells=int(O),
        rng=rng,
        focal_index=inputs.focal_index,
    )
    precursors = list(recruitment.sampled_state)

    precursor_mut = set()
    for c in precursors:
        precursor_mut.update(c.genotype)
    precursor_mut_tuple = tuple(sorted(int(x) for x in precursor_mut))

    quotas = allocate_equalized_quotas(O, recruitment.realized_cells)
    leaves, leaf_depths = expand_founders_to_exact_leaves(
        params,
        precursors,
        quotas,
        rng=(np_rng if np_rng is not None else np.random.default_rng(rng.randrange(2**31))),
    )

    if inputs.sequenced_cells is None:
        sampled = leaves
    else:
        n_req = int(inputs.sequenced_cells)
        if n_req <= 0:
            raise ValueError('sequenced_cells must be positive')
        if n_req >= len(leaves):
            sampled = leaves
        else:
            sampled = rng.sample(leaves, k=n_req)

    lineage_counts: Dict[int, int] = {}
    for c in sampled:
        lid = int(getattr(c, 'lineage_id'))
        lineage_counts[lid] = lineage_counts.get(lid, 0) + 1
    if len(lineage_counts) == 0:
        dominant_lineage_id = -1
        dominant_fraction = 0.0
        lineage_entropy = 0.0
    else:
        dominant_lineage_id = max(lineage_counts.items(), key=lambda kv: kv[1])[0]
        n_tot = float(sum(lineage_counts.values()))
        dominant_fraction = float(lineage_counts[dominant_lineage_id]) / n_tot if n_tot > 0 else 0.0
        lineage_entropy = 0.0
        for cnt in lineage_counts.values():
            pp = float(cnt) / n_tot
            if pp > 0:
                lineage_entropy -= pp * float(np.log(pp))

    allele_counts = _compute_allele_counts_from_genotypes([c.genotype for c in sampled])
    max_depth = max(leaf_depths) if leaf_depths else 0

    return OrganEvent(
        organ_id=str(inputs.organ_id),
        event_time=(int(inputs.event_age), int(inputs.start_dev) + int(max_depth)),
        precursor_number=int(recruitment.realized_cells),
        precursor_number_requested=int(recruitment.requested_cells),
        precursor_number_capped=bool(recruitment.capped),
        boundary_cells=int(N),
        organ_total_cells=int(len(leaves)),
        focal_index=int(recruitment.focal_index),
        precursor_indices=tuple(int(i) for i in recruitment.sampled_indices),
        precursor_quotas=tuple(int(q) for q in quotas),
        precursor_mutations=precursor_mut_tuple,
        sequenced_cells=int(len(sampled)),
        allele_counts_by_mutation={int(k): int(v) for k, v in allele_counts.items()},
        lineage_counts={int(k): int(v) for k, v in lineage_counts.items()},
        dominant_lineage_id=int(dominant_lineage_id),
        dominant_fraction=float(dominant_fraction),
        lineage_entropy=float(lineage_entropy),
        precursor_wraparound=bool(recruitment.wraparound),
        precursor_min_index=recruitment.min_index,
        precursor_max_index=recruitment.max_index,
        max_dev_depth=int(max_depth),
    )
