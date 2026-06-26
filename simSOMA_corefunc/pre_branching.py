"""pre_branching.py

Pre-branching amplification to a fixed peripheral SAM boundary size.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence, Tuple, List

import numpy as np

from self_renewal import CellState, SelfRenewalParams, Time2D
from recruitment_utils import allocate_equalized_quotas, expand_founders_to_exact_leaves


@dataclass(frozen=True)
class PreBranchingInputs:
    branch_id: str
    event_age: int
    boundary_cells_target: int
    start_dev: int = 0
    seed: Optional[int] = None


@dataclass(frozen=True)
class PreBranchingOutputs:
    ring_state: Tuple[CellState, ...]
    sector_bounds: Tuple[Tuple[int, int], ...]
    founder_leaf_counts: Tuple[int, ...]
    leaf_depths: Tuple[int, ...]
    boundary_cells_realized: int
    start_time: Time2D
    end_time: Time2D



def pre_branching_amplify(
    params: SelfRenewalParams,
    source_state: Sequence[CellState],
    inputs: PreBranchingInputs,
    rng: np.random.Generator,
) -> PreBranchingOutputs:
    m = int(params.m)
    if len(source_state) != m:
        raise ValueError(f'source_state must have length m={m}')
    B = int(inputs.boundary_cells_target)
    if B < m:
        raise ValueError(f'boundary_cells_target must be >= m; got B={B}, m={m}')

    local_rng = np.random.default_rng(int(inputs.seed)) if inputs.seed is not None else rng
    quotas = allocate_equalized_quotas(B, m)

    ring_state: List[CellState] = []
    leaf_depths: List[int] = []
    bounds: List[Tuple[int, int]] = []
    cur = 0
    for founder, quota in zip(source_state, quotas):
        leaves, depths = expand_founders_to_exact_leaves(params, [founder], [int(quota)], rng=local_rng)
        start = cur
        ring_state.extend(leaves)
        leaf_depths.extend(depths)
        cur += len(leaves)
        bounds.append((start, cur))

    max_depth = max(leaf_depths) if leaf_depths else 0
    start_time: Time2D = (int(inputs.event_age), int(inputs.start_dev))
    end_time: Time2D = (int(inputs.event_age), int(inputs.start_dev) + int(max_depth))

    return PreBranchingOutputs(
        ring_state=tuple(ring_state),
        sector_bounds=tuple((int(a), int(b)) for a, b in bounds),
        founder_leaf_counts=tuple(int(q) for q in quotas),
        leaf_depths=tuple(int(d) for d in leaf_depths),
        boundary_cells_realized=int(len(ring_state)),
        start_time=start_time,
        end_time=end_time,
    )
