"""branching.py

Branch founding from a fixed-size boundary ring into a child niche of fixed size m.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple
import random

from self_renewal import CellState, Time2D
from recruitment_utils import allocate_equalized_quotas, select_contiguous_precursors


@dataclass(frozen=True)
class BranchSelection:
    branch_id: str
    child_id: str
    event_time: Time2D
    niche_size_m: int
    precursor_number: int
    requested_precursor_number: int
    precursor_number_capped: bool
    focal_index: int
    sampled_indices: Tuple[int, ...]
    precursor_quotas: Tuple[int, ...]
    expanded_state: Tuple[CellState, ...]
    precursor_wraparound: bool
    precursor_min_index: int | None
    precursor_max_index: int | None



def sample_branch_niche(
    *,
    ring_state: Sequence[CellState],
    branch_id: str,
    child_id: str,
    event_time: Time2D,
    precursor_number: int,
    niche_size_m: int,
    rng: random.Random,
    focal_index: Optional[int] = None,
) -> BranchSelection:
    m = int(niche_size_m)
    if m < 1:
        raise ValueError('niche_size_m must be >= 1')
    recruitment = select_contiguous_precursors(
        ring_state,
        requested_cells=int(precursor_number),
        max_realized_cells=int(m),
        rng=rng,
        focal_index=focal_index,
    )
    quotas = allocate_equalized_quotas(m, recruitment.realized_cells)
    expanded: List[CellState] = []
    for cell, q in zip(recruitment.sampled_state, quotas):
        expanded.extend([cell] * int(q))
    return BranchSelection(
        branch_id=str(branch_id),
        child_id=str(child_id),
        event_time=(int(event_time[0]), int(event_time[1])),
        niche_size_m=m,
        precursor_number=int(recruitment.realized_cells),
        requested_precursor_number=int(recruitment.requested_cells),
        precursor_number_capped=bool(recruitment.capped),
        focal_index=int(recruitment.focal_index),
        sampled_indices=tuple(int(i) for i in recruitment.sampled_indices),
        precursor_quotas=tuple(int(q) for q in quotas),
        expanded_state=tuple(expanded),
        precursor_wraparound=bool(recruitment.wraparound),
        precursor_min_index=recruitment.min_index,
        precursor_max_index=recruitment.max_index,
    )
