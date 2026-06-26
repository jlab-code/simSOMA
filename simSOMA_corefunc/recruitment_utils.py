from __future__ import annotations

"""Shared helpers for fixed-size recruitment and exact founder expansion.

This module centralizes three operations now used across the pipeline:
- contiguous precursor selection from a ring
- exact integer quota allocation to a fixed target size
- exact balanced binary expansion from founders to a requested number of leaves
"""

from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple
import math
import random

import numpy as np

from self_renewal import CellState, SelfRenewalParams, SelfRenewalSimulator


@dataclass(frozen=True)
class ContiguousRecruitment:
    requested_cells: int
    realized_cells: int
    capped: bool
    focal_index: int
    sampled_indices: Tuple[int, ...]
    sampled_state: Tuple[CellState, ...]
    wraparound: bool
    min_index: int | None
    max_index: int | None


def _round_half_up(x: float) -> int:
    return int(math.floor(float(x) + 0.5))


def contiguous_block_indices(n: int, focal: int, k: int) -> List[int]:
    if n <= 0:
        raise ValueError('n must be positive')
    if k <= 0 or k > n:
        raise ValueError('k must be in 1..n')
    if focal < 0 or focal >= n:
        raise ValueError('focal out of bounds')
    if k == 1:
        return [focal]
    half = k // 2
    if k % 2 == 1:
        start = focal - half
    else:
        start = focal - (half - 1)
    return [((start + i) % n) for i in range(k)]


def _wraparound_from_indices(indices: Sequence[int], n: int) -> bool:
    if len(indices) <= 1:
        return False
    s = sorted(int(i) for i in indices)
    return (s[-1] - s[0] + 1) != len(s)


def select_contiguous_precursors(
    ring_state: Sequence[CellState],
    *,
    requested_cells: int,
    max_realized_cells: int,
    rng: random.Random,
    focal_index: Optional[int] = None,
) -> ContiguousRecruitment:
    N = len(ring_state)
    if N <= 0:
        raise ValueError('ring_state must be non-empty')
    req = int(requested_cells)
    if req < 1:
        raise ValueError('requested_cells must be >= 1')
    max_realized = int(max_realized_cells)
    if max_realized < 1:
        raise ValueError('max_realized_cells must be >= 1')
    realized = min(req, max_realized, N)
    focal = int(focal_index) if focal_index is not None else int(rng.randrange(N))
    idx = contiguous_block_indices(N, focal, realized)
    sampled = tuple(ring_state[i] for i in idx)
    wrap = _wraparound_from_indices(idx, N)
    sorted_idx = sorted(int(i) for i in idx)
    return ContiguousRecruitment(
        requested_cells=req,
        realized_cells=realized,
        capped=bool(realized < req),
        focal_index=focal,
        sampled_indices=tuple(int(i) for i in idx),
        sampled_state=sampled,
        wraparound=wrap,
        min_index=(sorted_idx[0] if sorted_idx else None),
        max_index=(sorted_idx[-1] if sorted_idx else None),
    )


def allocate_equalized_quotas(total_cells: int, n_founders: int) -> List[int]:
    total = int(total_cells)
    n = int(n_founders)
    if total < 1:
        raise ValueError('total_cells must be >= 1')
    if n < 1:
        raise ValueError('n_founders must be >= 1')
    if n > total:
        raise ValueError('n_founders cannot exceed total_cells when each founder must contribute at least one output cell')
    base = total // n
    rem = total % n
    return [base + (1 if i < rem else 0) for i in range(n)]


def expand_founders_to_exact_leaves(
    params: SelfRenewalParams,
    founders: Sequence[CellState],
    quotas: Sequence[int],
    *,
    rng: Optional[np.random.Generator] = None,
) -> tuple[list[CellState], list[int]]:
    if len(founders) != len(quotas):
        raise ValueError('founders and quotas must have the same length')
    sim = SelfRenewalSimulator(rng=rng if rng is not None else np.random.default_rng())

    def _grow(cell: CellState, n_leaves: int, depth: int) -> tuple[list[CellState], list[int]]:
        n = int(n_leaves)
        if n < 1:
            raise ValueError('quota leaves must be >= 1')
        if n == 1:
            return [cell], [depth]
        left_n = n // 2
        right_n = n - left_n
        left = sim.divide(params, cell)
        right = sim.divide(params, cell)
        leaves_l, depths_l = _grow(left, left_n, depth + 1)
        leaves_r, depths_r = _grow(right, right_n, depth + 1)
        return leaves_l + leaves_r, depths_l + depths_r

    leaves: list[CellState] = []
    depths: list[int] = []
    for founder, quota in zip(founders, quotas):
        sub_leaves, sub_depths = _grow(founder, int(quota), 0)
        leaves.extend(sub_leaves)
        depths.extend(sub_depths)
    return leaves, depths
