"""Deterministic summaries of realized founder-lineage diversity.

This module consumes already generated CellState objects and makes no random
calls.  The normalized diversity is

    d = (N_eff - 1) / (m - 1),

where N_eff is the inverse-Simpson effective number of represented ASC
lineages.  For m <= 1, d is zero.
"""
from __future__ import annotations

from collections import Counter
from typing import Any, Iterable


def summarize_founder_lineages(cells: Iterable[Any], niche_size_m: int) -> dict:
    counts = Counter(int(getattr(cell, "lineage_id")) for cell in cells)
    total = int(sum(counts.values()))
    if total <= 0:
        effective = 0.0
        diversity = 0.0
        dominant_fraction = 0.0
    else:
        proportions = [float(count) / float(total) for count in counts.values()]
        denom = sum(p * p for p in proportions)
        effective = 1.0 / denom if denom > 0.0 else 0.0
        m = int(niche_size_m)
        diversity = 0.0 if m <= 1 else (effective - 1.0) / (m - 1.0)
        diversity = min(1.0, max(0.0, float(diversity)))
        dominant_fraction = max(proportions)
    return {
        "founder_lineage_counts": {int(k): int(v) for k, v in sorted(counts.items())},
        "founder_sector_count": int(len(counts)),
        "founder_effective_sectors": float(effective),
        "founder_diversity": float(diversity),
        "founder_dominant_fraction": float(dominant_fraction),
    }
