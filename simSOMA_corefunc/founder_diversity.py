"""Clonal composition of branch and organ founders (deterministic; no random calls).

Definition (local_sector_v1, simSOMA >= 0.2.0)
----------------------------------------------
At every developmental event the current ASC niche (m ordered positions) is amplified
into the SAM-boundary ring of C cells. Niche position i seeds one contiguous clonal
sector S_i of the ring (sizes from the equalized quota rule, differing by at most one
cell). A founding event samples a contiguous block of P_eff precursor cells
(P_eff = min(P_b, m, C) for a branch, min(P_o, C, O) for an organ).

For one event, with n_i = number of founders from sector S_i and p_i = n_i / P_eff:

    founder_sector_count  K       = #{i : n_i > 0}     (clonal sectors spanned)
    founder_polyclonal            = K > 1              (1 = polyclonal, 0 = monoclonal)
    founder_effective_sectors     N_eff = 1 / sum_i p_i^2
    founder_diversity     d       = (N_eff - 1) / (m - 1)    (0 if m = 1)
    founder_dominant_fraction     = max_i p_i

"Sector" is local: it refers to the niche positions at the time of the event, not to the
founding cells of the root niche. Run-level targets are event means, e.g.
pi_B = mean founder_polyclonal over branch events, pi_O over organ events.

Under uniform focal placement (no phyllotaxy) the expected values depend only on
(P_eff, C, m). For C divisible by m and P_eff <= C/m + 1:

    Pr(polyclonal) = min(1, (P_eff - 1) m / C),   E[K] = 1 + Pr(polyclonal),

i.e. on the compound quantity phi = (P_eff - 1) m / C. polyclonal_founding_expectation()
computes the exact values for any (P_eff, C, m) by enumerating focal positions.

The previous definition (simSOMA <= 0.1.x) counted distinct ROOT lineage labels. Those
labels are assigned once at the root and never reset, so the statistic measured whether
founders descend from different root-niche cells; it collapses with turnover (rho > 0)
and branch order. It is retained as root_lineage_* fields for comparison only and must
not be used as an inferential target.
"""
from __future__ import annotations

import bisect
from collections import Counter
from typing import Any, Iterable, Sequence

from recruitment_utils import allocate_equalized_quotas, contiguous_block_indices

FOUNDER_DEFINITION_VERSION = "local_sector_v1"


def sector_upper_bounds(boundary_cells: int, niche_size_m: int) -> list[int]:
    """Cumulative upper bounds of the m clonal sectors on a ring of C cells."""
    quotas = allocate_equalized_quotas(int(boundary_cells), int(niche_size_m))
    out, acc = [], 0
    for q in quotas:
        acc += int(q)
        out.append(acc)
    return out


def sector_of_indices(indices: Sequence[int], boundary_cells: int, niche_size_m: int) -> list[int]:
    bounds = sector_upper_bounds(boundary_cells, niche_size_m)
    return [int(bisect.bisect_right(bounds, int(i))) for i in indices]


def _composition_stats(counts: Counter, niche_size_m: int) -> dict:
    total = int(sum(counts.values()))
    if total <= 0:
        return {"sector_count": 0, "effective": 0.0, "diversity": 0.0, "dominant": 0.0}
    props = [c / total for c in counts.values()]
    effective = 1.0 / sum(p * p for p in props)
    m = int(niche_size_m)
    diversity = 0.0 if m <= 1 else min(1.0, max(0.0, (effective - 1.0) / (m - 1.0)))
    return {"sector_count": len(counts), "effective": float(effective), "diversity": float(diversity),
            "dominant": float(max(props))}


def summarize_founder_sectors(sampled_indices: Sequence[int], boundary_cells: int, niche_size_m: int) -> dict:
    """Local clonal-sector composition of the realized founder block (definition above)."""
    sectors = sector_of_indices(sampled_indices, boundary_cells, niche_size_m)
    counts = Counter(sectors)
    s = _composition_stats(counts, niche_size_m)
    return {
        "founder_definition": FOUNDER_DEFINITION_VERSION,
        "founder_sector_counts": {int(k): int(v) for k, v in sorted(counts.items())},
        # kept for interface compatibility (fitSOMA <= 0.3.19 reads this key); same content
        "founder_lineage_counts": {int(k): int(v) for k, v in sorted(counts.items())},
        "founder_sector_count": int(s["sector_count"]),
        "founder_polyclonal": int(s["sector_count"] > 1),
        "founder_effective_sectors": s["effective"],
        "founder_diversity": s["diversity"],
        "founder_dominant_fraction": s["dominant"],
    }


def polyclonal_founding_expectation(precursors_eff: int, boundary_cells: int, niche_size_m: int) -> dict:
    """Exact expectation under a uniformly placed contiguous block (no phyllotaxy).

    Returns probability of polyclonal founding, expected sector count, and the compound
    quantity phi = (P_eff - 1) m / C.
    """
    P, C, m = int(precursors_eff), int(boundary_cells), int(niche_size_m)
    if P < 1 or C < 1 or m < 1 or P > C:
        raise ValueError("need 1 <= P_eff <= C and m >= 1")
    bounds = sector_upper_bounds(C, m)
    sector = [bisect.bisect_right(bounds, i) for i in range(C)]
    n_poly, k_sum = 0, 0
    for f in range(C):
        k = len({sector[i] for i in contiguous_block_indices(C, f, P)})
        n_poly += int(k > 1)
        k_sum += k
    return {"probability": n_poly / C, "expected_sector_count": k_sum / C, "phi": (P - 1) * m / C}


def summarize_founder_lineages(cells: Iterable[Any], niche_size_m: int) -> dict:
    """DEPRECATED root-lineage statistic (old definition); returned with root_lineage_* keys."""
    counts = Counter(int(getattr(cell, "lineage_id")) for cell in cells)
    s = _composition_stats(counts, niche_size_m)
    return {
        "root_lineage_counts": {int(k): int(v) for k, v in sorted(counts.items())},
        "root_lineage_sector_count": int(s["sector_count"]),
        "root_lineage_effective_sectors": s["effective"],
        "root_lineage_diversity": s["diversity"],
        "root_lineage_dominant_fraction": s["dominant"],
    }


FOUNDER_EVENT_FIELDS = (
    "founder_definition", "founder_sector_counts", "founder_lineage_counts", "founder_sector_count",
    "founder_polyclonal", "founder_effective_sectors", "founder_diversity", "founder_dominant_fraction",
    "root_lineage_counts", "root_lineage_sector_count", "root_lineage_effective_sectors",
    "root_lineage_diversity", "root_lineage_dominant_fraction",
)
