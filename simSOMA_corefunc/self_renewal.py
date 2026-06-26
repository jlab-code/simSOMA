"""self_renewal.py

Self-renewal (ASC niche) simulator.

Design goal
-----------
State-based niche simulator with no global lineage store. The niche state is a
list of `CellState` objects of length `m`.

Time semantics
--------------
We keep a 2D time stamp for bookkeeping: `Time2D = (t_age, t_dev)`.

- Self-renewal steps along a branch advance `t_age` and always use `t_dev = 0`.
- Fast developmental bursts (pre-branching, organ expansion) keep `t_age` fixed
  and advance `t_dev`.

Mutation model
--------------
- Infinite-sites-by-ID approximation: each division edge receives Poisson(mu_div)
  new mutation IDs.
- Mutation IDs are drawn as random 63-bit integers from the provided RNG.
- Genotypes are represented as a frozenset[int] of mutation IDs.

Competition model
-----------------
Each branch has one branch-local favored clone, seeded from a single niche
position at branch start. `comp_label=1` marks membership in that favored clone
within the current branch; `comp_label=0` otherwise. Descendants inherit the
label within the branch, but the label is reset at each new branch founding
by the wrapper.

With probability `rho` per self-renewal round, one displacement occurs:
- the displacer is chosen non-uniformly using the branch-local bias `w` in [0,1]
- the victim is chosen using a one-parameter nearest-neighbor mixture on the
  ring, controlled by `victim_locality` in [0,1].
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

Time2D = Tuple[int, int]  # (t_age, t_dev)


@dataclass(frozen=True)
class CellState:
    genotype: frozenset[int]
    lineage_id: int
    comp_label: int = 0  # 1 if in the favored branch-local clone, else 0


@dataclass(frozen=True)
class SelfRenewalParams:
    m: int
    rho: float

    # Mutation and time-mapping parameters (user-facing)
    mu_year: float        # mutations per topology unit (year or meter)
    kappa_sr: float       # self-renewal divisions per topology unit

    # Optional override: if provided, used directly instead of mu_year/kappa_sr
    mu_div: Optional[float] = None

    # Victim selection locality in [0,1]
    # 0.0 -> uniform among all non-displacer sites
    # 1.0 -> uniform among nearest neighbors of the displacer
    victim_locality: float = 0.0

    # Branch-local displacement advantage in [0,1]
    # 0.0 -> neutral displacer choice
    # 1.0 -> only favored-clone cells can displace (if any are present)
    branch_comp_bias: float = 0.0


@dataclass(frozen=True)
class SelfRenewalInputs:
    branch_id: str
    T: int
    age_offset: int = 0
    init_state: Optional[Sequence[CellState]] = None
    snapshot_times: Optional[Sequence[int]] = None  # local (0..T)


@dataclass(frozen=True)
class DisplacementEvent:
    branch_id: str
    time: Time2D
    displacer_site: int
    victim_site: int


@dataclass
class SelfRenewalOutputs:
    branch_id: str
    m: int
    age_offset: int
    final_state: List[CellState]
    snapshots: Dict[int, List[CellState]]  # absolute t_age -> state
    displacement_events: List[DisplacementEvent]
    num_displacements: int


def _ring_distance(i: int, j: int, m: int) -> int:
    d = abs(i - j)
    return min(d, m - d)


def _nearest_neighbors(i: int, m: int) -> List[int]:
    if m <= 1:
        raise ValueError("No victim exists when m <= 1")
    if m == 2:
        return [1 - i]
    candidates = [j for j in range(m) if j != i]
    dists = [_ring_distance(i, j, m) for j in candidates]
    dmin = min(dists)
    return [j for j, d in zip(candidates, dists) if d == dmin]


def _choose_victim_site(i: int, m: int, victim_locality: float, rng: np.random.Generator) -> int:
    if m <= 1:
        raise ValueError("No victim exists when m <= 1")
    ell = float(victim_locality)
    if (ell < 0.0) or (ell > 1.0):
        raise ValueError("victim_locality must be in [0,1]")

    candidates = [j for j in range(m) if j != i]
    if m == 2:
        return candidates[0]

    if ell <= 0.0:
        return int(rng.choice(candidates))
    if ell >= 1.0:
        return int(rng.choice(_nearest_neighbors(i, m)))

    if rng.random() < ell:
        return int(rng.choice(_nearest_neighbors(i, m)))
    return int(rng.choice(candidates))


def assign_branch_comp_seed(state: Sequence[CellState], favored_site: int) -> List[CellState]:
    if len(state) == 0:
        raise ValueError("state must be non-empty")
    if favored_site < 0 or favored_site >= len(state):
        raise ValueError("favored_site out of bounds")
    out: List[CellState] = []
    for i, cell in enumerate(state):
        out.append(
            CellState(
                genotype=cell.genotype,
                lineage_id=int(cell.lineage_id),
                comp_label=1 if i == favored_site else 0,
            )
        )
    return out


class SelfRenewalSimulator:
    """Stateless (history-free) simulator using a provided RNG."""

    def __init__(self, rng: Optional[np.random.Generator] = None):
        self.rng = rng if rng is not None else np.random.default_rng()

    def mu_div(self, params: SelfRenewalParams) -> float:
        if params.mu_div is not None:
            return float(params.mu_div)
        if params.kappa_sr <= 0:
            raise ValueError("kappa_sr must be > 0 to derive mu_div")
        return float(params.mu_year) / float(params.kappa_sr)

    def _draw_mutation_ids(self, k: int) -> List[int]:
        if k <= 0:
            return []
        arr = self.rng.integers(1, 2**63 - 1, size=int(k), dtype=np.int64)
        return [int(x) for x in arr.tolist()]

    def divide(self, params: SelfRenewalParams, parent: CellState) -> CellState:
        lam = self.mu_div(params)
        k = int(self.rng.poisson(lam))
        if k <= 0:
            return CellState(genotype=parent.genotype, lineage_id=parent.lineage_id, comp_label=parent.comp_label)
        mids = self._draw_mutation_ids(k)
        g = set(parent.genotype)
        g.update(mids)
        return CellState(genotype=frozenset(g), lineage_id=parent.lineage_id, comp_label=parent.comp_label)

    def initialize_founders(self, params: SelfRenewalParams) -> List[CellState]:
        m = int(params.m)
        if m <= 0:
            raise ValueError("m must be positive")
        return [CellState(genotype=frozenset(), lineage_id=int(i), comp_label=0) for i in range(m)]

    def _choose_displacer_site(self, children: Sequence[CellState], w: float) -> int:
        if (w < 0.0) or (w > 1.0):
            raise ValueError("branch_comp_bias must be in [0,1]")
        weights = np.array([1.0 if int(c.comp_label) == 1 else 1.0 - w for c in children], dtype=float)
        sw = float(weights.sum())
        if sw <= 0.0 or (not np.isfinite(sw)):
            return int(self.rng.integers(0, len(children)))
        return int(self.rng.choice(np.arange(len(children)), p=(weights / sw)))

    def simulate_segment(self, params: SelfRenewalParams, inputs: SelfRenewalInputs) -> SelfRenewalOutputs:
        m = int(params.m)
        T = int(inputs.T)
        if T < 0:
            raise ValueError("T must be >= 0")

        if inputs.init_state is None:
            state = self.initialize_founders(params)
        else:
            state = list(inputs.init_state)
            if len(state) != m:
                raise ValueError(f"init_state length must be m={m}")

        snapshot_set = set(int(t) for t in inputs.snapshot_times) if inputs.snapshot_times else set()
        snapshots: Dict[int, List[CellState]] = {}
        events: List[DisplacementEvent] = []
        num_disp = 0

        if 0 in snapshot_set:
            snapshots[int(inputs.age_offset)] = list(state)

        for t_local in range(1, T + 1):
            t_abs = int(inputs.age_offset) + t_local
            time: Time2D = (t_abs, 0)

            # every site divides (synchronous update)
            children = [self.divide(params, parent=state[i]) for i in range(m)]

            if float(params.rho) > 0 and self.rng.random() < float(params.rho):
                displacer = self._choose_displacer_site(children, float(params.branch_comp_bias))

                victim = _choose_victim_site(displacer, m, float(params.victim_locality), self.rng)

                children[victim] = children[displacer]
                events.append(
                    DisplacementEvent(
                        branch_id=str(inputs.branch_id),
                        time=time,
                        displacer_site=int(displacer),
                        victim_site=int(victim),
                    )
                )
                num_disp += 1

            state = children

            if t_local in snapshot_set:
                snapshots[t_abs] = list(state)

        return SelfRenewalOutputs(
            branch_id=str(inputs.branch_id),
            m=m,
            age_offset=int(inputs.age_offset),
            final_state=list(state),
            snapshots=snapshots,
            displacement_events=events,
            num_displacements=int(num_disp),
        )
