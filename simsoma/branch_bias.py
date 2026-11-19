"""branch_bias.py

Branch-local displacement-bias policy.

This module cleanly separates:
- how a realized branch-local displacement bias value w_b is assigned, and
- how the favored seed cell for a branch is chosen.

The biology is implemented elsewhere; this module only samples branch-specific
metadata used by the wrapper.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Literal, Optional
import json

import numpy as np

BiasMode = Literal["fixed", "draw", "map"]


def _clamp01(x: float) -> float:
    y = float(x)
    if y < 0.0 or y > 1.0:
        raise ValueError("bias values must lie in [0,1]")
    return y


def _beta_ab_from_mean_kappa(mean: float, kappa: float) -> tuple[float, float]:
    m = _clamp01(mean)
    if m <= 0.0:
        return 0.0, 0.0
    if m >= 1.0:
        return np.inf, np.inf
    if kappa <= 0:
        raise ValueError("draw_kappa must be > 0")
    return m * float(kappa), (1.0 - m) * float(kappa)


@dataclass(frozen=True)
class BranchBiasParams:
    seed: int
    bias_mode: BiasMode = "fixed"
    fixed_value: float = 0.0
    draw_mean: float = 0.0
    draw_kappa: float = 50.0
    map_values: Optional[Dict[str, float]] = None


class BranchBiasState:
    def __init__(self, params: BranchBiasParams):
        self.params = params
        self.rng = np.random.default_rng(int(params.seed))
        self._w_by_branch: Dict[str, float] = {}
        self._seed_by_branch: Dict[str, int] = {}

    @staticmethod
    def load_map_json(path: str | Path) -> Dict[str, float]:
        raw = json.loads(Path(path).read_text(encoding="utf-8"))
        if not isinstance(raw, dict):
            raise ValueError("branch bias map JSON must be a dict of branch_id -> value")
        out: Dict[str, float] = {}
        for k, v in raw.items():
            out[str(k)] = _clamp01(float(v))
        return out

    def branch_bias(self, branch_id: str) -> float:
        bid = str(branch_id)
        if bid in self._w_by_branch:
            return self._w_by_branch[bid]

        mode = str(self.params.bias_mode).lower()
        if mode == "fixed":
            w = _clamp01(self.params.fixed_value)
        elif mode == "draw":
            a, b = _beta_ab_from_mean_kappa(float(self.params.draw_mean), float(self.params.draw_kappa))
            if a == 0.0 and b == 0.0:
                w = 0.0
            elif np.isinf(a) and np.isinf(b):
                w = 1.0
            else:
                w = float(self.rng.beta(a, b))
        elif mode == "map":
            mv = self.params.map_values or {}
            if bid in mv:
                w = _clamp01(float(mv[bid]))
            elif "default" in mv:
                w = _clamp01(float(mv["default"]))
            else:
                raise KeyError(f"No branch bias value found for branch '{bid}' in map mode")
        else:
            raise ValueError(f"Unknown bias_mode: {self.params.bias_mode}")

        self._w_by_branch[bid] = float(w)
        return float(w)

    def branch_seed(self, branch_id: str, niche_size: int) -> int:
        bid = str(branch_id)
        n = int(niche_size)
        if n <= 0:
            raise ValueError("niche_size must be positive")
        if bid not in self._seed_by_branch:
            self._seed_by_branch[bid] = int(self.rng.integers(0, n))
        seed_idx = int(self._seed_by_branch[bid])
        if seed_idx >= n:
            raise ValueError(f"Cached seed index {seed_idx} out of range for branch {bid} with niche size {n}")
        return seed_idx

    def realize_all(self, *, branch_ids: list[str], niche_size_m: int) -> None:
        for bid in branch_ids:
            _ = self.branch_bias(bid)
            _ = self.branch_seed(bid, niche_size=niche_size_m)
