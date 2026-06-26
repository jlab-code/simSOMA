"""phyllotaxy.py

Topology annotation utilities for site-based phyllotaxy.

Design
------
- Phyllotaxy is assigned to *sites* along each axis/branch, not directly to
  individual events.
- Events tied at the same branch-local coordinate belong to the same site and
  inherit the same angle.
- Supported modes intentionally stay lightweight and topology-driven:
    * off        : do not annotate; downstream modules choose focal index freely
    * random     : one iid Uniform(0, 2pi) angle per site
    * spiral     : constant divergence angle between successive sites
    * distichous : spiral with 180 degree divergence
    * tristichous: spiral with 120 degree divergence
    * decussate  : site-based approximation using 90 degree divergence

Notes
-----
The site abstraction is intentionally generic. If a leaf and its axillary
branch are encoded as tied events on the same axis, they inherit one shared
site angle and therefore the same focal index on the parent ring.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple
import copy
import math
import random


_TWO_PI = 2.0 * math.pi
_DEFAULT_TIE_TOL = 1e-9


@dataclass(frozen=True)
class PhyllotaxyConfig:
    mode: str = "off"
    divergence_angle_deg: Optional[float] = None
    tie_tol: float = _DEFAULT_TIE_TOL


def _normalize_mode(mode: Optional[str]) -> str:
    s = str(mode or "off").strip().lower()
    aliases = {
        "none": "off",
        "disabled": "off",
        "uniform": "random",
        "iid": "random",
        "custom": "spiral",
        "custom_angle": "spiral",
        "alternate": "distichous",
        "tri": "tristichous",
    }
    return aliases.get(s, s)


def normalize_config(cfg: Optional[Dict[str, Any]]) -> PhyllotaxyConfig:
    if not cfg:
        return PhyllotaxyConfig(mode="off")

    mode = _normalize_mode(cfg.get("mode", "off"))
    tie_tol = float(cfg.get("tie_tol", _DEFAULT_TIE_TOL))
    if tie_tol < 0:
        raise ValueError("phyllotaxy.tie_tol must be >= 0")

    divergence_angle_deg = cfg.get("divergence_angle_deg")
    if mode == "spiral":
        if divergence_angle_deg is None:
            raise ValueError("phyllotaxy mode='spiral' requires divergence_angle_deg")
        divergence_angle_deg = float(divergence_angle_deg)
    elif mode == "distichous":
        divergence_angle_deg = 180.0
    elif mode == "tristichous":
        divergence_angle_deg = 120.0
    elif mode == "decussate":
        # Site-based approximation. Opposite members of a decussate pair would
        # typically be represented as separate sites in richer topologies.
        divergence_angle_deg = 90.0
    elif mode in {"off", "random"}:
        divergence_angle_deg = None
    else:
        raise ValueError(
            "Unknown phyllotaxy mode. Allowed: off, random, spiral, distichous, tristichous, decussate"
        )

    return PhyllotaxyConfig(
        mode=mode,
        divergence_angle_deg=divergence_angle_deg,
        tie_tol=tie_tol,
    )


def _site_coordinate(event: Dict[str, Any]) -> float:
    if "time_obs" in event:
        return float(event["time_obs"])
    return float(event["time"])


def _site_sort_key(event: Dict[str, Any]) -> Tuple[float, str, str]:
    return (
        float(_site_coordinate(event)),
        str(event.get("type", "")).upper(),
        str(event.get("target_id", "")),
    )


def _same_site(a: float, b: float, tol: float) -> bool:
    return math.isclose(a, b, rel_tol=0.0, abs_tol=float(tol))


def _normalize_angle_rad(x: float) -> float:
    y = float(x) % _TWO_PI
    if y < 0.0:
        y += _TWO_PI
    return y


def angle_rad_to_deg(angle_rad: float) -> float:
    return float(_normalize_angle_rad(angle_rad) * 180.0 / math.pi)


def angle_to_focal_index(angle_rad: float, ring_size: int) -> int:
    n = int(ring_size)
    if n <= 0:
        raise ValueError("ring_size must be positive")
    frac = _normalize_angle_rad(angle_rad) / _TWO_PI
    return int(math.floor(n * frac + 0.5)) % n


def _site_angle_rad(cfg: PhyllotaxyConfig, *, phase0: float, site_rank: int, rng: random.Random) -> float:
    if cfg.mode == "random":
        return _normalize_angle_rad(rng.random() * _TWO_PI)
    if cfg.mode == "off":
        raise ValueError("site angle requested in phyllotaxy mode='off'")
    assert cfg.divergence_angle_deg is not None
    alpha = float(cfg.divergence_angle_deg) * math.pi / 180.0
    return _normalize_angle_rad(float(phase0) + float(site_rank - 1) * alpha)


def annotate_topology_sites(
    topo: Dict[str, Any],
    *,
    config: Optional[Dict[str, Any]] = None,
    seed: Optional[int] = None,
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """Return a deep-copied topology with per-event phyllotaxy annotations.

    Per-event fields added when mode != 'off':
      - phyllo_site_id
      - phyllo_site_rank
      - phyllo_angle_rad
      - phyllo_angle_deg

    Branch-level metadata added when mode != 'off':
      - phyllo_phase_rad
      - phyllo_phase_deg
    """
    cfg = normalize_config(config)
    out = copy.deepcopy(topo)
    summary: Dict[str, Any] = {
        "mode": cfg.mode,
        "divergence_angle_deg": (None if cfg.divergence_angle_deg is None else float(cfg.divergence_angle_deg)),
        "tie_tol": float(cfg.tie_tol),
        "annotated": bool(cfg.mode != "off"),
        "site_counts_by_branch": {},
        "phase_deg_by_branch": {},
    }
    if cfg.mode == "off":
        return out, summary

    rng = random.Random(seed)

    for branch_id, spec in out.get("branches", {}).items():
        events = list(spec.get("events", []))
        events_sorted = sorted(events, key=_site_sort_key)
        phase0 = _normalize_angle_rad(rng.random() * _TWO_PI)
        spec["phyllo_phase_rad"] = float(phase0)
        spec["phyllo_phase_deg"] = float(angle_rad_to_deg(phase0))
        summary["phase_deg_by_branch"][str(branch_id)] = float(spec["phyllo_phase_deg"])

        site_rank = 0
        prev_coord: Optional[float] = None
        current_angle_rad: Optional[float] = None
        current_site_id: Optional[str] = None
        for e in events_sorted:
            coord = float(_site_coordinate(e))
            if prev_coord is None or not _same_site(coord, prev_coord, cfg.tie_tol):
                site_rank += 1
                current_site_id = f"{branch_id}:S{site_rank}"
                current_angle_rad = _site_angle_rad(cfg, phase0=phase0, site_rank=site_rank, rng=rng)
                prev_coord = coord
            assert current_angle_rad is not None and current_site_id is not None
            e["phyllo_site_id"] = str(current_site_id)
            e["phyllo_site_rank"] = int(site_rank)
            e["phyllo_angle_rad"] = float(current_angle_rad)
            e["phyllo_angle_deg"] = float(angle_rad_to_deg(current_angle_rad))

        summary["site_counts_by_branch"][str(branch_id)] = int(site_rank)

    return out, summary
