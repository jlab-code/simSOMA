"""Read-level observation model (see package docstring for provenance and versioning).

Data flow for one site i and one sample s
-----------------------------------------
    assay VAF   v_is   = phase_s * c_ks * f_is         (computed by the caller of observe())
    read prob.  p_is   = v_is (1 - e) + (1 - v_is) e   (symmetric sequencing error e)
    depth       D_is   ~ depth model
    alt reads   A_is   ~ Binomial(D_is, p_is)  or  Binomial(D_is, Beta(p k, (1-p) k))
    call        D_is >= min_depth  and  A_is >= min_alt_reads  and  A_is/D_is >= min_observed_vaf
    ascertain   keep site if called in >= 1 sample (optional)

Background artefact sites (optional) have v_is = q_i, the same in every sample, and pass
through the same depth / read / caller steps.

Depth models
------------
    fixed                 D_is = mean_s
    normal                D_is = round(Normal(mean_s, sd_s)), floored at `minimum`
                          (alias: sample_design; the historical simSOMA default)
    poisson               D_is ~ Poisson(mean_s)
    lognormal_site_sample D_is ~ Poisson(mean_s g_i e_is),
                          log g_i ~ N(-sg^2/2, sg) truncated at g_i <= max_site_factor,
                          log e_is ~ N(-se^2/2, se)
"""
from __future__ import annotations

import copy
import hashlib
import json
from typing import Any, Mapping, Optional, Sequence

import numpy as np

OBSERVATION_MODEL_VERSION = "1.0.0"

DEFAULTS: dict[str, Any] = {
    "depth": {
        "mode": "normal",
        "mean": 100.0,
        "sd": 15.0,
        "minimum": 1,
        "site_sdlog": 0.58,
        "sample_sdlog": 0.235,
        "max_site_factor": 3.0,
    },
    "reads": {
        "type": "beta_binomial",
        "concentration": 200.0,
        "sequencing_error": 0.001,
    },
    "background": {
        "n_sites": 0,
        "distribution": "exp",
        "mean": 0.015,
        "clip": 0.45,
    },
    "caller": {
        "min_depth": 20,
        "min_alt_reads": 3,
        "min_observed_vaf": 0.0,
        "retain_called_any": True,
    },
}

DEPTH_MODES = ("fixed", "normal", "poisson", "lognormal_site_sample")
_DEPTH_ALIASES = {"sample_design": "normal", "lognormal": "lognormal_site_sample"}
READ_TYPES = ("binomial", "beta_binomial")


def _merge(base: dict, upd: Optional[Mapping]) -> dict:
    out = copy.deepcopy(base)
    for k, v in (upd or {}).items():
        if isinstance(v, Mapping) and isinstance(out.get(k), dict):
            out[k] = _merge(out[k], v)
        else:
            out[k] = v
    return out


def normalize_config(config: Optional[Mapping[str, Any]] = None) -> dict[str, Any]:
    """Fill defaults and validate. Raises ValueError with the offending key."""
    cfg = _merge(DEFAULTS, config)
    unknown = set(cfg) - set(DEFAULTS)
    if unknown:
        raise ValueError(f"unknown observation-model sections: {sorted(unknown)}")
    for sec in DEFAULTS:
        extra = set(cfg[sec]) - set(DEFAULTS[sec])
        if extra:
            raise ValueError(f"unknown keys in observation-model section '{sec}': {sorted(extra)}")
    d = cfg["depth"]
    d["mode"] = _DEPTH_ALIASES.get(str(d["mode"]), str(d["mode"]))
    if d["mode"] not in DEPTH_MODES:
        raise ValueError(f"depth.mode must be one of {list(DEPTH_MODES)}, got {d['mode']!r}")
    if float(d["mean"]) <= 0 or float(d["sd"]) < 0 or int(d["minimum"]) < 0:
        raise ValueError("depth.mean must be > 0, depth.sd >= 0, depth.minimum >= 0")
    if float(d["max_site_factor"]) <= 0 or float(d["site_sdlog"]) < 0 or float(d["sample_sdlog"]) < 0:
        raise ValueError("depth.max_site_factor must be > 0 and sdlog values >= 0")
    r = cfg["reads"]
    if r["type"] not in READ_TYPES:
        raise ValueError(f"reads.type must be one of {list(READ_TYPES)}, got {r['type']!r}")
    if not (0.0 <= float(r["sequencing_error"]) < 0.5):
        raise ValueError("reads.sequencing_error must be in [0, 0.5)")
    if float(r["concentration"]) <= 0:
        raise ValueError("reads.concentration must be > 0 (beta-binomial: rho = 1 / (concentration + 1))")
    b = cfg["background"]
    if int(b["n_sites"]) < 0 or float(b["mean"]) <= 0 or not (0 < float(b["clip"]) <= 1):
        raise ValueError("background.n_sites >= 0, background.mean > 0, 0 < background.clip <= 1")
    _parse_background_distribution(str(b["distribution"]), float(b["mean"]))
    c = cfg["caller"]
    _canonical_types(cfg)
    if int(c["min_depth"]) < 0 or int(c["min_alt_reads"]) < 0 or not (0 <= float(c["min_observed_vaf"]) <= 1):
        raise ValueError("caller thresholds: min_depth >= 0, min_alt_reads >= 0, 0 <= min_observed_vaf <= 1")
    return cfg


_INT_KEYS = {("depth", "minimum"), ("background", "n_sites"), ("caller", "min_depth"), ("caller", "min_alt_reads")}
_BOOL_KEYS = {("caller", "retain_called_any")}
_STR_KEYS = {("depth", "mode"), ("reads", "type"), ("background", "distribution")}


def _canonical_types(cfg: dict) -> None:
    """Cast values to canonical types so equal settings hash equally (60 == 60.0)."""
    for sec, vals in cfg.items():
        for k, v in vals.items():
            if (sec, k) in _INT_KEYS:
                vals[k] = int(v)
            elif (sec, k) in _BOOL_KEYS:
                vals[k] = bool(v)
            elif (sec, k) in _STR_KEYS:
                vals[k] = str(v)
            else:
                vals[k] = float(v)


def config_sha256(config: Mapping[str, Any]) -> str:
    """Digest of the normalized configuration together with the model version."""
    payload = {"version": OBSERVATION_MODEL_VERSION, "config": normalize_config(config)}
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


# ---------------------------------------------------------------------------------------- depth
def draw_depths(n_sites: int, n_samples: int, depth_cfg: Mapping[str, Any], rng: np.random.Generator,
                sample_means: Optional[Sequence[float]] = None,
                sample_sds: Optional[Sequence[float]] = None) -> np.ndarray:
    """Depth matrix (n_sites x n_samples). sample_means / sample_sds override depth.mean / depth.sd
    per sample (NaN entries fall back to the defaults)."""
    d = depth_cfg
    mode = _DEPTH_ALIASES.get(str(d["mode"]), str(d["mode"]))
    means = _per_sample(sample_means, float(d["mean"]), n_samples)
    sds = _per_sample(sample_sds, float(d["sd"]), n_samples)
    minimum = int(d["minimum"])
    if mode == "lognormal_site_sample":
        sg, se, cap = float(d["site_sdlog"]), float(d["sample_sdlog"]), float(d["max_site_factor"])
        g = np.empty(n_sites)
        filled = 0
        while filled < n_sites:                       # truncation by resampling (max-depth filter)
            x = np.exp(rng.normal(-sg ** 2 / 2, sg, max(4 * (n_sites - filled), 16)))
            x = x[x <= cap][: n_sites - filled]
            g[filled:filled + len(x)] = x
            filled += len(x)
        e = np.exp(rng.normal(-se ** 2 / 2, se, (n_sites, n_samples)))
        out = rng.poisson(means[None, :] * g[:, None] * e)
        return np.maximum(out, minimum).astype(int)
    out = np.empty((n_sites, n_samples), dtype=int)
    for s in range(n_samples):                         # column order kept for backward compatibility
        if mode == "fixed":
            col = np.repeat(int(round(means[s])), n_sites)
        elif mode == "normal":
            col = (np.repeat(int(round(means[s])), n_sites) if sds[s] <= 0
                   else np.rint(rng.normal(means[s], sds[s], size=n_sites)).astype(int))
        else:  # poisson
            col = rng.poisson(means[s], size=n_sites)
        out[:, s] = np.maximum(col, minimum)
    return out


def _per_sample(values, default: float, n: int) -> np.ndarray:
    if values is None:
        return np.full(n, default, dtype=float)
    v = np.asarray(values, dtype=float)
    if v.shape != (n,):
        raise ValueError(f"per-sample depth parameters must have length {n}")
    return np.where(np.isnan(v), default, v)


# ----------------------------------------------------------------------------------- background
def _parse_background_distribution(spec: str, mean: float):
    kind, *args = spec.split(":")
    if kind == "exp" and not args:
        return lambda rng, n: rng.exponential(mean, n)
    if kind == "gamma" and len(args) == 1:
        k = float(args[0]); return lambda rng, n: rng.gamma(k, mean / k, n)
    if kind == "lnorm" and len(args) == 1:
        sd = float(args[0]); return lambda rng, n: rng.lognormal(np.log(mean) - sd ** 2 / 2, sd, n)
    if kind == "mix" and len(args) == 2:
        w2, m2 = float(args[0]), float(args[1]); m1 = (mean - w2 * m2) / (1 - w2)
        if not (0 < w2 < 1) or m1 <= 0:
            raise ValueError("background mix:<w2>:<mean2> needs 0 < w2 < 1 and a positive first-component mean")
        return lambda rng, n: np.where(rng.random(n) < w2, rng.exponential(m2, n), rng.exponential(m1, n))
    raise ValueError(f"background.distribution must be exp | gamma:<shape> | lnorm:<sdlog> | mix:<w2>:<mean2>, got {spec!r}")


def draw_background_vafs(background_cfg: Mapping[str, Any], rng: np.random.Generator) -> np.ndarray:
    b = background_cfg
    n = int(b["n_sites"])
    if n == 0:
        return np.zeros(0)
    q = _parse_background_distribution(str(b["distribution"]), float(b["mean"]))(rng, n)
    return np.clip(q, 0.0, float(b["clip"]))


# --------------------------------------------------------------------------------------- reads
def draw_alt_reads(assay_vaf: np.ndarray, depth: np.ndarray, reads_cfg: Mapping[str, Any],
                   rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    """Returns (alt counts, read probability incl. sequencing error)."""
    eps = float(reads_cfg["sequencing_error"])
    p = np.clip(assay_vaf * (1.0 - eps) + (1.0 - assay_vaf) * eps, 0.0, 1.0)
    if reads_cfg["type"] == "binomial":
        return rng.binomial(depth, p), p
    k = float(reads_cfg["concentration"])
    safe = np.clip(p, 1e-10, 1.0 - 1e-10)
    return rng.binomial(depth, rng.beta(safe * k, (1.0 - safe) * k)), p


def call_variants(alt: np.ndarray, depth: np.ndarray, caller_cfg: Mapping[str, Any]):
    """Returns (callable, called, status) arrays."""
    vaf = np.divide(alt, depth, out=np.zeros_like(alt, dtype=float), where=depth > 0)
    callable_ = depth >= int(caller_cfg["min_depth"])
    called = callable_ & (alt >= int(caller_cfg["min_alt_reads"])) & (vaf >= float(caller_cfg["min_observed_vaf"]))
    status = np.full(called.shape, "BELOW_CALL_THRESHOLD", dtype=object)
    status[~callable_] = "NO_COVERAGE"
    status[called] = "CALLED"
    return callable_, called, status


# ------------------------------------------------------------------------------------- observe
def observe(assay_vaf: np.ndarray, config: Optional[Mapping[str, Any]], rng: np.random.Generator,
            sample_means: Optional[Sequence[float]] = None,
            sample_sds: Optional[Sequence[float]] = None) -> dict[str, Any]:
    """Apply the observation model to an assay-VAF matrix (somatic sites x samples).

    Background sites, if configured, are appended after the somatic rows. Returns a dict with
    numpy arrays (all n_total x n_samples unless noted): assay_vaf, depth, alt, read_probability,
    observed_vaf, callable, called, status; plus is_background (n_total,), keep (n_total,, the
    ascertainment mask), version and config_sha256.

    Random-number order: depths of somatic sites, reads of somatic sites, then background VAFs,
    depths and reads (so configurations without background reproduce the historical simSOMA
    transform exactly).
    """
    cfg = normalize_config(config)
    v = np.clip(np.asarray(assay_vaf, dtype=float), 0.0, 1.0)
    n_sites, n_samples = v.shape
    depth = draw_depths(n_sites, n_samples, cfg["depth"], rng, sample_means, sample_sds)
    alt, p = draw_alt_reads(v, depth, cfg["reads"], rng)
    is_bg = np.zeros(n_sites, dtype=bool)
    q = draw_background_vafs(cfg["background"], rng)
    if len(q):
        vb = np.repeat(q[:, None], n_samples, axis=1)
        db = draw_depths(len(q), n_samples, cfg["depth"], rng, sample_means, sample_sds)
        ab, pb = draw_alt_reads(vb, db, cfg["reads"], rng)
        v, depth, alt, p = np.vstack([v, vb]), np.vstack([depth, db]), np.vstack([alt, ab]), np.vstack([p, pb])
        is_bg = np.r_[is_bg, np.ones(len(q), dtype=bool)]
    callable_, called, status = call_variants(alt, depth, cfg["caller"])
    keep = called.any(axis=1) if bool(cfg["caller"]["retain_called_any"]) else np.ones(len(v), dtype=bool)
    observed = np.divide(alt, depth, out=np.full(alt.shape, np.nan), where=depth > 0)
    return {
        "assay_vaf": v, "depth": depth, "alt": alt, "read_probability": p, "observed_vaf": observed,
        "callable": callable_, "called": called, "status": status, "is_background": is_bg, "keep": keep,
        "version": OBSERVATION_MODEL_VERSION, "config_sha256": config_sha256(cfg), "config": cfg,
    }
