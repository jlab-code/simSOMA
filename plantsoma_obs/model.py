"""Read-level observation model (see package docstring for provenance and versioning).

Data flow for one site i and one sample s
-----------------------------------------
    assay VAF   v_is   = phase_s * c_ks * f_is         (computed by the caller of observe())
    read prob.  p_is   = v_is (1 - e) + (1 - v_is) e   (symmetric sequencing error e)
    depth       D_is   ~ depth model
    alt reads   A_is   ~ Binomial(D_is, p_is)  or  Binomial(D_is, Beta(p k, (1-p) k))
    call        per sample: D_is >= min_depth  and  A_is >= min_alt_reads  and  A_is/D_is >= min_observed_vaf
    ascertain   caller.mode 'per_sample': keep site if called in >= 1 (tumour) sample
                caller.mode 'pooled_lod': keep site if the pooled tumour LOD >= min_pooled_lod and
                                          max_s A_is >= min_max_alt
                normal filter (optional): drop the site if a normal sample has > normal_max_alt alt reads

Background artefact sites (optional) have v_is = q_i, the same in every sample they occur in, and pass
through the same depth / read / caller steps. Either one component for all samples (background.n_sites,
.distribution, .mean) or several components (background.components), each restricted to one sample type
and optionally absent from the normal samples.

Depth models
------------
    fixed                 D_is = mean_s
    normal                D_is = round(Normal(mean_s, sd_s)), floored at `minimum`
                          (alias: sample_design; the historical simSOMA default)
    poisson               D_is ~ Poisson(mean_s)
    lognormal_site_sample D_is ~ Poisson(mean_s g_i e_is h_s)
                          site factor g_i: 'lognormal'  log g ~ N(-sg^2/2, sg)
                                           'mixture'    two log-normals (main, low-depth tail), E[g] = 1
                          truncated at g_i <= max_site_factor;
                          log e_is ~ N(-se^2/2, se); log h_s ~ N(-so^2/2, so), ONE h_s per sample shared by
                          somatic and background sites (1.1.1; 1.1.0 drew it twice)
"""
from __future__ import annotations

import copy
import hashlib
import json
from typing import Any, Mapping, Optional, Sequence

import numpy as np

OBSERVATION_MODEL_VERSION = "1.1.1"

DEFAULTS: dict[str, Any] = {
    "depth": {
        "mode": "normal",
        "mean": 100.0,
        "sd": 15.0,
        "minimum": 1,
        "site_sdlog": 0.58,            # site factor g_i (main component if site_factor = mixture)
        "sample_sdlog": 0.17,          # site-by-sample factor e_is
        "sample_factor_sdlog": 0.27,   # sample (organ / library) factor h_s
        "max_site_factor": 3.0,
        "site_factor": "lognormal",    # 1.1.1: lognormal | mixture
        "site_low_weight": 0.0,        # mixture: weight of the low-depth component
        "site_low_shift": 0.0,         # mixture: log-mean of the low component minus that of the main one
        "site_low_sdlog": 0.5,         # mixture: sdlog of the low component
    },
    "reads": {
        "type": "beta_binomial",
        "concentration": 200.0,
        "sequencing_error": 0.001,
        "overdispersion": "read_probability",   # 1.1.1: read_probability (1.1.0) | allele_fraction
    },
    "background": {
        "n_sites": 0,
        "distribution": "exp",
        "mean": 0.015,
        "clip": 0.45,
        "components": [],              # 1.1.1: list of {sample_type, n_sites, distribution, mean, in_normals}
        "error_sites_genome_mb": 0.0,  # 1.1.1: callable genome (Mb) scanned for error-only candidate sites; 0 = off
        "pools": [],                   # 1.1.1: empirical background pools: list of {label, sha256, sample_type, mode, scale}
    },
    "caller": {
        "min_depth": 20,
        "min_alt_reads": 3,
        "min_observed_vaf": 0.0,
        "retain_called_any": True,
        "mode": "per_sample",          # 1.1.1: per_sample | pooled_lod
        "min_pooled_lod": 7.0,         # pooled_lod: threshold (log10 likelihood ratio), used if lod_depth_curve is empty
        "lod_depth_curve": [],         # pooled_lod: [[pooled tumour depth, threshold], ...], interpolated in log depth
        "edge_curve": [],              # edge: [[pooled tumour depth, min pooled tumour alt reads], ...] (empirical edge)
        "min_max_alt": 2,              # pooled_lod: at least this many alt reads in one tumour sample
        "lod_error": -1.0,             # pooled_lod: error rate in the LOD; < 0 -> reads.sequencing_error
        "normal_max_alt": -1,          # normal filter: drop sites with > this many alt reads in a normal; < 0 = off
    },
}

COMPONENT_KEYS = {"sample_type": str, "n_sites": int, "distribution": str, "mean": float, "in_normals": bool}
COMPONENT_DEFAULTS = {"in_normals": True}

DEPTH_MODES = ("fixed", "normal", "poisson", "lognormal_site_sample")
_DEPTH_ALIASES = {"sample_design": "normal", "lognormal": "lognormal_site_sample"}
READ_TYPES = ("binomial", "beta_binomial")
OVERDISPERSION = ("read_probability", "allele_fraction")
SITE_FACTORS = ("lognormal", "mixture")
CALLER_MODES = ("per_sample", "pooled_lod", "edge")
POOL_KEYS = {"label": str, "sha256": str, "sample_type": str, "mode": str, "scale": float}
POOL_DEFAULTS = {"mode": "exact", "scale": 1.0}

# ------------------------------------------------------------------------------------------ presets
# Provisional calibration on the apricot data (Goel et al. 2024). Only the observation model is set here; per-sample
# depth means, the normal samples and the sample types belong to the design and are passed to observe().
PRESETS: dict[str, dict[str, Any]] = {
    "apricot_v1": {
        "depth": {"mode": "lognormal_site_sample", "site_factor": "mixture", "site_sdlog": 0.39,
                  "site_low_weight": 0.22, "site_low_shift": -1.45, "site_low_sdlog": 0.69, "max_site_factor": 2.7,
                  "sample_sdlog": 0.17, "sample_factor_sdlog": 0.26},
        "reads": {"type": "beta_binomial", "concentration": 140.0, "sequencing_error": 3e-4,
                  "overdispersion": "allele_fraction"},
        "background": {"clip": 0.45, "error_sites_genome_mb": 0.0, "components": [
            {"sample_type": "fruit", "n_sites": 330, "distribution": "gamma:9", "mean": 0.027, "in_normals": False},
            {"sample_type": "leaf", "n_sites": 1250, "distribution": "exp", "mean": 0.037, "in_normals": False}]},
        # calling edge of the full-depth fruit pseudobulk run (minimum pooled tumour alt reads by pooled tumour depth);
        # parametric background valid at full depth only; below it use the design presets of the data pack (pools)
        "caller": {"mode": "edge", "edge_curve": [[60.5, 4.39], [133.5, 4.39], [224.0, 5.0], [322.0, 6.0], [364.0, 7.0], [389.0, 7.4], [416.5, 7.4], [460.0, 7.4], [523.0, 7.4], [640.0, 8.4], [734.0, 9.39], [855.0, 9.78]], "min_max_alt": 1, "normal_max_alt": 1, "retain_called_any": True,
                   "lod_error": 3e-4},
    },
}
# sensitivity presets for the read error
for _e, _tag in [(1e-4, "e1e-4"), (1e-3, "e1e-3")]:
    _p = copy.deepcopy(PRESETS["apricot_v1"]); _p["reads"]["sequencing_error"] = _e
    PRESETS[f"apricot_v1_{_tag}"] = _p


def get_preset(name: str) -> dict[str, Any]:
    if name not in PRESETS:
        raise ValueError(f"unknown preset {name!r}; available: {sorted(PRESETS)}")
    return copy.deepcopy(PRESETS[name])


def _merge(base: dict, upd: Optional[Mapping]) -> dict:
    out = copy.deepcopy(base)
    for k, v in (upd or {}).items():
        if isinstance(v, Mapping) and isinstance(out.get(k), dict):
            out[k] = _merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def normalize_config(config: Optional[Mapping[str, Any]] = None) -> dict[str, Any]:
    """Fill defaults and validate. A top-level key 'preset' expands that preset first; the other keys of the
    configuration override it. Raises ValueError with the offending key."""
    config = dict(config or {})
    preset = config.pop("preset", None)
    base = _merge(DEFAULTS, get_preset(preset)) if preset else DEFAULTS
    cfg = _merge(base, config)
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
    if float(d["max_site_factor"]) <= 0 or min(float(d["site_sdlog"]), float(d["sample_sdlog"]),
                                               float(d["sample_factor_sdlog"]), float(d["site_low_sdlog"])) < 0:
        raise ValueError("depth.max_site_factor must be > 0 and sdlog values >= 0")
    if str(d["site_factor"]) not in SITE_FACTORS:
        raise ValueError(f"depth.site_factor must be one of {list(SITE_FACTORS)}, got {d['site_factor']!r}")
    if not (0.0 <= float(d["site_low_weight"]) < 1.0):
        raise ValueError("depth.site_low_weight must be in [0, 1)")
    r = cfg["reads"]
    if r["type"] not in READ_TYPES:
        raise ValueError(f"reads.type must be one of {list(READ_TYPES)}, got {r['type']!r}")
    if not (0.0 <= float(r["sequencing_error"]) < 0.5):
        raise ValueError("reads.sequencing_error must be in [0, 0.5)")
    if str(r["overdispersion"]) not in OVERDISPERSION:
        raise ValueError(f"reads.overdispersion must be one of {list(OVERDISPERSION)}, got {r['overdispersion']!r}")
    if float(r["concentration"]) <= 0:
        raise ValueError("reads.concentration must be > 0 (beta-binomial: rho = 1 / (concentration + 1))")
    b = cfg["background"]
    if int(b["n_sites"]) < 0 or float(b["mean"]) <= 0 or not (0 < float(b["clip"]) <= 1):
        raise ValueError("background.n_sites >= 0, background.mean > 0, 0 < background.clip <= 1")
    _parse_background_distribution(str(b["distribution"]), float(b["mean"]))
    if not isinstance(b["components"], (list, tuple)):
        raise ValueError("background.components must be a list")
    if b["components"] and int(b["n_sites"]) > 0:
        raise ValueError("use either background.n_sites (one component for all samples) or background.components")
    comps = []
    for j, c in enumerate(b["components"]):
        c = {**COMPONENT_DEFAULTS, **dict(c)}
        if set(c) != set(COMPONENT_KEYS):
            raise ValueError(f"background.components[{j}] needs keys {sorted(COMPONENT_KEYS)}, got {sorted(c)}")
        c = {k: COMPONENT_KEYS[k](v) for k, v in c.items()}
        if c["n_sites"] < 0 or c["mean"] <= 0:
            raise ValueError(f"background.components[{j}]: n_sites >= 0 and mean > 0")
        _parse_background_distribution(c["distribution"], c["mean"])
        comps.append(dict(sorted(c.items())))
    b["components"] = comps
    pools = []
    if not isinstance(b["pools"], (list, tuple)):
        raise ValueError("background.pools must be a list")
    for j, pspec in enumerate(b["pools"]):
        pspec = {**POOL_DEFAULTS, **dict(pspec)}
        if set(pspec) != set(POOL_KEYS):
            raise ValueError(f"background.pools[{j}] needs keys {sorted(POOL_KEYS)}, got {sorted(pspec)}")
        pspec = {k: POOL_KEYS[k](v) for k, v in pspec.items()}
        if pspec["mode"] not in ("exact", "profile") or pspec["scale"] <= 0:
            raise ValueError(f"background.pools[{j}]: mode exact | profile, scale > 0")
        pools.append(dict(sorted(pspec.items())))
    b["pools"] = pools
    if float(b["error_sites_genome_mb"]) < 0:
        raise ValueError("background.error_sites_genome_mb must be >= 0")
    c = cfg["caller"]
    if str(c["mode"]) not in CALLER_MODES:
        raise ValueError(f"caller.mode must be one of {list(CALLER_MODES)}, got {c['mode']!r}")
    curve = c["lod_depth_curve"]
    if not isinstance(curve, (list, tuple)) or any(len(pt) != 2 for pt in curve):
        raise ValueError("caller.lod_depth_curve must be a list of [pooled depth, threshold] pairs")
    curve = sorted([float(x), float(y)] for x, y in curve)
    if any(x <= 0 for x, _ in curve):
        raise ValueError("caller.lod_depth_curve depths must be > 0")
    c["lod_depth_curve"] = curve
    ec = c["edge_curve"]
    if not isinstance(ec, (list, tuple)) or any(len(pt) != 2 for pt in ec):
        raise ValueError("caller.edge_curve must be a list of [pooled depth, min alt] pairs")
    c["edge_curve"] = sorted([float(x), float(y)] for x, y in ec)
    if c["mode"] == "edge" and not c["edge_curve"]:
        raise ValueError("caller.mode 'edge' needs caller.edge_curve")
    _canonical_types(cfg)
    if float(b["error_sites_genome_mb"]) > 0:
        if c["mode"] != "pooled_lod":
            raise ValueError("background.error_sites_genome_mb needs caller.mode 'pooled_lod'")
        if r["type"] == "beta_binomial" and r["overdispersion"] == "read_probability":
            raise ValueError("error-only sites need binomial error reads: reads.type 'binomial' or "
                             "reads.overdispersion 'allele_fraction'")
    if int(c["min_depth"]) < 0 or int(c["min_alt_reads"]) < 0 or not (0 <= float(c["min_observed_vaf"]) <= 1):
        raise ValueError("caller thresholds: min_depth >= 0, min_alt_reads >= 0, 0 <= min_observed_vaf <= 1")
    if int(c["min_max_alt"]) < 0 or not (float(c["lod_error"]) < 0.5):
        raise ValueError("caller.min_max_alt >= 0 and caller.lod_error < 0.5")
    if c["mode"] == "pooled_lod" and not (0 < (float(c["lod_error"]) if c["lod_error"] >= 0 else float(r["sequencing_error"]))):
        raise ValueError("caller.mode 'pooled_lod' needs an error rate > 0 (caller.lod_error or reads.sequencing_error)")
    return cfg


_INT_KEYS = {("depth", "minimum"), ("background", "n_sites"), ("caller", "min_depth"), ("caller", "min_alt_reads"),
             ("caller", "min_max_alt"), ("caller", "normal_max_alt")}
_BOOL_KEYS = {("caller", "retain_called_any")}
_STR_KEYS = {("depth", "mode"), ("reads", "type"), ("background", "distribution"), ("depth", "site_factor"),
             ("caller", "mode"), ("reads", "overdispersion")}
_LIST_KEYS = {("background", "components"), ("caller", "lod_depth_curve"), ("background", "pools"), ("caller", "edge_curve")}


def _canonical_types(cfg: dict) -> None:
    """Cast values to canonical types so equal settings hash equally (60 == 60.0)."""
    for sec, vals in cfg.items():
        for k, v in vals.items():
            if (sec, k) in _LIST_KEYS:
                continue
            if (sec, k) in _INT_KEYS:
                vals[k] = int(v)
            elif (sec, k) in _BOOL_KEYS:
                vals[k] = bool(v)
            elif (sec, k) in _STR_KEYS:
                vals[k] = str(v)
            else:
                vals[k] = float(v)


def config_sha256(config: Mapping[str, Any]) -> str:
    """Digest of the normalized configuration (preset expanded) together with the model version."""
    payload = {"version": OBSERVATION_MODEL_VERSION, "config": normalize_config(config)}
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


# ---------------------------------------------------------------------------------------- depth
def _site_factors(n_sites: int, d: Mapping[str, Any], rng: np.random.Generator) -> np.ndarray:
    sg, cap = float(d["site_sdlog"]), float(d["max_site_factor"])
    g = np.empty(n_sites)
    filled = 0
    if str(d["site_factor"]) == "lognormal":
        while filled < n_sites:                       # truncation by resampling (max-depth filter)
            x = np.exp(rng.normal(-sg ** 2 / 2, sg, max(4 * (n_sites - filled), 16)))
            x = x[x <= cap][: n_sites - filled]
            g[filled:filled + len(x)] = x
            filled += len(x)
        return g
    w, dl, sl = float(d["site_low_weight"]), float(d["site_low_shift"]), float(d["site_low_sdlog"])
    # main log-mean m0 chosen so that the untruncated mixture has E[g] = 1
    m0 = -np.log((1 - w) * np.exp(sg ** 2 / 2) + w * np.exp(dl + sl ** 2 / 2))
    while filled < n_sites:
        k = max(4 * (n_sites - filled), 16)
        low = rng.random(k) < w
        x = np.exp(np.where(low, rng.normal(m0 + dl, sl, k), rng.normal(m0, sg, k)))
        x = x[x <= cap][: n_sites - filled]
        g[filled:filled + len(x)] = x
        filled += len(x)
    return g


def draw_depths(n_sites: int, n_samples: int, depth_cfg: Mapping[str, Any], rng: np.random.Generator,
                sample_means: Optional[Sequence[float]] = None,
                sample_sds: Optional[Sequence[float]] = None,
                sample_factor: Optional[np.ndarray] = None,
                return_sample_factor: bool = False):
    """Depth matrix (n_sites x n_samples). sample_means / sample_sds override depth.mean / depth.sd
    per sample (NaN entries fall back to the defaults). lognormal_site_sample: pass sample_factor to reuse
    the sample factors h_s of an earlier call (one h_s per sample for somatic and background sites);
    with return_sample_factor the call returns (depths, h)."""
    d = depth_cfg
    mode = _DEPTH_ALIASES.get(str(d["mode"]), str(d["mode"]))
    means = _per_sample(sample_means, float(d["mean"]), n_samples)
    sds = _per_sample(sample_sds, float(d["sd"]), n_samples)
    minimum = int(d["minimum"])
    h = None
    if mode == "lognormal_site_sample":
        # D_is ~ Poisson(mean_s * g_i * e_is * h_s); draw order g, e, h as in the vafSOMA benchmark generator v4.
        se, so = float(d["sample_sdlog"]), float(d["sample_factor_sdlog"])
        g = _site_factors(n_sites, {**DEFAULTS["depth"], **d}, rng)
        e = np.exp(rng.normal(-se ** 2 / 2, se, (n_sites, n_samples)))
        if sample_factor is not None:
            h = np.asarray(sample_factor, dtype=float)
            if h.shape != (n_samples,):
                raise ValueError(f"sample_factor must have length {n_samples}")
        else:
            h = np.exp(rng.normal(-so ** 2 / 2, so, n_samples)) if so > 0 else np.ones(n_samples)
        out = np.maximum(rng.poisson(means[None, :] * g[:, None] * e * h[None, :]), minimum).astype(int)
        return (out, h) if return_sample_factor else out
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
    return (out, None) if return_sample_factor else out


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
    """VAFs of the single all-sample background component (background.n_sites)."""
    b = background_cfg
    n = int(b["n_sites"])
    if n == 0:
        return np.zeros(0)
    q = _parse_background_distribution(str(b["distribution"]), float(b["mean"]))(rng, n)
    return np.clip(q, 0.0, float(b["clip"]))


def background_vaf_matrix(background_cfg: Mapping[str, Any], n_samples: int, rng: np.random.Generator,
                          sample_types: Optional[Sequence[str]] = None,
                          normal_mask: Optional[np.ndarray] = None) -> np.ndarray:
    """Assay-VAF matrix of all background sites (n_bg x n_samples). Single component: q_i in every sample.
    Components: q_i only in the samples of its sample_type (and not in normals if in_normals is False)."""
    b = background_cfg
    comps = b.get("components") or []
    if not comps:
        q = draw_background_vafs(b, rng)
        return np.repeat(q[:, None], n_samples, axis=1)
    if sample_types is None or len(sample_types) != n_samples:
        raise ValueError("background.components need sample_types (one label per sample) in observe()")
    types = np.asarray(sample_types, dtype=object)
    nm = np.zeros(n_samples, bool) if normal_mask is None else np.asarray(normal_mask, bool)
    blocks = []
    for c in comps:
        n = int(c["n_sites"])
        if n == 0:
            continue
        cols = types == c["sample_type"]
        if not bool(c["in_normals"]):
            cols = cols & ~nm
        if not cols.any():                     # sample type absent from this design: the component adds no sites
            continue
        q = np.clip(_parse_background_distribution(c["distribution"], float(c["mean"]))(rng, n), 0.0, float(b["clip"]))
        blocks.append(q[:, None] * cols[None, :])
    return np.vstack(blocks) if blocks else np.zeros((0, n_samples))


# --------------------------------------------------------------------------------------- reads
def draw_alt_reads(assay_vaf: np.ndarray, depth: np.ndarray, reads_cfg: Mapping[str, Any],
                   rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    """Returns (alt counts, read probability incl. sequencing error)."""
    eps = float(reads_cfg["sequencing_error"])
    p = np.clip(assay_vaf * (1.0 - eps) + (1.0 - assay_vaf) * eps, 0.0, 1.0)
    if reads_cfg["type"] == "binomial":
        return rng.binomial(depth, p), p
    k = float(reads_cfg["concentration"])
    if str(reads_cfg.get("overdispersion", "read_probability")) == "allele_fraction":
        # overdispersion of the true allele fraction only; error reads stay binomial (v = 0 -> p = e)
        v = np.clip(assay_vaf, 0.0, 1.0)
        pos = v > 0
        vs = np.zeros_like(v, dtype=float)
        sv = np.clip(v[pos], 1e-10, 1.0 - 1e-10)
        vs[pos] = rng.beta(sv * k, (1.0 - sv) * k)
        return rng.binomial(depth, np.clip(vs * (1.0 - eps) + (1.0 - vs) * eps, 0.0, 1.0)), p
    safe = np.clip(p, 1e-10, 1.0 - 1e-10)
    return rng.binomial(depth, rng.beta(safe * k, (1.0 - safe) * k)), p


def call_variants(alt: np.ndarray, depth: np.ndarray, caller_cfg: Mapping[str, Any]):
    """Per-sample calls. Returns (callable, called, status) arrays."""
    vaf = np.divide(alt, depth, out=np.zeros_like(alt, dtype=float), where=depth > 0)
    callable_ = depth >= int(caller_cfg["min_depth"])
    called = callable_ & (alt >= int(caller_cfg["min_alt_reads"])) & (vaf >= float(caller_cfg["min_observed_vaf"]))
    status = np.full(called.shape, "BELOW_CALL_THRESHOLD", dtype=object)
    status[~callable_] = "NO_COVERAGE"
    status[called] = "CALLED"
    return callable_, called, status


def pooled_lod(alt: np.ndarray, depth: np.ndarray, error: float) -> np.ndarray:
    """log10 likelihood ratio, per site, of the pooled tumour reads: binomial at the observed fraction vs the
    error rate (a proxy of Mutect2's tumour LOD on joint evidence; Mutect2 also uses base qualities)."""
    y = np.asarray(alt, float).sum(axis=1); n = np.asarray(depth, float).sum(axis=1)
    f = np.clip(np.divide(y, n, out=np.zeros_like(y), where=n > 0), error, 1 - 1e-12)
    out = y * np.log10(f / error) + (n - y) * np.log10((1 - f) / (1 - error))
    return np.where(n > 0, out, 0.0)


def lod_threshold(caller_cfg: Mapping[str, Any], pooled_depth: float) -> float:
    """Pooled-LOD threshold for a design whose tumour samples have this summed mean depth: caller.lod_depth_curve
    interpolated linearly in log depth and held constant outside its range, else caller.min_pooled_lod."""
    curve = caller_cfg.get("lod_depth_curve") or []
    if not curve:
        return float(caller_cfg["min_pooled_lod"])
    x = np.log([pt[0] for pt in curve]); y = np.array([pt[1] for pt in curve], float)
    return float(np.interp(np.log(max(float(pooled_depth), 1e-9)), x, y))


def design_pooled_depth(sample_means: np.ndarray, sample_factor: Optional[np.ndarray], normal_mask: np.ndarray) -> float:
    """Summed expected depth of the tumour samples (per-sample mean x sample factor), the x of lod_depth_curve."""
    m = np.asarray(sample_means, float) * (1.0 if sample_factor is None else np.asarray(sample_factor, float))
    return float(m[~np.asarray(normal_mask, bool)].sum())


def ascertain(alt: np.ndarray, depth: np.ndarray, called: np.ndarray, caller_cfg: Mapping[str, Any],
              reads_cfg: Mapping[str, Any], normal_mask: Optional[np.ndarray] = None,
              lod_thr: Optional[float] = None) -> np.ndarray:
    """Site-level ascertainment mask (see module docstring). Normals never count as evidence for a call.
    lod_thr: pooled-LOD threshold of the design (default: caller.min_pooled_lod)."""
    n_samples = alt.shape[1]
    nm = np.zeros(n_samples, bool) if normal_mask is None else np.asarray(normal_mask, bool)
    tum = ~nm
    c = caller_cfg
    if c["mode"] == "pooled_lod":
        err = float(c["lod_error"]) if float(c["lod_error"]) >= 0 else float(reads_cfg["sequencing_error"])
        thr = float(c["min_pooled_lod"]) if lod_thr is None else float(lod_thr)
        keep = (pooled_lod(alt[:, tum], depth[:, tum], err) >= thr) & \
               (alt[:, tum].max(axis=1, initial=0) >= int(c["min_max_alt"]))
    elif c["mode"] == "edge":
        y, n = alt[:, tum].sum(axis=1), depth[:, tum].sum(axis=1)
        ec = np.asarray(c["edge_curve"], float)
        keep = (y >= np.ceil(np.interp(n, ec[:, 0], ec[:, 1]) - 1e-9)) & \
               (alt[:, tum].max(axis=1, initial=0) >= int(c["min_max_alt"]))
    elif bool(c["retain_called_any"]):
        keep = called[:, tum].any(axis=1)
    else:
        keep = np.ones(alt.shape[0], dtype=bool)
    if int(c["normal_max_alt"]) >= 0 and nm.any():
        keep &= (alt[:, nm] <= int(c["normal_max_alt"])).all(axis=1)
    return keep


def draw_error_sites(cfg: Mapping[str, Any], n_samples: int, rng: np.random.Generator, sample_means: np.ndarray,
                     sample_sds: Optional[Sequence[float]], sample_factor: Optional[np.ndarray], normal_mask: np.ndarray,
                     lod_thr: float, n_depth_vectors: int = 20000) -> tuple[np.ndarray, np.ndarray]:
    """Error-only candidate sites: positions without mutation or artefact whose error reads pass the pooled-LOD caller.
    The genome has 3 x background.error_sites_genome_mb x 1e6 (position, allele) pairs; error reads are binomial
    with rate e, so the pooled tumour count is Bin(sum D, e) and, given the total, the split over tumour samples is
    multivariate hypergeometric. Sampled exactly by conditioning on the smallest total that can pass:
    K ~ Poisson(N * mean_i P_i), depth vectors i chosen with probability ~ P_i, totals from the truncated binomial.
    Returns (depth, alt) of the sampled sites, which then go through the caller like all other sites."""
    b, c, r = cfg["background"], cfg["caller"], cfg["reads"]
    N = 3.0 * float(b["error_sites_genome_mb"]) * 1e6
    e = float(r["sequencing_error"])
    if N <= 0 or e <= 0:
        return np.zeros((0, n_samples), int), np.zeros((0, n_samples), int)
    err = float(c["lod_error"]) if float(c["lod_error"]) >= 0 else e
    D = draw_depths(n_depth_vectors, n_samples, cfg["depth"], rng, sample_means, sample_sds, sample_factor=sample_factor)
    tum = ~np.asarray(normal_mask, bool)
    n = D[:, tum].sum(1)
    # smallest pooled count passing the LOD (LOD increases with the count for y/n > err)
    ymin = np.full(len(n), np.iinfo(np.int64).max, dtype=np.int64)
    lo = max(int(c["min_max_alt"]), 1)
    for y in range(lo, 200):
        okk = (ymin == np.iinfo(np.int64).max) & (pooled_lod(np.full((len(n), 1), y), n[:, None], err) >= lod_thr)
        ymin[okk] = y
        if (ymin < np.iinfo(np.int64).max).all():
            break
    P = np.where(ymin < np.iinfo(np.int64).max, _binom_sf(ymin - 1, n, e), 0.0)
    K = rng.poisson(N * P.mean())
    if K == 0 or P.sum() == 0:
        return np.zeros((0, n_samples), int), np.zeros((0, n_samples), int)
    idx = rng.choice(len(n), size=K, p=P / P.sum())
    dep = D[idx]
    alt = np.zeros_like(dep)
    for j, i in enumerate(idx):
        T = _binom_trunc(rng, int(n[i]), e, int(ymin[i]))
        alt[j, tum] = rng.multivariate_hypergeometric(D[i, tum], T)
    if (~tum).any():
        alt[:, ~tum] = rng.binomial(dep[:, ~tum], e)
    return dep, alt


def _binom_sf(k: np.ndarray, n: np.ndarray, p: float) -> np.ndarray:
    """P(Bin(n, p) > k), vectorised (regularised incomplete beta)."""
    from scipy.special import betainc  # scipy only needed for error-only sites
    k = np.asarray(k, float); n = np.asarray(n, float)
    out = np.where(k < 0, 1.0, 0.0)
    ok = (k >= 0) & (k < n)
    out[ok] = betainc(k[ok] + 1, n[ok] - k[ok], p)
    return out


def _binom_trunc(rng: np.random.Generator, n: int, p: float, kmin: int) -> int:
    """One draw from Bin(n, p) conditioned on >= kmin (inverse CDF on the tail)."""
    from scipy.stats import binom
    ks = np.arange(kmin, min(n, kmin + 200) + 1)
    w = binom.pmf(ks, n, p)
    if w.sum() <= 0:
        return kmin
    return int(rng.choice(ks, p=w / w.sum()))


def callable_fraction(config: Optional[Mapping[str, Any]], cutoffs: Sequence[int], rng: np.random.Generator,
                      sample_means: Optional[Sequence[float]] = None, sample_sds: Optional[Sequence[float]] = None,
                      n_positions: int = 200000) -> dict[int, float]:
    """Simulated callable fraction per depth cutoff d: P(every sample has depth >= d) at a random position under the
    configured depth model (one sample factor per sample, as for the sites). Simulated G(d) = genome size x this
    fraction. This is the definition to match with the real G (computed from the reads Mutect2 counts)."""
    cfg = normalize_config(config)
    n_samples = len(sample_means) if sample_means is not None else 1
    D = draw_depths(n_positions, n_samples, cfg["depth"], rng, sample_means, sample_sds)
    m = D.min(1)
    return {int(d): float(np.mean(m >= int(d))) for d in cutoffs}


# ------------------------------------------------------------------------------------- observe
def observe(assay_vaf: np.ndarray, config: Optional[Mapping[str, Any]], rng: np.random.Generator,
            sample_means: Optional[Sequence[float]] = None,
            sample_sds: Optional[Sequence[float]] = None,
            sample_types: Optional[Sequence[str]] = None,
            normal_samples: Optional[Sequence[int]] = None,
            pools: Optional[Mapping[str, Mapping[str, Any]]] = None) -> dict[str, Any]:
    """Apply the observation model to an assay-VAF matrix (somatic sites x samples).

    pools: empirical background pools by label (plantsoma_obs.empirical.load_pools), as referenced by
    background.pools (label and sha256 must match); their rows are appended with site_class 'pool' and kept.

    sample_types: one label per sample (needed by background.components); normal_samples: indices of the
    matched-normal samples (simulated like the other samples; not used as evidence for calls; normal filter
    caller.normal_max_alt).

    Background artefact sites and error-only sites, if configured, are appended after the somatic rows. Returns a
    dict with numpy arrays (all n_total x n_samples unless noted): assay_vaf, depth, alt, read_probability,
    observed_vaf, callable, called, status; plus site_class (n_total,: 'somatic' | 'background' | 'error'),
    is_background (artefact sites), keep (n_total,, the ascertainment mask), normal_mask, sample_factor,
    pooled_depth and lod_threshold (pooled_lod caller), version and config_sha256.

    Random-number order: depths of somatic sites, reads of somatic sites, then background VAFs, depths and reads,
    then error-only sites. Configurations without background reproduce 1.1.0 (and the historical simSOMA
    transform) exactly; with background in the lognormal_site_sample mode the sample factor is no longer redrawn
    for the background sites (1.1.1).
    """
    cfg = normalize_config(config)
    v = np.clip(np.asarray(assay_vaf, dtype=float), 0.0, 1.0)
    n_sites, n_samples = v.shape
    nm = np.zeros(n_samples, bool)
    if normal_samples is not None:
        idx = np.asarray(list(normal_samples), dtype=int)
        if idx.size and (idx.min() < 0 or idx.max() >= n_samples):
            raise ValueError("normal_samples must be sample indices")
        nm[idx] = True
    means = _per_sample(sample_means, float(cfg["depth"]["mean"]), n_samples)
    depth, h = draw_depths(n_sites, n_samples, cfg["depth"], rng, sample_means, sample_sds, return_sample_factor=True)
    alt, p = draw_alt_reads(v, depth, cfg["reads"], rng)
    cls = np.full(n_sites, "somatic", dtype=object)
    vb = background_vaf_matrix(cfg["background"], n_samples, rng, sample_types, nm)
    if len(vb):
        db = draw_depths(len(vb), n_samples, cfg["depth"], rng, sample_means, sample_sds, sample_factor=h)
        ab, pb = draw_alt_reads(vb, db, cfg["reads"], rng)
        v, depth, alt, p = np.vstack([v, vb]), np.vstack([depth, db]), np.vstack([alt, ab]), np.vstack([p, pb])
        cls = np.r_[cls, np.full(len(vb), "background", dtype=object)]
    pooled = design_pooled_depth(means, h, nm)
    thr = lod_threshold(cfg["caller"], pooled) if cfg["caller"]["mode"] == "pooled_lod" else float("nan")
    if float(cfg["background"]["error_sites_genome_mb"]) > 0:
        de, ae = draw_error_sites(cfg, n_samples, rng, sample_means, sample_sds, h, nm, thr)
        if len(de):
            e = float(cfg["reads"]["sequencing_error"])
            v, depth, alt = np.vstack([v, np.zeros(de.shape)]), np.vstack([depth, de]), np.vstack([alt, ae])
            p = np.vstack([p, np.full(de.shape, e)])
            cls = np.r_[cls, np.full(len(de), "error", dtype=object)]
    n_called_rows = len(cls)
    for pspec in cfg["background"]["pools"]:
        from .empirical import draw_pool_rows
        pool = (pools or {}).get(pspec["label"])
        if pool is None:
            raise ValueError(f"background pool {pspec['label']!r} not passed to observe(pools=...)")
        if pool["sha256"] != pspec["sha256"]:
            raise ValueError(f"background pool {pspec['label']!r}: sha256 differs from the configuration")
        if sample_types is None:
            raise ValueError("background.pools need sample_types in observe()")
        types = np.asarray(sample_types, dtype=object)
        tcols = np.where((types == pspec["sample_type"]) & ~nm)[0]
        ncols = np.where(nm)[0]
        dp_, ap_ = draw_pool_rows(pool, pspec, n_samples, tcols, ncols, cfg, rng, sample_means, sample_sds, h)
        if len(dp_):
            e = float(cfg["reads"]["sequencing_error"])
            v, depth, alt = np.vstack([v, np.zeros(dp_.shape)]), np.vstack([depth, dp_]), np.vstack([alt, ap_])
            p = np.vstack([p, np.full(dp_.shape, e)])
            cls = np.r_[cls, np.full(len(dp_), "pool", dtype=object)]
    callable_, called, status = call_variants(alt, depth, cfg["caller"])
    keep = ascertain(alt, depth, called, cfg["caller"], cfg["reads"], nm, None if np.isnan(thr) else thr)
    keep[n_called_rows:] = True                    # pool rows were called in the real design
    observed = np.divide(alt, depth, out=np.full(alt.shape, np.nan), where=depth > 0)
    return {
        "assay_vaf": v, "depth": depth, "alt": alt, "read_probability": p, "observed_vaf": observed,
        "callable": callable_, "called": called, "status": status, "site_class": cls,
        "is_background": cls == "background", "keep": keep, "normal_mask": nm,
        "sample_means": means, "sample_factor": h, "pooled_depth": pooled, "lod_threshold": thr,
        "version": OBSERVATION_MODEL_VERSION, "config_sha256": config_sha256(cfg), "config": cfg,
    }


def subset_samples(obs: Mapping[str, Any], columns: Sequence[int]) -> dict[str, Any]:
    """Re-ascertain an observe() result on a subset of samples (e.g. 1, 2 or 3 leaves per tip), keeping the same reads,
    as a Mutect2 run on those samples would. Normals among the columns stay normals. The pooled-LOD threshold is
    recomputed for the subset's depth. Limitation: error-only sites were sampled for the full design; for a subset
    with a lower threshold some error-only sites are missing (negligible at full apricot depth)."""
    cols = np.asarray(list(columns), dtype=int)
    cfg = obs["config"]
    out = {k: (v[:, cols] if isinstance(v, np.ndarray) and v.ndim == 2 else v) for k, v in obs.items()}
    nm = np.asarray(obs["normal_mask"], bool)[cols]
    h = None if obs.get("sample_factor") is None else np.asarray(obs["sample_factor"])[cols]
    means = np.asarray(obs["sample_means"])[cols]
    pooled = design_pooled_depth(means, h, nm)
    thr = lod_threshold(cfg["caller"], pooled) if cfg["caller"]["mode"] == "pooled_lod" else float("nan")
    keep = ascertain(out["alt"], out["depth"], out["called"], cfg["caller"], cfg["reads"], nm,
                     None if np.isnan(thr) else thr)
    keep |= np.asarray(obs.get("site_class", np.array([], dtype=object))) == "pool"   # pool rows stay (called in the real design)
    out.update(normal_mask=nm, sample_factor=h, sample_means=means, pooled_depth=pooled, lod_threshold=thr, keep=keep,
               subset_columns=cols)
    return out
