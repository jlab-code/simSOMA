"""Empirical background pools: the non-somatic candidates of a real reference dataset, processed with the same upstream
pipeline (caller and settings), used as the background of simulated tables.

Why: what the caller reports besides somatic mutations (artefacts, error-driven sites, sub-clonal sites) depends on the
caller's filters in ways that cannot be identified from passed sites alone. Taking it from real data processed the same
way gives simulated tables the background a downstream tool (vafSOMA, fitSOMA) meets in real data. When rejected sites
are available (plantSOMA `calibrate` stage), a parametric background can be fitted and checked against the pool.

A pool holds, per site, the alt reads and depths of every sample of one real design (tumour samples of one tissue and
the matched normal(s)). Two ways to use it:
  exact    the simulated design has the same samples as the real one (same tissue, number of tumour samples and normals,
           per-sample depth means taken from the real design): real rows are inserted as they are (columns in the real
           order, which follows the tree's branches).
  profile  any other design: each site keeps its tissue-level VAF (tumour samples pooled) and the normal's VAF, drawn
           from Beta(y + 1/2, n - y + 1/2); depths and reads for the simulated samples come from the observation model.
           Relies on the background VAF being shared across organs of one tissue.
Pool rows are kept without re-applying the caller (they were called in the real design).

Which sites are excluded as somatic: sites carried as layer-fixed in the deepest design of the same tissue
(`somatic_reference_sites`: a tumour sample with >= 5 alt reads at VAF >= 0.25, no reads in the normal, the site in the
main depth range, not carried by every tumour sample). Real sub-clonal mutations stay in the pool (conservative:
slightly more background).
"""
from __future__ import annotations

import hashlib
import io
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

import numpy as np

POOL_FORMAT = "plantsoma_obs.pool/1"


def _site_ids(df) -> np.ndarray:
    return (df["chrom"].astype(str) + ":" + df["pos"].astype(str) + ":" + df["ref"].astype(str) + ">" + df["alt"].astype(str)).values


def _wide(df, samples: Sequence[str]):
    """Long table (sample, chrom, pos, ref, alt, alt_depth, DP) -> site ids, A, D (sites x samples, given order)."""
    x = df.copy(); x["site"] = _site_ids(x)
    A = x.pivot_table(index="site", columns="sample", values="alt_depth", aggfunc="first")
    D = x.pivot_table(index="site", columns="sample", values="DP", aggfunc="first")
    missing = [s for s in samples if s not in A.columns]
    if missing:
        raise ValueError(f"samples not in the table: {missing}")
    A, D = A[list(samples)], D[list(samples)]
    ok = A.notna().all(axis=1) & D.notna().all(axis=1)
    return A.index.values[ok.values], A.values[ok.values].astype(np.int64), D.values[ok.values].astype(np.int64)


def somatic_reference_sites(table, tumour_samples: Sequence[str], normal_samples: Sequence[str] = (), min_alt: int = 5,
                            min_vaf: float = 0.25, min_site_depth_factor: float = float(np.exp(-1)),
                            exclude_all_tumours: bool = True) -> set:
    """Sites taken as somatic (layer-fixed) in a (deepest) reference design, to be excluded from the pool:
    a carrier tumour sample (alt >= min_alt and VAF >= min_vaf), no alt reads in the normal(s), the site in the main depth
    range (mean of D / per-sample mean >= min_site_depth_factor; low-depth sites give high VAF from few reads), and, with
    exclude_all_tumours, not carried by every tumour sample (a fixed-like site in all tumour samples but absent from the
    normal is an artefact unless the samples share a lineage segment after the normal's divergence; set False then)."""
    cols = list(tumour_samples) + list(normal_samples)
    sites, A, D = _wide(table, cols)
    nt = len(tumour_samples)
    At, Dt = A[:, :nt], D[:, :nt]
    V = At / np.maximum(Dt, 1)
    car = (At >= min_alt) & (V >= min_vaf)
    ok = car.any(axis=1)
    if len(normal_samples):
        ok &= (A[:, nt:] == 0).all(axis=1)
    g = (D / np.maximum(D.mean(axis=0), 1e-9)).mean(axis=1)
    ok &= g >= min_site_depth_factor
    if exclude_all_tumours:
        ok &= car.sum(axis=1) < nt
    return set(sites[ok])


def build_pool(table, tumour_samples: Sequence[str], normal_samples: Sequence[str], tissue: str, label: str,
               source: str = "", exclude_sites: Optional[set] = None, meta: Optional[Mapping[str, Any]] = None,
               somatic_means: Optional[Sequence[float]] = None) -> dict:
    """Pool from the lowest-tier table of one real design (vafSOMA input format, one row per site and sample)."""
    cols = list(tumour_samples) + list(normal_samples)
    sites, A, D = _wide(table, cols)
    keep = ~np.isin(sites, list(exclude_sites or set()))
    pool = dict(format=POOL_FORMAT, label=str(label), tissue=str(tissue), source=str(source),
                tumour_samples=list(tumour_samples), normal_samples=list(normal_samples),
                sample_means=D.mean(axis=0).round(3).tolist(), n_sites=int(keep.sum()), n_excluded=int((~keep).sum()),
                somatic_means=(list(np.round(somatic_means, 3)) if somatic_means is not None else D.mean(axis=0).round(3).tolist()),
                sites=sites[keep], A=A[keep], D=D[keep], meta=dict(meta or {}))
    pool["sha256"] = pool_sha256(pool)
    return pool


def pool_sha256(pool: Mapping[str, Any]) -> str:
    h = hashlib.sha256()
    h.update(json.dumps({k: pool[k] for k in ["format", "label", "tissue", "tumour_samples", "normal_samples"]}, sort_keys=True).encode())
    h.update(json.dumps(pool.get("somatic_means")).encode())
    for k in ["A", "D"]:
        h.update(np.ascontiguousarray(pool[k], dtype=np.int64).tobytes())
    h.update("\n".join(map(str, pool["sites"])).encode())
    return h.hexdigest()


def save_pool(pool: Mapping[str, Any], path) -> None:
    meta = {k: v for k, v in pool.items() if k not in ("A", "D", "sites")}
    np.savez_compressed(path, A=pool["A"], D=pool["D"], sites=np.asarray(pool["sites"], dtype=str),
                        meta=np.frombuffer(json.dumps(meta).encode(), dtype=np.uint8))


def load_pool(path) -> dict:
    z = np.load(path, allow_pickle=False)
    meta = json.loads(bytes(z["meta"]).decode())
    pool = dict(meta, A=z["A"].astype(np.int64), D=z["D"].astype(np.int64), sites=z["sites"].astype(object))
    if pool_sha256(pool) != meta["sha256"]:
        raise ValueError(f"pool {path}: checksum mismatch")
    return pool


def load_pools(directory) -> dict:
    """All pools (*.npz) of a data pack directory, by label."""
    return {p["label"]: p for p in (load_pool(f) for f in sorted(Path(directory).glob("*.npz")))}


def edge_from_candidates(A: np.ndarray, D: np.ndarray, tumour_mask: np.ndarray, quantile: float = 0.01,
                         n_bins: int = 12, min_per_bin: int = 30) -> list:
    """Calling edge [[pooled tumour depth, min pooled alt reads], ...] of a real design: the given quantile of the pooled
    tumour alt count among its candidates, in bins of pooled tumour depth (equal-count bins), made non-decreasing.
    The candidates are what the caller passed, so the lower envelope is the caller's rule in these two quantities."""
    t = np.asarray(tumour_mask, bool)
    y, n = A[:, t].sum(1), D[:, t].sum(1)
    qs = np.unique(np.quantile(n, np.linspace(0, 1, n_bins + 1)))
    pts = []
    for lo, hi in zip(qs[:-1], qs[1:]):
        m = (n >= lo) & (n <= hi)
        if m.sum() >= min_per_bin:
            pts.append([float(np.median(n[m])), float(np.quantile(y[m], quantile))])
    ys = np.maximum.accumulate([p[1] for p in pts])
    return [[p[0], float(v)] for p, v in zip(pts, ys)]


def draw_pool_rows(pool: Mapping[str, Any], spec: Mapping[str, Any], n_samples: int, tumour_cols: np.ndarray,
                   normal_cols: np.ndarray, cfg: Mapping[str, Any], rng: np.random.Generator,
                   sample_means: Optional[Sequence[float]], sample_sds: Optional[Sequence[float]],
                   sample_factor: Optional[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    """(depth, alt) rows of one pool for the simulated design. tumour_cols: simulated tumour samples of the pool's tissue;
    normal_cols: simulated normal samples (in order)."""
    from .model import draw_depths, draw_alt_reads
    mode, scale = str(spec["mode"]), float(spec["scale"])
    nt, nn = len(tumour_cols), len(normal_cols)
    N = pool["n_sites"]
    k = N if scale == 1.0 else int(rng.poisson(scale * N))
    idx = np.arange(N) if k == N else rng.integers(0, N, k)
    A0, D0 = pool["A"][idx], pool["D"][idx]
    pt, pn = len(pool["tumour_samples"]), len(pool["normal_samples"])
    depth = np.zeros((k, n_samples), np.int64); alt = np.zeros((k, n_samples), np.int64)
    if mode == "exact":
        if nt != pt or nn != pn:
            raise ValueError(f"pool {pool['label']}: exact mode needs {pt} tumour and {pn} normal samples of type "
                             f"{pool['tissue']} in the design (got {nt} and {nn}); use mode 'profile'")
        other = _other_columns(n_samples, tumour_cols, normal_cols)
        if other.size:                     # samples of other tissues: the site is absent there (error reads only)
            Do = draw_depths(k, n_samples, cfg["depth"], rng, sample_means, sample_sds, sample_factor=sample_factor)
            Ao, _ = draw_alt_reads(np.zeros((k, n_samples)), Do, cfg["reads"], rng)
            depth[:, other], alt[:, other] = Do[:, other], Ao[:, other]
        depth[:, tumour_cols], alt[:, tumour_cols] = D0[:, :pt], A0[:, :pt]
        if nn:
            depth[:, normal_cols], alt[:, normal_cols] = D0[:, pt:], A0[:, pt:]
        return depth, alt
    if mode != "profile":
        raise ValueError(f"pool mode must be exact or profile, got {mode!r}")
    yt, ntot = A0[:, :pt].sum(1), D0[:, :pt].sum(1)
    q = yt / np.maximum(ntot, 1)                     # observed tissue-level fraction (keeps the real selection)
    V = np.zeros((k, n_samples)); V[:, tumour_cols] = q[:, None]
    if nn and pn:
        yn, nnd = A0[:, pt:].sum(1), D0[:, pt:].sum(1)
        V[:, normal_cols] = (yn / np.maximum(nnd, 1))[:, None]
    # depths: keep each site's own depth level (its site factor in the real design: mean of D / real sample means), so
    # that pool sites stay in the depth classes where the caller found them; per-sample variation from the model
    from .model import _per_sample
    means = _per_sample(sample_means, float(cfg["depth"]["mean"]), n_samples)
    hs = np.ones(n_samples) if sample_factor is None else np.asarray(sample_factor, float)
    pm = np.asarray(pool["sample_means"], float)                        # candidate-site means of the real design
    ref = np.asarray(pool.get("somatic_means", pool["sample_means"]), float)   # the means the design preset expects
    g = (D0 / np.maximum(pm, 1e-9)[None, :]).mean(axis=1)
    # the simulated per-sample means are given on the scale of `somatic_means`; background sites keep the scale of the
    # real candidates: depth = (sim mean / reference mean) * candidate mean * site level
    pool_scale = float(np.mean(pm / np.maximum(ref, 1e-9)))
    means = means * pool_scale
    se = float(cfg["depth"]["sample_sdlog"]) if str(cfg["depth"]["mode"]) == "lognormal_site_sample" else 0.0
    eps = np.exp(rng.normal(-se ** 2 / 2, se, (k, n_samples))) if se > 0 else np.ones((k, n_samples))
    Dd = np.maximum(rng.poisson(means[None, :] * hs[None, :] * g[:, None] * eps), int(cfg["depth"]["minimum"]))
    Ad, _ = draw_alt_reads(V, Dd, cfg["reads"], rng)
    return Dd.astype(np.int64), Ad.astype(np.int64)         # other tissues: V = 0 (error reads only)


def _other_columns(n_samples: int, tumour_cols, normal_cols) -> np.ndarray:
    used = set(np.asarray(tumour_cols, int).tolist()) | set(np.asarray(normal_cols, int).tolist())
    return np.array([j for j in range(n_samples) if j not in used], dtype=int)


def load_data_pack(directory) -> tuple[dict, dict]:
    """A calibration data pack (directory): pools (*.npz) and presets.json. Returns (pools by label, presets by name);
    the presets are also registered in plantsoma_obs.PRESETS so that {"preset": name} works."""
    from .model import PRESETS
    d = Path(directory)
    pools = load_pools(d)
    presets = json.loads((d / "presets.json").read_text())
    for name, p in presets.items():
        for spec in p.get("background", {}).get("pools", []):
            if spec["label"] not in pools or pools[spec["label"]]["sha256"] != spec["sha256"]:
                raise ValueError(f"preset {name}: pool {spec['label']} missing or checksum differs")
        PRESETS[name] = p
    return pools, presets
