"""vafSOMA input layout: depth tiers as filters of one read set.

vafSOMA reads one VAF table per depth cutoff ("tier"). In real data the tiers are filters of
the same reads, and a site enters a tier only if every sample reaches the cutoff. This module
reproduces that layout from an observe() result (benchmark generator v4 conventions).
"""
from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd


def depth_tier_tables(obs: dict, samples: Sequence[str], site_ids: Sequence[str], *, first_cut: int,
                      n_tiers: int = 4, site_filter_all_samples: bool = True,
                      chrom: str = "chr1", pos_start: int = 10_000, pos_step: int = 3) -> dict[int, pd.DataFrame]:
    """Return {cut: table} with columns sample, chrom, pos, ref, alt, ref_depth, alt_depth, DP, VAF,
    site_id. Only ascertained sites (obs['keep']) are written; with site_filter_all_samples a site
    is in a tier only if DP >= cut in all samples."""
    dp, alt, keep = obs["depth"], obs["alt"], obs["keep"]
    n_sites, n_samples = dp.shape
    if len(samples) != n_samples or len(site_ids) != n_sites:
        raise ValueError("samples / site_ids do not match the observation matrices")
    pos = pos_start + pos_step * np.arange(n_sites)
    out = {}
    for t in range(int(n_tiers)):
        cut = int(first_cut) + t
        ok = keep & ((dp >= cut).all(axis=1) if site_filter_all_samples else np.ones(n_sites, bool))
        idx = np.flatnonzero(ok)
        d = pd.DataFrame({
            "sample": np.tile(np.asarray(samples, dtype=object), len(idx)),
            "chrom": chrom, "pos": np.repeat(pos[idx], n_samples), "ref": "A", "alt": "T",
            "ref_depth": (dp - alt)[idx].ravel(), "alt_depth": alt[idx].ravel(), "DP": dp[idx].ravel(),
            "site_id": np.repeat(np.asarray(site_ids, dtype=object)[idx], n_samples),
        })
        if not site_filter_all_samples:
            d = d[d.DP >= cut]
        d["VAF"] = np.where(d.DP > 0, d.alt_depth / d.DP.where(d.DP > 0, 1), 0.0)
        out[cut] = d.reset_index(drop=True)
    return out
