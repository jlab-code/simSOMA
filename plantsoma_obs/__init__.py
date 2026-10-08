"""plantsoma_obs: shared, versioned sequencing observation model of the plantSOMA tools.

One implementation of how true (latent) somatic VAFs become read counts, used by simSOMA
(forward simulation) and intended to be pinned by vafSOMA and fitSOMA so all tools generate
and interpret data under the identical observation model.

Every output records ``OBSERVATION_MODEL_VERSION`` and ``config_sha256(config)``. Any change
that alters generated data for a given configuration (including a default) requires a new
version and a CHANGELOG entry.

Version history
---------------
1.1.0 (2026-10-08)  Read-depth model of the vafSOMA benchmark generator v4 (simSOMA_extensions
    v4): D_is ~ Poisson(mean * g_i * e_is * h_s) with log-normal site, site-by-sample and
    sample (organ/library) factors, sdlog 0.58 / 0.17 / 0.27, means 1, fitted to candidate sites of
    the apricot fruit pseudobulk input tables (Goel et al. 2024 data; lowest depth cutoff, Poisson
    variance subtracted); site factors > 3 removed (max-depth filter). New key
    depth.sample_factor_sdlog; default of depth.sample_sdlog 0.235 -> 0.17. Data from the
    lognormal_site_sample mode differ from 1.0.0; other modes unchanged.
1.0.0 (2026-10-08)  First packaged version: simSOMA mutation-level transform (fixed / normal
    per-sample depth, binomial or beta-binomial reads, symmetric sequencing error, caller
    thresholds, called-in-any-sample ascertainment) merged with the generator's background
    artefact sites (exponential / gamma / log-normal / mixture VAF distributions; apricot layer
    data: gamma, shape ~2, mean ~0.023), depth tiers as nested filters of one read set, and the
    all-samples site-depth filter of vafSOMA input tables. (Its log-normal depth defaults are
    withdrawn.)
"""
from .model import (  # noqa: F401
    OBSERVATION_MODEL_VERSION,
    DEFAULTS,
    normalize_config,
    config_sha256,
    draw_depths,
    draw_background_vafs,
    draw_alt_reads,
    call_variants,
    observe,
)
from .vafsoma_io import depth_tier_tables  # noqa: F401

__version__ = OBSERVATION_MODEL_VERSION
