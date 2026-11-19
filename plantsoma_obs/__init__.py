"""plantsoma_obs: shared, versioned sequencing observation model of the plantSOMA tools.

One implementation of how true (latent) somatic VAFs become read counts, used by simSOMA
(forward simulation) and intended to be pinned by vafSOMA and fitSOMA so all tools generate
and interpret data under the identical observation model.

Every output records ``OBSERVATION_MODEL_VERSION`` and ``config_sha256(config)``. Any change
that alters generated data for a given configuration (including a default) requires a new
version and a CHANGELOG entry.

Version history
---------------
1.1.1 (2026-10-26)  Packaged on its own for the plantSOMA tools (simSOMA bundles a copy). Fix: in the
    lognormal_site_sample mode the sample factor h_s is drawn once per sample and shared by somatic and
    background sites (1.1.0 drew it again for the background sites). New options, all off by default:
    depth.site_factor 'mixture' (main + low-depth log-normal components); background.components (one
    background per sample type, optionally absent from normals); matched normals (observe(normal_samples=...))
    with caller.normal_max_alt; caller.mode 'pooled_lod' (pooled tumour evidence) and 'edge'; empirical
    background pools (plantsoma_obs.empirical). New preset 'apricot_v1' (provisional calibration on the apricot
    data of Goel et al. 2024). Configurations without background, and background configurations in the other
    depth modes, give the same data as 1.1.0.
1.1.0 (2026-10-08)  Read-depth model of the vafSOMA benchmark generator:
    D_is ~ Poisson(mean * g_i * e_is * h_s) with log-normal site, site-by-sample and
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
    ascertain,
    pooled_lod,
    background_vaf_matrix,
    PRESETS,
    get_preset,
    lod_threshold,
    draw_error_sites,
    callable_fraction,
    subset_samples,
)
from .vafsoma_io import depth_tier_tables  # noqa: F401
from . import empirical  # noqa: F401
from .empirical import build_pool, load_pool, load_pools, save_pool, edge_from_candidates, somatic_reference_sites, load_data_pack  # noqa: F401

__version__ = OBSERVATION_MODEL_VERSION
