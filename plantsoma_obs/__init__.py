"""plantsoma_obs: shared, versioned sequencing observation model of the plantSOMA tools.

One implementation of how true (latent) somatic VAFs become read counts, used by simSOMA
(forward simulation) and intended to be pinned by vafSOMA and fitSOMA so all tools generate
and interpret data under the identical observation model.

Every output records ``OBSERVATION_MODEL_VERSION`` and ``config_sha256(config)``. Any change
that alters generated data for a given configuration (including a default) requires a new
version and a CHANGELOG entry.

Version history
---------------
1.0.0 (2026-10-08)  First packaged version. Read model = union of
    * simSOMA mutation-level transform (fixed / normal per-sample depth, binomial or
      beta-binomial reads, symmetric sequencing error, caller thresholds, called-in-any-sample
      ascertainment), unchanged numerically for existing configurations, and
    * vafSOMA benchmark generator v3 (simSOMA_extensions, 2026-10-07): log-normal site x
      sample depth effects (sd 0.58 / 0.235 on the log scale, fitted to beech candidate
      sites; sites > 3x median depth removed), shared background artefact sites
      (exponential / gamma / log-normal / mixture VAF distributions; apricot layer data:
      gamma, shape ~2, mean ~0.023), depth tiers as filters of one read set, and the
      all-samples site-depth filter of vafSOMA input tables.
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
