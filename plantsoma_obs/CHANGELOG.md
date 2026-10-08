# plantsoma_obs changelog

Any change that alters generated data for a given configuration (including defaults) gets a new
version here and in `OBSERVATION_MODEL_VERSION`.

## 1.1.0 (2026-10-08)
Read-depth model of the vafSOMA benchmark generator v4: D_is ~ Poisson(mean * g_i * e_is * h_s),
log-normal site, site-by-sample and sample (organ/library) factors with sdlog 0.58 / 0.17 / 0.27
(means 1), fitted to candidate sites of the apricot fruit pseudobulk input tables (Goel et al. 2024
data; lowest depth cutoff, Poisson variance subtracted); site factors > 3 removed. New key
`depth.sample_factor_sdlog` (0.27); `depth.sample_sdlog` default 0.235 -> 0.17. Changes data of the
`lognormal_site_sample` depth mode only.

## 1.0.0 (2026-10-08)
First packaged version: simSOMA mutation-level read transform (numerically unchanged) merged with
the benchmark generator's background artefact sites, depth tiers as nested filters of one read set,
and the all-sample site filter. Its log-normal depth defaults are withdrawn (superseded by 1.1.0).
