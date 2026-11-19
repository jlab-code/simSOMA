# plantsoma_obs changelog

Any change that alters generated data for a given configuration (including defaults) gets a new
version here and in `OBSERVATION_MODEL_VERSION`.

## 1.1.1 (2026-10-26)

Packaged on its own for the plantSOMA tools (1.0.0–1.1.0 shipped inside simSOMA, which now bundles this version).

Changes data:
- `lognormal_site_sample`: the sample (organ / library) factor h_s is drawn once per sample and shared by somatic and
  background sites. 1.1.0 drew it a second time for the background sites. Configurations without background sites, and
  background configurations in the other depth modes, are unchanged (same numbers for the same seed; tested against
  1.1.0).

New options (all off by default; the defaults give 1.1.0 behaviour apart from the fix):
- `depth.site_factor`: `lognormal` (default) or `mixture` — a main log-normal (`site_sdlog`) plus a low-depth
  component (`site_low_weight`, `site_low_shift` = difference of log-means, `site_low_sdlog`); E[g] = 1 before the
  `max_site_factor` truncation.
- `reads.overdispersion`: `read_probability` (default, 1.1.0: beta-binomial on p = v(1 − e) + (1 − v)e, so error reads
  are overdispersed too) or `allele_fraction` (beta-binomial on the true allele fraction only; error reads binomial).
- `background.components`: list of `{sample_type, n_sites, distribution, mean, in_normals}`; a component's sites carry
  one VAF q_i (drawn once per site) in every sample of that type (`observe(sample_types=...)`), are independent between
  components, and are absent from the normals if `in_normals` is false. Exclusive with `background.n_sites`.
- `background.error_sites_genome_mb`: error-only candidate sites (no mutation, no artefact) from a callable genome of
  this size (3 alternative alleles per position), sampled exactly under binomial error given the depth model and the
  pooled-LOD caller (tested against brute-force simulation). Needs `caller.mode` `pooled_lod` and binomial error reads.
  Output `site_class` = somatic | background | error.
- Matched normals: `observe(normal_samples=[...])` — simulated like the other samples (same depth and read model);
  never evidence for a call; `caller.normal_max_alt` (≥ 0) drops sites with more alt reads in a normal (−1 = off).
- `caller.mode`: `per_sample` (default; called in ≥ 1 tumour sample), `pooled_lod` (pooled tumour log10 likelihood
  ratio ≥ threshold and ≥ `min_max_alt` alt reads in one tumour sample), or `edge`. The `pooled_lod` threshold is
  `caller.min_pooled_lod`, or, if `caller.lod_depth_curve` ([[pooled depth, threshold], ...]) is given, that curve
  interpolated in log depth at the design's summed tumour depth (Σ per-sample mean × sample factor), held constant
  outside the curve. `edge` uses `caller.edge_curve` ([[pooled tumour depth, min pooled tumour alt], ...]), a calling
  edge read off the candidates of a real design (`edge_from_candidates`).
- Empirical background pools (`plantsoma_obs.empirical`): the non-somatic candidates of a real reference design
  (`build_pool`, `save_pool`, `load_pool`, `load_data_pack`, `somatic_reference_sites`, `edge_from_candidates`). Config
  `background.pools` (`{label, sha256, sample_type, mode exact | profile, scale}`), `observe(..., pools=...)`; rows get
  `site_class` 'pool' and are kept. exact = real rows as they are (simulated design = real sample set); profile = site's
  tissue-level VAF and depth level from the pool, per-sample depths and reads from the model (any design). Sites
  excluded from pools as somatic (`somatic_reference_sites`): a carrier tumour sample, no alt reads in the normal, the
  site in the main depth range, not carried by every tumour sample. Pools carry `somatic_means`, the per-sample depth
  means at the reference mutations; pass these as `sample_means`, and profile mode scales background depths from them
  back to the candidate scale.
- `subset_samples(obs, columns)`: re-ascertainment on a subset of samples with the same reads (threshold recomputed).
- `callable_fraction(config, cutoffs, rng, sample_means)`: P(every sample ≥ d) under the depth model; simulated
  G(d) = genome size × this fraction.
- Presets: `{"preset": "apricot_v1", ...overrides}`, `get_preset(name)`; sensitivity presets `apricot_v1_e1e-4` and
  `apricot_v1_e1e-3` (read error 1e-4 / 1e-3). The digest is computed on the expanded configuration. Design presets
  (pools and calling edge of one real design) are read from a calibration data pack (`load_data_pack`).

Preset `apricot_v1` (provisional calibration on the apricot data of Goel et al. 2024):
depth mixture (w_low 0.22, shift −1.45, sdlog 0.39 / 0.69, cap 2.7; site × sample 0.17; sample factor 0.26);
beta-binomial on the allele fraction, concentration 140 (ψ ≈ 0.007); error 3e-4; backgrounds fruit gamma:9 mean
0.027 (330 sites) and leaf exp mean 0.037 (1250 sites), both absent from the normals; no error-only sites; caller
`edge` with the calling edge of the full-depth fruit pseudobulk data, min_max_alt 1, normal_max_alt 1. The parametric
background is meant for full depth; at lower depths use design presets (pools).

## 1.1.0 (2026-10-08)
Read-depth model of the vafSOMA benchmark generator: D_is ~ Poisson(mean * g_i * e_is * h_s),
log-normal site, site-by-sample and sample (organ/library) factors with sdlog 0.58 / 0.17 / 0.27
(means 1), fitted to candidate sites of the apricot fruit pseudobulk input tables (Goel et al. 2024
data; lowest depth cutoff, Poisson variance subtracted); site factors > 3 removed. New key
`depth.sample_factor_sdlog` (0.27); `depth.sample_sdlog` default 0.235 -> 0.17. Changes data of the
`lognormal_site_sample` depth mode only.

## 1.0.0 (2026-10-08)
First packaged version: simSOMA mutation-level read transform (numerically unchanged) merged with
the benchmark generator's background artefact sites, depth tiers as nested filters of one read set,
and the all-sample site filter. Its log-normal depth defaults are withdrawn (superseded by 1.1.0).
