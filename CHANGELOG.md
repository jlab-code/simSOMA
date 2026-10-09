# Changelog

## 0.2.1 — 2026-10-09 (usability fixes found while writing the new tutorial)

- `simsoma topology-from-csv` works (the CLI passed positional arguments to a converter that
  expects `--branches/--organs/--out`, so the documented command always failed).
- Topology plot: start/end coordinates are derived from branch lengths and event positions
  when a hand-written topology omits them (the plot was empty); rotated tip labels and the
  title are no longer cut off.
- A missing `topology_json` gives a one-line error naming the resolved path (paths are read
  relative to the config file) instead of a traceback.
- `simsoma layers`: topology and output paths follow the same rules as `simsoma run`.
- `--splits`: config-relative paths (topology, output folder) now work; split subconfigs
  previously resolved them relative to `grid_splits/subconfigs/` and failed.
- 48 tests. No change to simulation results.

## 0.2.0 — 2026-10-09 (release accompanying the revised paper)

- **TLS segment tables as a topology input mode.** `topology.topology_tls` in a config converts a
  TreeQSM-style segment table at run time (same rules as `simsoma topology-from-tls`; organ
  selection, minimum axis length, pruning, seed). Works for `check`, `run` and `run --splits`
  (converted once, split runs read the converted JSON). Converted topology and report (with the
  segment table's SHA-256) go to `<outdir>/<experiment>/topology_input/`. Example config uses a
  synthetic 50-segment table (`simSOMA_inputs/examples/tls_synthetic_segments.txt`). 43 tests.
- **TLS converter** `simsoma topology-from-tls` (topology_tls 1.0.0; `simSOMA_docs/tls_topology.md`).

### Revision code fixes (branch `revision/code-fixes`)

Fixes agreed 2026-10-07 during the manuscript revision. Paper-relevant notes in brackets.

- **Removed stochastic (Poisson) topology mapping.** Only `mapping_mode: deterministic`
  (`T_steps = round(kappa_sr * L)`) is supported; `poisson` now fails with a clear error.
  The removed sampler (Knuth algorithm) saturated at ~745 steps per branch because
  `exp(-mean)` underflows, and it was drawn once per run (topology-check step, fixed seed)
  rather than per replicate. [Paper runs used deterministic mapping: results unaffected.
  SI Sec. 2.5 and Appendix Table 1 (`mapping mode`) and the tutorial must drop the Poisson option.]
- **Displacement now uses an independent daughter.** In a displacement event the displacer
  divides once more and its second daughter, with its own Poisson mutation draw, replaces the
  victim. Previously the displacer's first daughter was copied, so mutations of that round
  started at two niche positions; this was inconsistent with the amplification modules.
  Changes outputs for rho > 0 (RNG stream). Check on topology 02 (m in {2,4}, rho in
  {0.1,0.5,1}, 40 reps): fixed/intermediate/private fractions changed by < 0.005, within
  Monte-Carlo error. [SI Sec. 4.5 wording to update.]
- **Naming.** `mu_unit` is the canonical name of the mutation rate per lineage per topology
  unit (config `simulation.modules.self_renewal.mu_unit`, `SelfRenewalParams.mu_unit`,
  worker `--mu_unit`). `mu_year` remains a deprecated alias everywhere (input) and is still
  written as an extra output column. Realized/configured truth exports `P_b_eff` (branch
  precursor number realized in the child niche); `P_a_eff` is still written as a deprecated
  alias for fitSOMA <= 0.3.19. Example configs updated to `mu_unit`.
- **Founder clonal composition redefined (`local_sector_v1`).** `founder_sector_count`,
  `founder_diversity`, `founder_effective_sectors`, `founder_dominant_fraction` and
  `founder_lineage_counts` now refer to the local clonal sectors of the boundary ring at the time
  of the event (definition in `simSOMA_docs/founder_composition.md`). New: `founder_polyclonal`,
  `founder_sector_counts`, `founder_definition`; configured truth `phi_B/O`, `pi_B/O_expected`,
  `sector_count_B/O_expected` (exact expectation). The old values counted root-niche lineage labels
  and collapsed with turnover and branch order; they remain as `root_lineage_*` (comparison only).
  Field names consumed by fitSOMA are unchanged, so fitSOMA's pi_B, d_B, pi_O, d_O now follow the
  corrected definition (checked: fitSOMA targets give pi_O = 0.443 at rho = 0 and 0.5, exact 0.438).
  [Manuscript figures unaffected: none use founder/lineage statistics.]
- **Shared observation model `plantsoma_obs` (1.1.0).** Depth / read / error / caller /
  ascertainment steps of the mutation-level transform now call one versioned implementation,
  merged with the vafSOMA benchmark generator v4 (simSOMA_extensions v4; v3 withdrawn): read depth
  Poisson(D g_i e_is h_s) with log-normal site / site-by-sample / sample factors (sdlog 0.58 / 0.17 /
  0.27) fitted to apricot fruit pseudobulk input tables (Goel et al. 2024), background artefact
  sites, depth tiers as nested filters of one read set, all-sample site filter. Existing configurations are numerically
  identical (tested against the previous implementation). Outputs record model version and a
  settings digest.
- **Per-layer simulation (`layered.py`, `simsoma layers`).** Independent L1/L2/L3 lineage
  histories on one topology (own seed and mu_unit per layer; shared developmental parameters and
  phyllotactic event positions), bulk combination with global or per-organ layer contributions,
  read-level observation, vafSOMA-format depth-tier tables. Follows the conventions of the
  simSOMA_extensions v4 generator (depth draws identical for the same seed; tested).
- **Organ mutation-rate multiplier** `organ_mu_multiplier` (optional organ-module parameter,
  default 1): scales the per-division rate during organ amplification only. Historical set
  labels unchanged when 1.
- **Packaging.** `pyproject.toml` (package `simsoma` with CLI `simsoma`, modules installed as
  `simsoma_corefunc`, plus `plantsoma_obs`); file layout of `simSOMA_corefunc/` unchanged
  (fitSOMA interface). Topology plotter and version lookup work from an installed package.
  CSV topology converter moved to `simSOMA_corefunc/topology_csv.py` (old script path kept).
  VERSION -> simSOMA-v0.2.0-dev. NOTE: fitSOMA 0.3.19 accepts simSOMA >=0.1,<0.2 only; its
  adapter range needs a one-line update for 0.2.x (otherwise compatible: validate + adapter
  simulation pass with the version string patched).
- **TLS topologies** (`simsoma topology-from-tls`, `simSOMA_corefunc/topology_tls.py`): TreeQSM-style
  segment tables -> topology in meters (same-order chains merged into axes, branch positions from base
  distances, organ selection, pruning of unsampled sub-trees). Doc `simSOMA_docs/tls_topology.md`,
  example config `simSOMA_configs/example_tls_tree.json`, tests `tests/test_topology_tls.py`.
- Tests: `tests/test_observation_layers.py` (observation model, layers, organ multiplier, CLI),
  `tests/test_founder_sectors.py` (closed form; realized vs exact; independence of rho and
  branch order), `tests/test_code_fixes.py` (mapping, displacement, aliases). Checked compatible with
  fitSOMA 0.3.19.1 (`validate_simsoma` -> compatible; adapter simulation runs).

## Unreleased — config-driven observation model and fitSOMA handoff

- Added an optional top-level `observation_model` block to the ordinary simSOMA config.
- Added two modes: deterministic layer/phasing transformations and stochastic read-count observations.
- Preserved backward compatibility: configs without an observation block follow the legacy developmental pipeline.
- Added deterministic identity output for reproducing the original noiseless simSOMA VAFs.
- Added fixed-depth binomial/beta-binomial read sampling, symmetric sequencing error, callability, and caller emulation.
- Added complete mutation-by-sample fitSOMA handoffs for every parameter set and biological replicate.
- Added configured and realized truth exports, including founder-sector summaries.
- Added compact `realized_event_truth.csv.gz` output when required for handoff.
- Fixed top-level `observation_model.retain_called_any` parsing and added canonical handoff observation-contract provenance.
- Integrated observation processing with local and split runs; split jobs transform only after merged output is complete.
- Added deterministic, read-count, and split example configs plus focused tests and documentation.

## simSOMA-v0.1.0-beta — GitHub clean beta release

Initial GitHub-ready beta release for the accompanying paper:

> simSOMA: a cell-lineage based simulator of the somatic VAF spectrum in plants

Repository cleanup and release-hygiene changes:

- Added MIT `LICENSE`.
- Added `CITATION.cff` with repository citation metadata.
- Added short GitHub landing-page `README.md`.
- Added `tests/test_quick_run.sh` as a fresh-clone smoke-test script.
- Expanded `.gitignore` for virtual environments, runtime outputs, caches, build products, logs, generated documentation, and local archives.
- Removed local/generated artifacts from the release archive:
  - `.venv/`
  - `simSOMA_output/`
  - `__pycache__/`
  - generated tutorial HTML/PDF and tutorial build files from the core repository
- Kept runnable simulator code, example configs, example inputs/topologies, and scripts.
- Kept the full tutorial outside the core repository by design.

The previous internal changelog is retained in `simSOMA_docs/changelog.txt`.
