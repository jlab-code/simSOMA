# Changelog

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
