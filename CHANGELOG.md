# Changelog

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
