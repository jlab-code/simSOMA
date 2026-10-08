# Observation-model update notes

This is a candidate source update based on the authoritative working copy supplied on
2026-07-24. The `VERSION` file remains `simSOMA-v0.1.0-beta`; no official release number
has been assigned by this update.

The source preserves the working copy's uncommitted observation-transform and
founder-diversity additions. It adds config-driven orchestration, fitSOMA handoff exports,
examples, documentation, tests, and validation reports.

## Run the tests

```bash
python -m unittest -v tests/test_observation_model.py
python -m compileall -q simSOMA_corefunc
```

## Run the deterministic identity example

```bash
bash simSOMA_scripts/run_config_local.sh \
  simSOMA_configs/example_observation_deterministic.json
```

## Run the read-count example

```bash
bash simSOMA_scripts/run_config_local.sh \
  simSOMA_configs/example_observation_read_counts.json
```

The standard topology check and confirmation behavior remains active. The resulting
fitSOMA datasets are placed under each scenario's `fitSOMA/` directory.

See `simSOMA_docs/observation_model_config.md`. (The validation report referenced in earlier
versions of this note was not part of the repository; the corresponding checks are in
`tests/test_observation_model.py` and `tests/test_observation_layers.py`.)

Since 0.2.0 the depth / read / caller steps are delegated to the shared `plantsoma_obs` model
(numerically identical for existing configurations).
