# Optional observation model and fitSOMA handoff

simSOMA retains its existing developmental simulation workflow. An optional top-level
`observation_model` block can be added to an ordinary simulation config. The block is
processed only after the normal `grid_parameter/` output has been completed and merged.
The primary developmental outputs are never modified.

## Backward compatibility

If `observation_model` is absent, simSOMA executes exactly the legacy developmental
pipeline and creates no observation-transform folder. Existing configs therefore retain
their original behavior.

There are two observation modes:

- `deterministic`: exact layer, tissue-composition, and phasing transformation without
  stochastic read noise.
- `read_counts`: the same deterministic assay transformation followed by depth/read
  sampling, sequencing error, callability, and caller emulation.

The deterministic identity configuration reproduces the original layer-equivalent
simSOMA VAFs:

```json
"observation_model": {
  "mode": "deterministic",
  "scenario_name": "legacy_noiseless",
  "layers": ["layer_equivalent"],
  "layer_weights": [1.0],
  "sampling": "layer_specific",
  "target_layer": "layer_equivalent",
  "phase": "phased",
  "export_fitsoma": true
}
```

With this configuration, the exact transformed VAF equals the original simSOMA VAF.
The optional fitSOMA handoff also stores an integer read-count encoding at
`deterministic_depth` (default 1,000,000), while preserving the unrounded value in
`exact_vaf`.

## Deterministic transformations

For a bulk sample, the exact transformed VAF is

```text
source_vaf * layer_weight * phase_factor
```

The phase factor is 1.0 for `phased` and 0.5 for `unphased`. For
`layer_specific` sampling, the target layer has contribution 1.0 and non-target layers
are omitted. `layer_weights` must still list the modeled source layers and sum to 1.

Example:

```json
"observation_model": {
  "mode": "deterministic",
  "scenario_name": "three_layer_bulk_unphased",
  "layers": ["L1", "L2", "L3"],
  "layer_weights": [0.10, 0.70, 0.20],
  "sampling": "bulk",
  "phase": "unphased",
  "export_fitsoma": true
}
```

## Read-count observations

Example:

```json
"observation_model": {
  "mode": "read_counts",
  "scenario_name": "L2_phased_WGS",
  "layers": ["L2"],
  "layer_weights": [1.0],
  "sampling": "layer_specific",
  "target_layer": "L2",
  "phase": "phased",
  "export_fitsoma": true,
  "seed": 2026072400,
  "depth": {
    "mode": "fixed",
    "value": 100
  },
  "read_sampling": {
    "distribution": "beta_binomial",
    "concentration": 200
  },
  "sequencing_error": {
    "reference_to_alternate": 0.001,
    "alternate_to_reference": 0.001
  },
  "caller": {
    "minimum_depth": 20,
    "minimum_alt_reads": 3,
    "minimum_vaf": 0.02
  }
}
```

The current read-count implementation supports fixed or sample-design depth modes,
binomial or beta-binomial read sampling, and symmetric reference/alternate sequencing
error. Asymmetric error rates are rejected rather than silently approximated.

## Output layout

Observation summaries are written under:

```text
simSOMA_output/<experiment>/observation_model_transforms/<scenario>/
```

If `export_fitsoma` is true, each parameter set and biological replicate receives a
separate handoff directory:

```text
fitSOMA/set_00001/replicate_00000/
├── variants.tsv.gz
├── sample_info.tsv
├── configured_truth.json
├── realized_truth.json
├── realized_event_truth.tsv
└── provenance.json
```

`variants.tsv.gz` is a complete mutation-by-sample evidence matrix. simSOMA does not
apply fitSOMA's confidently-ubiquitous filter. That analysis filter is applied later by
fitSOMA identically to empirical and simulated input.

Configured truth records input parameters. Realized truth also records stochastic
founder-sector outcomes, including `pi_B`, `d_B`, `pi_O`, and `d_O` where available.

## Replicates

Use the existing setting:

```json
"simulation": {
  "n_sim": 50
}
```

to generate 50 independent developmental replicates for every parameter combination.
The observation transform exports each replicate separately. For fitSOMA validation,
fit each replicate independently; do not average VAFs or read counts across replicates.

## Local and split runs

No new runner is required:

```bash
bash simSOMA_scripts/run_config_local.sh simSOMA_configs/example_observation_read_counts.json
```

For split runs, the developmental jobs run normally, their outputs are merged, and the
observation transform runs once on the merged result. It is not duplicated inside each
split job.
