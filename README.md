# simSOMA

**simSOMA** is a cell-lineage based simulator of somatic variant allele-frequency (VAF) spectra in plants. It models how shoot apical meristem dynamics, cell-lineage turnover, branch founding, organ formation and plant topology jointly shape observed VAF spectra and variant sharing across sampled organs.

This repository accompanies the paper:

> **simSOMA: a cell-lineage based simulator of the somatic VAF spectrum in plants**

## Status

Version 0.2.0 (release accompanying the revised paper). Changes since 0.1.0-beta are listed in
`CHANGELOG.md`.

## Install and command line

```bash
pip install .            # from the repository root (Python >= 3.10; numpy, pandas, matplotlib)
simsoma version
simsoma run simSOMA_configs/quick_test_2organs.json        # topology check + grid simulation
simsoma run CONFIG --splits 8 --jobs 8                     # split the parameter grid
simsoma check CONFIG                                       # topology check / plots only
simsoma layers CONFIG                                      # per-layer simulation + read-level observation
simsoma template layered > layered.json                    # template for `simsoma layers`
simsoma topology-from-csv branches.csv organs.csv topology.json --unit meters
simsoma topology-from-tls segments.txt topology.json --organs all    # TLS / TreeQSM segment table
```

## Topology input modes

The `topology` block of a config takes exactly one of:

| Input | Config key | Notes |
|---|---|---|
| Topology JSON (user or internal format) | `topology_json` | units `steps`, `years` or `meters` (`mapping_unit`) |
| Branch/organ CSV tables | `topology_json` after `simsoma topology-from-csv` | |
| TLS / TreeQSM segment table | `topology_tls` | converted at run time; meters; see below |

```json
"topology": {
  "topology_tls": {"segments": "simSOMA_inputs/examples/tls_synthetic_segments.txt",
                   "organs": "all", "min_axis_length": 0.0, "prune": true, "seed": 0},
  "mapping_unit": "meters", "mapping_rate": 5.0
}
```

`topology_tls` may also be a plain path (all tips as organs). The converted topology and the
conversion report (with the SHA-256 of the segment table) are written to
`<outdir>/<experiment>/topology_input/`. Example: `simSOMA_configs/example_tls_tree.json`
(synthetic 50-segment tree).

The shell wrappers in `simSOMA_scripts/` keep working unchanged.

## Layers and sequencing observation

* `simsoma layers` simulates the SAM layers (e.g. L1, L2, L3) as independent cell-lineage
  histories on the same topology (own seed and mutation rate per layer, shared developmental
  parameters and phyllotactic event positions), combines them into bulk organ samples with
  layer contributions, and applies the read-level observation model. It writes per-layer
  carrier fractions (truth), read evidence, and optional vafSOMA-format depth-tier tables.
  See `simSOMA_corefunc/layered.py`.
* The sequencing observation model (depth, read sampling, sequencing error, background
  artefacts, caller thresholds, ascertainment) is the separate, versioned package
  `plantsoma_obs`, shared with the other plantSOMA tools. Every output records its version and a
  SHA-256 digest of its settings. See `plantsoma_obs/__init__.py` and `plantsoma_obs/CHANGELOG.md`.

## Repository layout

```text
simSOMA/
├── simsoma/               # package entry point and command-line interface
├── plantsoma_obs/         # shared, versioned sequencing observation model
├── simSOMA_corefunc/      # simulator source code and topology plotting/checking code
├── simSOMA_scripts/       # setup, quick-test, local-run, cluster-run, and utility wrappers
├── simSOMA_configs/       # example simulation configuration files
├── simSOMA_inputs/        # example topologies, templates, and CSV-topology examples
├── simSOMA_docs/          # stable design notes and internal changelog
├── docs/                  # pointer to external tutorial documentation
├── tests/                 # smoke-test scripts
├── pyproject.toml
├── requirements.txt
├── requirements.lock.txt
├── VERSION
├── CITATION.cff
└── LICENSE
```

Runtime outputs are written to `simSOMA_output/`. This folder is intentionally ignored by Git.

## Installation

A Python 3.11 environment is recommended. From the repository root, run:

```bash
bash simSOMA_scripts/01_setup_env.sh
```

The setup script creates a project-local `.venv/` and installs the pinned dependencies from `requirements.lock.txt` when available.

## Quick test

Run the tested smoke workflow:

```bash
SIMSOMA_SKIP_TOPOLOGY_CONFIRM=1 bash simSOMA_scripts/run_quick_test.sh
```

Expected output appears under:

```text
simSOMA_output/quick_test_2organs/
simSOMA_output/topology_check/quick_test_2organs_topology/
```

To run the same check through the test wrapper:

```bash
bash tests/test_quick_run.sh
```

## Check a topology

Before running a larger simulation, plot and inspect the topology:

```bash
bash simSOMA_scripts/check_topology.sh simSOMA_inputs/examples/quick_test_2organs_topology.json
```

The topology check writes PNG/PDF plots, a layout table, and a report JSON to:

```text
simSOMA_output/topology_check/<topology_name>/
```

## Run a configuration locally

```bash
bash simSOMA_scripts/run_config_local.sh simSOMA_configs/quick_test_2organs.json
```

Optional split run:

```bash
N_SPLITS=4 N_PARALLEL_JOBS=2 RUN_CHECK=yes \
  bash simSOMA_scripts/run_config_local.sh simSOMA_configs/simulation_04_organ_formation_rough.json
```

## Optional observation model

An ordinary simSOMA config may include a top-level `observation_model` block. The
developmental simulation still runs through the standard local or split workflow; the
observation transform is applied only after the normal outputs have been completed.

Two modes are available:

- `deterministic`: exact layer/tissue/phasing transformation without stochastic read noise.
- `read_counts`: deterministic assay transformation followed by depth/read sampling,
  sequencing error, callability, and caller emulation.

If the block is absent, simSOMA retains the legacy Git behavior and creates no observation
outputs. A deterministic identity configuration (`layer_equivalent`, weight 1,
layer-specific, phased) reproduces the original noiseless VAFs while optionally exporting
a complete fitSOMA handoff.

Runnable examples:

```text
simSOMA_configs/example_observation_deterministic.json
simSOMA_configs/example_observation_read_counts.json
simSOMA_configs/example_observation_split.json
```

Run them through the existing wrapper:

```bash
bash simSOMA_scripts/run_config_local.sh simSOMA_configs/example_observation_deterministic.json
```

Detailed schema and output documentation are in
`simSOMA_docs/observation_model_config.md`.

## Full tutorial

The full tutorial is maintained outside the core code repository because it is updated frequently.

Suggested tutorial repository/site:

```text
https://github.com/jlab-code/simSOMA-tutorial
```

## Citation

If you use simSOMA, please cite the accompanying paper and the software repository. Repository citation metadata are provided in `CITATION.cff`.

## License

simSOMA is released under the MIT License. See `LICENSE`.

## Topologies from terrestrial laser scans

Segment tables can be used directly as a config input (`topology.topology_tls`, above) or converted
once with `simsoma topology-from-tls SEGMENTS.txt OUT.json --organs random:30` (topology in meters).
Conversion rules: `simSOMA_docs/tls_topology.md`.
