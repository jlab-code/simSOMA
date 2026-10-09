# Topologies from terrestrial laser scans (TLS)

`simsoma topology-from-tls` converts a TreeQSM-style segment table (tab-separated, one row per
segment) into a simSOMA topology in meters. Required columns: `segment`, `parent_segment`
(0 = root), `branch_order`, `length_m`, `base_distance_m`. Implementation:
`simSOMA_corefunc/topology_tls.py` (converter version in the output, `conversion.converter_version`).

```bash
simsoma topology-from-tls simSOMA_inputs/examples/tls_synthetic_segments.txt tls_r10.json \
    --organs random:10 --min-axis-length 0.1 --seed 1 \
    --report tls_r10_report.json --axis-table tls_axes.csv
```

## Rules
1. **Axes.** Same-order segment chains are merged into one branch axis; a new axis starts at the
   root and wherever the branch order increases. If a segment has several same-order children, the
   thickest (then longest) continues the axis and the others start new axes (warning).
2. **Branch positions** use the path distance from the base:
   `pos = (base_distance(child) - base_distance(parent axis start)) / length(parent axis)`.
   The named parent segment is not used, because child base distances often lie beyond the end of
   that segment (common in real scans). Values outside [0, 1] are clipped and reported.
3. **Organs.** Candidates are the tips of all axes (on a whole crown typically 10^4 or more). `--organs` selects them:
   `all`, `random:<n>`, `min_order:<k>`, `orders:<k1,k2>`, or `@file` with tip segment IDs.
   `--min-axis-length` excludes short axes. Organ IDs are `O_<tip segment>`.
4. **Pruning** (default; `--no-prune` to keep everything). Axes without a selected organ downstream
   are dropped. Organ distributions are unchanged (branch events leave the parent niche unchanged),
   but random draws differ from an unpruned run.

The output also carries `start_age` / `end_age` (path coordinate from the base, m) for the topology
plotter, and a `conversion` block with counts and settings. `--axis-table` maps every scan segment
to its axis, for mapping results back onto the scan.

## Configuration (input mode `topology_tls`)
A segment table can be given directly in the config; it is converted at run time with the same
rules (both for `simsoma check` and `simsoma run`, also with `--splits`):

```json
"topology": {
  "topology_tls": {"segments": "simSOMA_inputs/examples/tls_synthetic_segments.txt",
                   "organs": "all", "min_axis_length": 0.0, "prune": true, "seed": 0},
  "mapping_unit": "meters", "mapping_mode": "deterministic", "mapping_rate": 5.0
}
```

* `segments` (required): path, resolved like `topology_json`. Shorthand: `"topology_tls": "<path>"`.
* `organs`, `min_axis_length`, `prune`, `seed`: as the CLI options (`--organs`, `--min-axis-length`,
  `--no-prune`, `--seed`); `organs` also accepts a list of tip segment IDs or `"@file"`.
* `topology_json` and `topology_tls` are mutually exclusive; unknown keys are rejected;
  `mapping_unit` defaults to (and must be) `meters`.
* Output: `<outdir>/<experiment>/topology_input/topology_from_tls.json` and
  `tls_conversion_report.json` (counts, settings, warnings, SHA-256 of the segment table).

`mapping_rate` = self-renewal rounds per meter of axis. The per-meter rate is a modelling choice:
rounds per year divided by annual extension growth. Example: `simSOMA_configs/example_tls_tree.json`
with the synthetic table `simSOMA_inputs/examples/tls_synthetic_segments.txt` (50 segments, 17 tips).

## Path distances
Branch positions and the plotter coordinates (`start_age`, `end_age`) use TreeQSM's
`base_distance_m`; the root-to-tip path of a tip is the base distance of its axis plus the axis length.
In TreeQSM output a child's base distance usually differs slightly from its parent segment's end
(median a few cm), so summing segment lengths along the parent chain (TreeQSM `twig_distance_m`)
gives 4-9% longer paths on real crowns. Total branch length is identical under both conventions.

## Provenance
Based on a conversion script written with ChatGPT for earlier TLS runs (2026); reviewed and
extended for simSOMA 0.2.0 (identical axes and positions on a test crown; added organ selection, pruning,
several-same-order-children handling, plotter coordinates, tests).
