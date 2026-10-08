# Topologies from terrestrial laser scans (TLS)

`simsoma topology-from-tls` converts a TreeQSM-style segment table (tab-separated, one row per
segment) into a simSOMA topology in meters. Required columns: `segment`, `parent_segment`
(0 = root), `branch_order`, `length_m`, `base_distance_m`. Implementation:
`simSOMA_corefunc/topology_tls.py` (converter version in the output, `conversion.converter_version`).

```bash
simsoma topology-from-tls segment_tree10.txt simSOMA_inputs/tls/tree10_r30.json \
    --organs random:30 --min-axis-length 0.1 --seed 1 \
    --report tree10_r30_report.json --axis-table tree10_axes.csv
```

## Rules
1. **Axes.** Same-order segment chains are merged into one branch axis; a new axis starts at the
   root and wherever the branch order increases. If a segment has several same-order children, the
   thickest (then longest) continues the axis and the others start new axes (warning).
2. **Branch positions** use the path distance from the base:
   `pos = (base_distance(child) - base_distance(parent axis start)) / length(parent axis)`.
   The named parent segment is not used, because child base distances often lie beyond the end of
   that segment (about a third of branches in the tree10 example). Values outside [0, 1] would be
   clipped and reported (none in tree10).
3. **Organs.** Candidates are the tips of all axes (tree10: 13,614). `--organs` selects them:
   `all`, `random:<n>`, `min_order:<k>`, `orders:<k1,k2>`, or `@file` with tip segment IDs.
   `--min-axis-length` excludes short axes. Organ IDs are `O_<tip segment>`.
4. **Pruning** (default; `--no-prune` to keep everything). Axes without a selected organ downstream
   are dropped. Organ distributions are unchanged (branch events leave the parent niche unchanged),
   but random draws differ from an unpruned run.

The output also carries `start_age` / `end_age` (path coordinate from the base, m) for the topology
plotter, and a `conversion` block with counts and settings. `--axis-table` maps every scan segment
to its axis, for mapping results back onto the scan.

## Configuration
Use `mapping_unit: "meters"` and `mapping_rate` = self-renewal rounds per meter of axis
(`simSOMA_configs/example_tls_tree.json`; set `topology_json` to the converted file). The
per-meter rate is a modelling choice: rounds per year divided by annual extension growth.

## Provenance
Based on the conversion script written with ChatGPT for the TLS tree04 runs (2026); reviewed and
extended for simSOMA 0.2.0 (identical axes and positions on tree10; added organ selection, pruning,
several-same-order-children handling, plotter coordinates, tests).
