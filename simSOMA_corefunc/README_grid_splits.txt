# simSOMA split-grid launcher

`launch_grid_splits.py` takes one grouped `grid_parameter` master config, splits the full Cartesian parameter grid into contiguous subsets, writes child configs, launches them, and merges the outputs back into the master run directory.

## Usage

```bash
cd simSOMA_corefunc
python launch_grid_splits.py \
  --master-config ../../simSOMA_configs/master_config_template.json \
  --n-splits 8 \
  --jobs 8
```

## Notes

- `--n-splits 1 --jobs 1` is the limiting non-parallel case.
- Child configs inherit absolute `outdir_root` and `topology_json` paths resolved from the master config, so split runs remain stable regardless of where the child config files are written.
- The script reuses an existing master `topology_check/` directory when present.
- Child subconfigs are stored under `<master_run>/grid_splits/subconfigs/`.
- Merged outputs are written to the master run directory under `grid_parameter/`.
- Full-result JSONs are merged when `store_full_results=true`.
