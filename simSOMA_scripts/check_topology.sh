#!/usr/bin/env bash
set -euo pipefail

# Check and plot a simSOMA topology.
#
# Preferred usage with a topology JSON:
#   TOPO=simSOMA_inputs/examples/quick_test_2organs_topology.json
#   bash simSOMA_scripts/check_topology.sh "$TOPO"
#
# Also allowed: pass a simulation config. In that case the script reads the
# topology_json entry from the config, but the output folder is still named by
# the topology JSON filename.
#   bash simSOMA_scripts/check_topology.sh simSOMA_configs/quick_test_2organs.json
#
# Default output:
#   simSOMA_output/topology_check/<topology_json_stem>/
#
# Existing files with the same names are overwritten. This is intentional: the
# folder is a live check of the current topology file.
#
# Optional override:
#   OUTDIR=simSOMA_output/my_topology_check bash simSOMA_scripts/check_topology.sh <config-or-topology.json>

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"

INPUT_JSON="${1:-simSOMA_inputs/examples/quick_test_2organs_topology.json}"
OUT_NAME="${OUT_NAME:-topology_plot}"
PY="${SIMSOMA_PYTHON:-$PROJECT_ROOT/.venv/bin/python}"

if [[ ! -x "$PY" ]]; then
  echo "Local .venv not found. Building it with simSOMA_scripts/01_setup_env.sh." >&2
  echo "The setup script will use PYTHON_BIN if set; otherwise it will search for Python 3.11." >&2
  if [[ -n "${PYTHON_BIN:-}" ]]; then
    PYTHON_BIN="$PYTHON_BIN" bash simSOMA_scripts/01_setup_env.sh
  else
    bash simSOMA_scripts/01_setup_env.sh
  fi
  PY="$PROJECT_ROOT/.venv/bin/python"
fi

# Resolve whether INPUT_JSON is a simulation config or a raw topology JSON.
# The Python block prints two tab-separated fields: topology_json and default_outdir.
RESOLVED="$($PY - "$PROJECT_ROOT" "$INPUT_JSON" <<'PY'
import json
import re
import sys
from pathlib import Path

project_root = Path(sys.argv[1]).resolve()
input_path = Path(sys.argv[2])
if not input_path.is_absolute():
    input_path = (project_root / input_path).resolve()

if not input_path.exists():
    raise SystemExit(f"Input JSON not found: {input_path}")
if input_path.is_dir():
    raise SystemExit(
        "Input is a directory, not a JSON file: "
        f"{input_path}\n"
        "Pass the complete .json path, for example:\n"
        "  TOPO=simSOMA_inputs/examples/quick_test_2organs_topology.json\n"
        "  bash simSOMA_scripts/check_topology.sh \"$TOPO\""
    )

with input_path.open() as f:
    data = json.load(f)

# Case 1: simulation config with run/topology blocks. Use it only to find the
# topology JSON; output is still keyed by topology filename, not experiment name.
if isinstance(data, dict) and "topology" in data and "topology_json" in data.get("topology", {}):
    topo_path = Path(data["topology"]["topology_json"])
    if not topo_path.is_absolute():
        # Config paths are interpreted relative to the config file location first.
        candidate = (input_path.parent / topo_path).resolve()
        if candidate.exists():
            topo_path = candidate
        else:
            topo_path = (project_root / topo_path).resolve()
else:
    # Case 2: raw topology JSON.
    topo_path = input_path

if not topo_path.exists():
    raise SystemExit(f"Topology JSON not found: {topo_path}")

stem = re.sub(r"[^A-Za-z0-9_.-]+", "_", topo_path.stem).strip("_")
default_outdir = project_root / "simSOMA_output" / "topology_check" / stem
print(f"{topo_path}\t{default_outdir}")
PY
)"

TOPOLOGY_JSON="$(printf '%s' "$RESOLVED" | cut -f1)"
DEFAULT_OUTDIR="$(printf '%s' "$RESOLVED" | cut -f2)"
OUTDIR="${OUTDIR:-$DEFAULT_OUTDIR}"

"$PY" simSOMA_corefunc/plot_topology_json.py \
  --topology_json "$TOPOLOGY_JSON" \
  --outdir "$OUTDIR" \
  --out_name "$OUT_NAME"

echo
echo "Topology input:"
echo "  $TOPOLOGY_JSON"
echo "Topology check written to:"
echo "  $OUTDIR/$OUT_NAME.png"
echo "  $OUTDIR/$OUT_NAME.pdf"
echo "  $OUTDIR/${OUT_NAME}_layout.csv"
echo "  $OUTDIR/${OUT_NAME}_report.json"
echo
echo "Existing files with these names are overwritten on rerun."
