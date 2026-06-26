#!/usr/bin/env bash
set -euo pipefail

# Convert branch/organ CSV topology tables to the standard simSOMA
# user-friendly topology JSON format.
#
# Usage:
#   bash simSOMA_scripts/run_convert_topology_csv_to_json.sh \
#     simSOMA_inputs/examples/example_topology_branches.csv \
#     simSOMA_inputs/examples/example_topology_organs.csv \
#     simSOMA_inputs/examples/example_topology_from_csv.json \
#     years
#
# Positional arguments:
#   1  Branch CSV path
#   2  Organ CSV path
#   3  Output topology JSON path
#   4  Unit written into JSON: years, meters, or steps. Default: years
#
# Optional environment variable:
#   REPORT_JSON=/path/to/report.json

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"

if [[ $# -lt 3 || $# -gt 4 ]]; then
  echo "Usage:" >&2
  echo "  bash simSOMA_scripts/run_convert_topology_csv_to_json.sh <branches.csv> <organs.csv> <out.json> [years|meters|steps]" >&2
  exit 2
fi

BRANCHES_CSV="$1"
ORGANS_CSV="$2"
OUT_JSON="$3"
UNIT="${4:-years}"

case "$UNIT" in
  years|meters|steps) ;;
  *)
    echo "ERROR: unit must be one of: years, meters, steps" >&2
    exit 2
    ;;
esac

PY="${SIMSOMA_PYTHON:-}"
if [[ -z "$PY" ]]; then
  if [[ -x "$PROJECT_ROOT/.venv/bin/python" ]]; then
    PY="$PROJECT_ROOT/.venv/bin/python"
  else
    PY="$(command -v python3)"
  fi
fi

CMD=("$PY" simSOMA_scripts/convert_topology_csv_to_json.py)
CMD+=(--branches "$BRANCHES_CSV")
CMD+=(--organs "$ORGANS_CSV")
CMD+=(--out "$OUT_JSON")
CMD+=(--unit "$UNIT")

if [[ -n "${REPORT_JSON:-}" ]]; then
  CMD+=(--report "$REPORT_JSON")
fi

"${CMD[@]}"
