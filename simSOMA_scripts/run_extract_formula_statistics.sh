#!/usr/bin/env bash
set -euo pipefail

# Wrapper for extracting fixed/intermediate/private/sharedness formula statistics.
#
# Usage:
#   bash run_extract_formula_statistics.sh <EXPERIMENT_OR_GRID_OR_OBSERVED_SCENARIO_DIR>
#
# Most common use:
#   bash run_extract_formula_statistics.sh simSOMA_output/quick_test_2organs
#
# This automatically detects:
#   - grid_parameter/vaf_count_spectrum_aggregated_summaries.csv
#   - observation_model_transforms/*/observed_vaf_count_spectrum_aggregated_summaries.tsv
#
# Optional environment variables:
#   FORMULA_TARGET=all|primary|observed  default: all
#   LOW_THRESHOLD=0.05                default: 0.05
#   FIXED_FILTERS="--fixed rho=0"      default: empty
#   OVERWRITE=yes|no                  default: yes
#   PYTHON_BIN=/path/to/python        default: .venv/bin/python if present, else python3

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

if [[ -x "${PROJECT_ROOT}/.venv/bin/python" ]]; then
  PYTHON_BIN="${PYTHON_BIN:-${PROJECT_ROOT}/.venv/bin/python}"
else
  PYTHON_BIN="${PYTHON_BIN:-python3}"
fi

LOW_THRESHOLD="${LOW_THRESHOLD:-0.05}"
TARGET="${FORMULA_TARGET:-all}"
FIXED_FILTERS="${FIXED_FILTERS:-}"
OVERWRITE="${OVERWRITE:-yes}"

if [[ $# -lt 1 ]]; then
  echo "Usage: bash $0 <EXPERIMENT_OR_GRID_OR_OBSERVED_SCENARIO_DIR>" >&2
  exit 2
fi

RESULT_DIR="$1"

CMD=(
  "$PYTHON_BIN"
  "$PROJECT_ROOT/simSOMA_corefunc/extract_formula_statistics.py"
  --result-dir "$RESULT_DIR"
  --target "$TARGET"
  --low-threshold "$LOW_THRESHOLD"
)

if [[ "$OVERWRITE" == "no" || "$OVERWRITE" == "NO" || "$OVERWRITE" == "false" || "$OVERWRITE" == "FALSE" ]]; then
  CMD+=(--no-overwrite)
fi

# FIXED_FILTERS is intentionally split so users can pass e.g.:
#   FIXED_FILTERS="--fixed rho=0 --fixed m=min,max"
# shellcheck disable=SC2206
EXTRA_ARGS=( $FIXED_FILTERS )
CMD+=("${EXTRA_ARGS[@]}")

printf 'Running:'
printf ' %q' "${CMD[@]}"
printf '\n'
"${CMD[@]}"
