#!/usr/bin/env bash
set -euo pipefail

# Wrapper for applying one observation-model transform to an existing simSOMA run.
# No observation-model JSON file is used.
#
# Required environment variables or positional arguments:
#   RUN_DIR          or $1  simSOMA experiment dir, or grid_parameter dir
#   SCENARIO_NAME   or $2  output scenario name
#   LAYER_NAMES     or $3  comma-separated layer names, e.g. L1,L2,L3
#   LAYER_WEIGHTS   or $4  comma-separated layer weights, e.g. 0.10,0.70,0.20
#   SAMPLING_MODE   or $5  bulk or layer_specific
#   PHASE_MODE      or $6  phased or unphased
#
# Optional environment variables:
#   TARGET_LAYERS=all
#   LAYER_WEIGHT_MODE=deterministic
#   LAYER_MUTATION_RATE_MULTIPLIERS=1,1,1
#   NORMALIZE_LAYER_WEIGHTS=no
#   OBSERVED_NBINS=50
#   OVERWRITE=no
#   OUT_BASE_DIR=/optional/output/base
#   SKIP_REPLICATE=no
#   SKIP_AGGREGATED=no
#   SOURCE_VAF_TABLE=/optional/source/table
#   SOURCE_VAF_REPLICATE_TABLE=/optional/source/table
#   PYTHON_BIN=/path/to/python

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-${PROJECT_ROOT}/.venv/bin/python}"

RUN_DIR="${RUN_DIR:-${1:-}}"
SCENARIO_NAME="${SCENARIO_NAME:-${2:-}}"
LAYER_NAMES="${LAYER_NAMES:-${3:-}}"
LAYER_WEIGHTS="${LAYER_WEIGHTS:-${4:-}}"
SAMPLING_MODE="${SAMPLING_MODE:-${5:-}}"
PHASE_MODE="${PHASE_MODE:-${6:-}}"

TARGET_LAYERS="${TARGET_LAYERS:-all}"
LAYER_WEIGHT_MODE="${LAYER_WEIGHT_MODE:-deterministic}"
LAYER_MUTATION_RATE_MULTIPLIERS="${LAYER_MUTATION_RATE_MULTIPLIERS:-}"
NORMALIZE_LAYER_WEIGHTS="${NORMALIZE_LAYER_WEIGHTS:-no}"
OBSERVED_NBINS="${OBSERVED_NBINS:-50}"
OVERWRITE="${OVERWRITE:-no}"
OUT_BASE_DIR="${OUT_BASE_DIR:-}"
SKIP_REPLICATE="${SKIP_REPLICATE:-no}"
SKIP_AGGREGATED="${SKIP_AGGREGATED:-no}"
SOURCE_VAF_TABLE="${SOURCE_VAF_TABLE:-}"
SOURCE_VAF_REPLICATE_TABLE="${SOURCE_VAF_REPLICATE_TABLE:-}"

if [[ -z "$RUN_DIR" || -z "$SCENARIO_NAME" || -z "$LAYER_NAMES" || -z "$LAYER_WEIGHTS" || -z "$SAMPLING_MODE" || -z "$PHASE_MODE" ]]; then
  cat >&2 <<'EOF'
Usage:
  bash run_observation_model_transform.sh <RUN_OR_GRID_DIR> <SCENARIO_NAME> <LAYER_NAMES> <LAYER_WEIGHTS> <SAMPLING_MODE> <PHASE_MODE>

Example:
  bash simSOMA_scripts/run_observation_model_transform.sh \
    simSOMA_output/quick_test_2organs \
    three_layer_bulk_unphased_L1_0p1_L2_0p7_L3_0p2 \
    L1,L2,L3 \
    0.10,0.70,0.20 \
    bulk \
    unphased
EOF
  exit 2
fi

if [[ ! -x "$PYTHON_BIN" ]]; then
  if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="$(command -v python3)"
  else
    echo "ERROR: Could not find Python. Set PYTHON_BIN=/path/to/python." >&2
    exit 1
  fi
fi

CMD=(
  "$PYTHON_BIN"
  "$PROJECT_ROOT/simSOMA_corefunc/transform_observation_model.py"
  --run-dir "$RUN_DIR"
  --scenario-name "$SCENARIO_NAME"
  --layer-names "$LAYER_NAMES"
  --layer-weights "$LAYER_WEIGHTS"
  --sampling-mode "$SAMPLING_MODE"
  --phase-mode "$PHASE_MODE"
  --layer-weight-mode "$LAYER_WEIGHT_MODE"
  --target-layers "$TARGET_LAYERS"
  --observed-nbins "$OBSERVED_NBINS"
)

if [[ -n "$LAYER_MUTATION_RATE_MULTIPLIERS" ]]; then
  CMD+=(--layer-mutation-rate-multipliers "$LAYER_MUTATION_RATE_MULTIPLIERS")
fi
if [[ "$NORMALIZE_LAYER_WEIGHTS" == "yes" || "$NORMALIZE_LAYER_WEIGHTS" == "true" || "$NORMALIZE_LAYER_WEIGHTS" == "1" ]]; then
  CMD+=(--normalize-layer-weights)
fi
if [[ -n "$OUT_BASE_DIR" ]]; then
  CMD+=(--out-base-dir "$OUT_BASE_DIR")
fi
if [[ "$OVERWRITE" == "yes" || "$OVERWRITE" == "true" || "$OVERWRITE" == "1" ]]; then
  CMD+=(--overwrite)
fi
if [[ "$SKIP_REPLICATE" == "yes" || "$SKIP_REPLICATE" == "true" || "$SKIP_REPLICATE" == "1" ]]; then
  CMD+=(--skip-replicate)
fi
if [[ "$SKIP_AGGREGATED" == "yes" || "$SKIP_AGGREGATED" == "true" || "$SKIP_AGGREGATED" == "1" ]]; then
  CMD+=(--skip-aggregated)
fi
if [[ -n "$SOURCE_VAF_TABLE" ]]; then
  CMD+=(--source-vaf-table "$SOURCE_VAF_TABLE")
fi
if [[ -n "$SOURCE_VAF_REPLICATE_TABLE" ]]; then
  CMD+=(--source-vaf-replicate-table "$SOURCE_VAF_REPLICATE_TABLE")
fi

printf 'Running:'
printf ' %q' "${CMD[@]}"
printf '\n'
"${CMD[@]}"
