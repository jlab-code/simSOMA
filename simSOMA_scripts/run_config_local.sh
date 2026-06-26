#!/usr/bin/env bash
set -euo pipefail

# Generic local simSOMA launcher.
# Usage:
#   bash simSOMA_scripts/run_config_local.sh simSOMA_configs/simulation_01_topology_depth_rough.json
# Optional:
#   N_SPLITS=4 N_PARALLEL_JOBS=2 RUN_CHECK=yes bash simSOMA_scripts/run_config_local.sh <config>

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"

CONFIG="${1:-${MASTER_CONFIG:-}}"
if [[ -z "$CONFIG" ]]; then
  echo "ERROR: provide a config path or set MASTER_CONFIG." >&2
  exit 2
fi

N_SPLITS="${N_SPLITS:-1}"
N_PARALLEL_JOBS="${N_PARALLEL_JOBS:-1}"
CLEANUP_SPLITS="${CLEANUP_SPLITS:-yes}"
SPLIT_AXIS="${SPLIT_AXIS:-parameter}"
RUN_CHECK="${RUN_CHECK:-yes}"

unset LD_PRELOAD 2>/dev/null || true
export SIMSOMA_PYTHON="${SIMSOMA_PYTHON:-$PROJECT_ROOT/.venv/bin/python}"

if [[ ! -x "$SIMSOMA_PYTHON" ]]; then
  echo "Local .venv not found or SIMSOMA_PYTHON is not executable. Building .venv." >&2
  echo "The setup script will use PYTHON_BIN if set; otherwise it will search for Python 3.11." >&2
  if [[ -n "${PYTHON_BIN:-}" ]]; then
    PYTHON_BIN="$PYTHON_BIN" bash simSOMA_scripts/01_setup_env.sh
  else
    bash simSOMA_scripts/01_setup_env.sh
  fi
  export SIMSOMA_PYTHON="$PROJECT_ROOT/.venv/bin/python"
fi

if [[ "$RUN_CHECK" == "yes" ]]; then
  bash simSOMA_scripts/00_pipeline.sh inspect "$CONFIG" 1 1 yes parameter
fi

bash simSOMA_scripts/00_pipeline.sh run "$CONFIG" "$N_SPLITS" "$N_PARALLEL_JOBS" "$CLEANUP_SPLITS" "$SPLIT_AXIS"
