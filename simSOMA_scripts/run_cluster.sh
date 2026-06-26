#!/usr/bin/env bash
set -euo pipefail

# Generic simSOMA cluster run launcher.
# Usage, from project root:
#   MASTER_CONFIG="simSOMA_configs/my_config.json" bash simSOMA_scripts/run_cluster.sh
# or:
#   bash simSOMA_scripts/run_cluster.sh simSOMA_configs/my_config.json
# Optional overrides:
#   N_SPLITS=300 N_PARALLEL_JOBS=20 CLEANUP_SPLITS=yes SPLIT_AXIS=parameter bash simSOMA_scripts/run_cluster.sh <config>

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PROJECT_DIR"

if [[ $# -ge 1 ]]; then
  MASTER_CONFIG="$1"
else
  MASTER_CONFIG="${MASTER_CONFIG:-}"
fi

if [[ -z "$MASTER_CONFIG" ]]; then
  echo "ERROR: MASTER_CONFIG is not set."
  echo "Use one of:"
  echo "  MASTER_CONFIG=\"simSOMA_configs/my_config.json\" bash simSOMA_scripts/run_cluster.sh"
  echo "  bash simSOMA_scripts/run_cluster.sh simSOMA_configs/my_config.json"
  exit 2
fi

N_SPLITS="${N_SPLITS:-300}"
N_PARALLEL_JOBS="${N_PARALLEL_JOBS:-20}"
CLEANUP_SPLITS="${CLEANUP_SPLITS:-yes}"
SPLIT_AXIS="${SPLIT_AXIS:-parameter}"
SIMSOMA_PYTHON="${SIMSOMA_PYTHON:-$PROJECT_DIR/.venv/bin/python}"

export CLEANUP_SPLITS
export SPLIT_AXIS
export SIMSOMA_PYTHON

echo "== simSOMA generic cluster run =="
echo "Project:          $PROJECT_DIR"
echo "MASTER_CONFIG:    $MASTER_CONFIG"
echo "N_SPLITS:         $N_SPLITS"
echo "N_PARALLEL_JOBS:  $N_PARALLEL_JOBS"
echo "CLEANUP_SPLITS:   $CLEANUP_SPLITS"
echo "SPLIT_AXIS:       $SPLIT_AXIS"
echo "SIMSOMA_PYTHON:   $SIMSOMA_PYTHON"
echo

if [[ ! -f "$MASTER_CONFIG" ]]; then
  echo "ERROR: Config file does not exist: $MASTER_CONFIG"
  exit 3
fi

if [[ ! -x "$SIMSOMA_PYTHON" ]]; then
  echo "ERROR: SIMSOMA_PYTHON is not executable: $SIMSOMA_PYTHON"
  echo "Run simSOMA_scripts/01_setup_env.sh first, or use cluster_prepare_and_start.sh."
  exit 4
fi

bash simSOMA_scripts/00_pipeline.sh run "$MASTER_CONFIG" "$N_SPLITS" "$N_PARALLEL_JOBS" "$CLEANUP_SPLITS" "$SPLIT_AXIS"
