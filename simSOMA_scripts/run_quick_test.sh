#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${PROJECT_ROOT}"

unset LD_PRELOAD 2>/dev/null || true

if [[ ! -x "${PROJECT_ROOT}/.venv/bin/python" ]]; then
  echo ".venv not found. Building it with simSOMA_scripts/01_setup_env.sh."
  echo "For strict reproducibility, activate conda environment simsoma311 first:"
  echo '  source "$(conda info --base)/etc/profile.d/conda.sh"'
  echo "  conda activate simsoma311"
  echo
  bash simSOMA_scripts/01_setup_env.sh
fi

PY="${PROJECT_ROOT}/.venv/bin/python"
PY_VER="$(${PY} -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')"
if [[ "${PY_VER}" != "3.11" && "${SIMSOMA_ALLOW_OTHER_PYTHON:-0}" != "1" ]]; then
  echo "Refusing to run quick test with Python ${PY_VER}; expected 3.11." >&2
  echo "Rebuild .venv with:" >&2
  echo '  source "$(conda info --base)/etc/profile.d/conda.sh"' >&2
  echo "  conda activate simsoma311" >&2
  echo '  PYTHON_BIN="$(which python)" bash simSOMA_scripts/01_setup_env.sh' >&2
  exit 1
fi

export SIMSOMA_PYTHON="${PY}"

# Remove previous quick-test output so the run is repeatable.
rm -rf simSOMA_output/quick_test_2organs
rm -rf simSOMA_output/quick_test_2organs__psplit_*

echo
echo "Plotting quick-test topology..."
MASTER_CONFIG="simSOMA_configs/quick_test_2organs.json"
TOPOLOGY_JSON="simSOMA_inputs/examples/quick_test_2organs_topology.json"
SIMSOMA_PYTHON="${PY}" bash simSOMA_scripts/check_topology.sh "${TOPOLOGY_JSON}"
echo

if [[ -t 0 && "${SIMSOMA_SKIP_TOPOLOGY_CONFIRM:-0}" != "1" ]]; then
  echo "Open the topology plot/report now. Continue only if the topology is correct."
  read -r -p "Run the quick simulation now? [y/N] " REPLY
  case "${REPLY}" in
    y|Y|yes|YES) ;;
    *)
      echo "Stopped before simulation."
      exit 0
      ;;
  esac
fi

echo
echo "Running quick simulation through 00_pipeline.sh..."
SIMSOMA_PYTHON="${PY}" bash simSOMA_scripts/00_pipeline.sh run "${MASTER_CONFIG}" 1 1 yes parameter

echo
echo "Quick test finished. Output folder:"
echo "  simSOMA_output/quick_test_2organs/"
