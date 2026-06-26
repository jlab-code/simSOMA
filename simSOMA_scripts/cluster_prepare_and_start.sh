#!/usr/bin/env bash
set -euo pipefail

# Cluster preflight + detached tmux starter.
# Same environment workflow as local use:
#   bash simSOMA_scripts/01_setup_env.sh
# Then this wrapper starts the selected config inside tmux.
# If .venv is missing, this script calls 01_setup_env.sh once.

PROJECT_DIR="${PROJECT_DIR:-$PWD}"
SESSION="${SESSION:-simsoma}"
LOG_FILE="${LOG_FILE:-$PROJECT_DIR/cluster_run.log}"

if [[ $# -ge 1 ]]; then
  MASTER_CONFIG="$1"
else
  MASTER_CONFIG="${MASTER_CONFIG:-}"
fi

if [[ -z "$MASTER_CONFIG" ]]; then
  echo "ERROR: MASTER_CONFIG is not set." >&2
  echo "Use:" >&2
  echo "  bash simSOMA_scripts/cluster_prepare_and_start.sh simSOMA_configs/my_config.json" >&2
  exit 2
fi

cd "$PROJECT_DIR"

N_SPLITS="${N_SPLITS:-300}"
N_PARALLEL_JOBS="${N_PARALLEL_JOBS:-20}"
CLEANUP_SPLITS="${CLEANUP_SPLITS:-yes}"
SPLIT_AXIS="${SPLIT_AXIS:-parameter}"
RUN_CLUSTER_CHECK="${RUN_CLUSTER_CHECK:-no}"
SIMSOMA_PYTHON="${SIMSOMA_PYTHON:-$PWD/.venv/bin/python}"

chmod +x simSOMA_scripts/00_pipeline.sh simSOMA_scripts/01_setup_env.sh simSOMA_scripts/run_cluster.sh 2>/dev/null || true

if [[ ! -x "$SIMSOMA_PYTHON" ]]; then
  echo "Local .venv not found. Building it with simSOMA_scripts/01_setup_env.sh."
  bash simSOMA_scripts/01_setup_env.sh
  SIMSOMA_PYTHON="$PWD/.venv/bin/python"
fi

if [[ ! -x "$SIMSOMA_PYTHON" ]]; then
  echo "ERROR: SIMSOMA_PYTHON is not executable: $SIMSOMA_PYTHON" >&2
  exit 3
fi

export MASTER_CONFIG N_SPLITS N_PARALLEL_JOBS CLEANUP_SPLITS SPLIT_AXIS RUN_CLUSTER_CHECK SIMSOMA_PYTHON

echo "== simSOMA cluster preflight =="
echo "Project:          $PROJECT_DIR"
echo "Config:           $MASTER_CONFIG"
echo "tmux session:     $SESSION"
echo "Log:              $LOG_FILE"
echo "N_SPLITS:         $N_SPLITS"
echo "N_PARALLEL_JOBS:  $N_PARALLEL_JOBS"
echo "CLEANUP_SPLITS:   $CLEANUP_SPLITS"
echo "SPLIT_AXIS:       $SPLIT_AXIS"
echo "RUN_CLUSTER_CHECK:$RUN_CLUSTER_CHECK"
echo "SIMSOMA_PYTHON:   $SIMSOMA_PYTHON"
echo

echo "== Python checks =="
"$SIMSOMA_PYTHON" --version
"$SIMSOMA_PYTHON" - <<'PY'
import numpy, matplotlib
print('numpy/matplotlib import: ok')
PY

echo
echo "== Version checks =="
if [[ -f simSOMA_scripts/02_version_report.sh ]]; then
  bash simSOMA_scripts/02_version_report.sh || true
fi
[[ -f VERSION ]] && cat VERSION || true

ACTUAL_CHECKSUM="$(sha256sum simSOMA_corefunc/run_from_config.py | awk '{print $1}')"
echo "run_from_config.py checksum: $ACTUAL_CHECKSUM"

if [[ ! -f "$MASTER_CONFIG" ]]; then
  echo "ERROR: MASTER_CONFIG does not exist: $MASTER_CONFIG" >&2
  exit 5
fi

echo
if [[ "$RUN_CLUSTER_CHECK" == "yes" ]]; then
  echo "== Optional cluster-side topology/config check =="
  SIMSOMA_PYTHON="$SIMSOMA_PYTHON" bash simSOMA_scripts/00_pipeline.sh inspect "$MASTER_CONFIG"
else
  echo "== Skipping cluster-side topology/config check =="
  echo "RUN_CLUSTER_CHECK is set to '$RUN_CLUSTER_CHECK'."
fi
echo

if tmux has-session -t "$SESSION" 2>/dev/null; then
  echo "ERROR: tmux session '$SESSION' already exists." >&2
  echo "Use: tmux attach -t $SESSION" >&2
  exit 6
fi

cat > simSOMA_scripts/_run_cluster_inside_tmux.sh <<EOF2
#!/usr/bin/env bash
set -euo pipefail
cd "$PROJECT_DIR"
unset LD_PRELOAD 2>/dev/null || true
unset PYTHONPATH 2>/dev/null || true
unset PYTHONHOME 2>/dev/null || true
export PYTHONNOUSERSITE=1
export MASTER_CONFIG="$MASTER_CONFIG"
export N_SPLITS="$N_SPLITS"
export N_PARALLEL_JOBS="$N_PARALLEL_JOBS"
export CLEANUP_SPLITS="$CLEANUP_SPLITS"
export SPLIT_AXIS="$SPLIT_AXIS"
export RUN_CLUSTER_CHECK="$RUN_CLUSTER_CHECK"
export SIMSOMA_PYTHON="\$PWD/.venv/bin/python"
bash simSOMA_scripts/run_cluster.sh "\$MASTER_CONFIG" 2>&1 | tee "$LOG_FILE"
EOF2
chmod +x simSOMA_scripts/_run_cluster_inside_tmux.sh

echo "== Starting long run inside detached tmux session =="
tmux new-session -d -s "$SESSION" "bash simSOMA_scripts/_run_cluster_inside_tmux.sh"

echo
echo "Run started. You may log out."
echo "Follow log:       tail -f $LOG_FILE"
echo "Attach tmux:      tmux attach -t $SESSION"
echo "List processes:   pgrep -a -u \"\$USER\" -f '00_pipeline.sh|launch_grid_splits.py|run_from_config.py'"
echo "Stop run:         pkill -u \"\$USER\" -f '00_pipeline.sh|launch_grid_splits.py|run_from_config.py'"
