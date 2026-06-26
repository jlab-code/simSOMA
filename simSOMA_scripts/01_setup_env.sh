#!/usr/bin/env bash
set -euo pipefail

# Build the project-local simSOMA virtual environment.
# This is the single setup entry point for both laptop and cluster use.
#
# Default behavior:
#   1. Use PYTHON_BIN if explicitly provided.
#   2. Otherwise use python3.11 if available.
#   3. Otherwise use active python/python3 if it is Python 3.11.
#   4. Otherwise, if conda is available, create/activate CONDA_ENV=simsoma311
#      and build .venv from that Python 3.11 environment.
#
# The resulting runtime is always:
#   .venv/bin/python

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VENV_DIR="${PROJECT_ROOT}/.venv"

PREFERRED_PYTHON="${PREFERRED_PYTHON:-python3.11}"
PYTHON_BIN="${PYTHON_BIN:-}"
CONDA_ENV="${CONDA_ENV:-simsoma311}"
CONDA_SH="${CONDA_SH:-}"
ALLOW_OTHER_PYTHON="${SIMSOMA_ALLOW_OTHER_PYTHON:-0}"

unset LD_PRELOAD 2>/dev/null || true
unset PYTHONPATH 2>/dev/null || true
unset PYTHONHOME 2>/dev/null || true
export PYTHONNOUSERSITE=1
hash -r

python_minor_version() {
  "$1" -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")' 2>/dev/null
}

is_python_311() {
  local candidate="$1"
  local version
  version="$(python_minor_version "$candidate" || true)"
  [[ "$version" == "3.11" ]]
}

find_conda_sh() {
  if [[ -n "${CONDA_SH}" && -f "${CONDA_SH}" ]]; then
    echo "${CONDA_SH}"
    return 0
  fi

  if command -v conda >/dev/null 2>&1; then
    local base
    base="$(conda info --base 2>/dev/null || true)"
    if [[ -n "$base" && -f "$base/etc/profile.d/conda.sh" ]]; then
      echo "$base/etc/profile.d/conda.sh"
      return 0
    fi
  fi

  for candidate in \
    "$HOME/miniconda3/etc/profile.d/conda.sh" \
    "$HOME/anaconda3/etc/profile.d/conda.sh" \
    "/opt/miniconda3/etc/profile.d/conda.sh" \
    "/opt/anaconda3/etc/profile.d/conda.sh"; do
    if [[ -f "$candidate" ]]; then
      echo "$candidate"
      return 0
    fi
  done

  return 1
}

choose_python() {
  # 1. Explicit interpreter.
  if [[ -n "${PYTHON_BIN}" ]]; then
    command -v "${PYTHON_BIN}" >/dev/null 2>&1 || {
      echo "Requested Python executable not found: ${PYTHON_BIN}" >&2
      exit 1
    }
    echo "${PYTHON_BIN}"
    return
  fi

  # 2. Preferred executable name.
  if command -v "${PREFERRED_PYTHON}" >/dev/null 2>&1; then
    echo "${PREFERRED_PYTHON}"
    return
  fi

  # 3. Active shell Python, if already 3.11.
  for candidate in python python3; do
    if command -v "${candidate}" >/dev/null 2>&1 && is_python_311 "${candidate}"; then
      echo "${candidate}"
      return
    fi
  done

  # 4. Conda fallback: make Python 3.11 available without requiring conda init.
  local conda_setup
  if conda_setup="$(find_conda_sh)"; then
    # shellcheck disable=SC1090
    source "$conda_setup"
    if ! conda env list | awk '{print $1}' | grep -qx "$CONDA_ENV"; then
      echo "Creating conda environment: $CONDA_ENV" >&2
      conda create -n "$CONDA_ENV" python=3.11 -y >&2
    fi
    conda activate "$CONDA_ENV"
    if is_python_311 python; then
      echo "$(command -v python)"
      return
    fi
  fi

  # 5. Deliberate non-3.11 fallback.
  if [[ "${ALLOW_OTHER_PYTHON}" == "1" ]]; then
    for candidate in python3 python; do
      if command -v "${candidate}" >/dev/null 2>&1; then
        echo "${candidate}"
        return
      fi
    done
  fi

  echo "Could not find Python 3.11." >&2
  echo "The simplest fix is:" >&2
  echo "  bash simSOMA_scripts/01_setup_env.sh" >&2
  echo "" >&2
  echo "If conda is installed but not found automatically, provide CONDA_SH:" >&2
  echo "  CONDA_SH=\$HOME/miniconda3/etc/profile.d/conda.sh bash simSOMA_scripts/01_setup_env.sh" >&2
  echo "" >&2
  echo "Or provide an explicit interpreter:" >&2
  echo "  PYTHON_BIN=/path/to/python3.11 bash simSOMA_scripts/01_setup_env.sh" >&2
  echo "" >&2
  echo "Not recommended, but possible:" >&2
  echo "  SIMSOMA_ALLOW_OTHER_PYTHON=1 bash simSOMA_scripts/01_setup_env.sh" >&2
  exit 1
}

PYTHON_BIN="$(choose_python)"
PY_VER="$("${PYTHON_BIN}" -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')"

if [[ "${PY_VER}" != "3.11" && "${ALLOW_OTHER_PYTHON}" != "1" ]]; then
  echo "Refusing to build .venv with Python ${PY_VER}; expected 3.11 for reproducibility." >&2
  echo "Use:" >&2
  echo "  bash simSOMA_scripts/01_setup_env.sh" >&2
  echo "or:" >&2
  echo "  PYTHON_BIN=/path/to/python3.11 bash simSOMA_scripts/01_setup_env.sh" >&2
  exit 1
fi

REQ_FILE="${PROJECT_ROOT}/requirements.lock.txt"
if [[ ! -f "${REQ_FILE}" ]]; then
  REQ_FILE="${PROJECT_ROOT}/requirements.txt"
fi

echo "Project root:      ${PROJECT_ROOT}"
echo "Virtual env dir:   ${VENV_DIR}"
echo "Python selected:   $(command -v "${PYTHON_BIN}" || echo "${PYTHON_BIN}")"
echo "Python version:    $("${PYTHON_BIN}" --version 2>&1)"
echo "Requirements file: ${REQ_FILE}"

REBUILD_VENV=0
if [[ ! -x "${VENV_DIR}/bin/python" ]]; then
  REBUILD_VENV=1
else
  EXISTING_PREFIX="$("${VENV_DIR}/bin/python" -c 'import sys; print(sys.base_prefix)')"
  NEW_PREFIX="$("${PYTHON_BIN}" -c 'import sys; print(sys.base_prefix)')"
  EXISTING_VER="$("${VENV_DIR}/bin/python" -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')"
  if [[ "${EXISTING_PREFIX}" != "${NEW_PREFIX}" || "${EXISTING_VER}" != "${PY_VER}" ]]; then
    REBUILD_VENV=1
  fi
fi

if [[ "${REBUILD_VENV}" == "1" ]]; then
  rm -rf "${VENV_DIR}"
  "${PYTHON_BIN}" -m venv "${VENV_DIR}"
fi

"${VENV_DIR}/bin/python" -m pip install --upgrade pip
"${VENV_DIR}/bin/python" -m pip install -r "${REQ_FILE}"

echo
echo "Virtual environment ready at:"
echo "  ${VENV_DIR}"
echo
echo "Runtime Python:"
echo "  ${VENV_DIR}/bin/python"
