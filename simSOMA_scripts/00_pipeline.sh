#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MODE="${1:-inspect}"
MASTER_CONFIG_REL="${2:-simSOMA_configs/master_config_template.json}"
N_SPLITS="${3:-8}"
N_PARALLEL_JOBS="${4:-8}"
CLEANUP_SPLITS="${5:-no}"
SPLIT_AXIS="${6:-parameter}"
ALLOW_OTHER_PYTHON="${SIMSOMA_ALLOW_OTHER_PYTHON:-0}"

case "${MODE}" in
  inspect|run|all) ;;
  *)
    echo "Invalid mode: ${MODE}" >&2
    echo "Use one of: inspect | run | all" >&2
    exit 1
    ;;
esac

MASTER_CONFIG="${PROJECT_ROOT}/${MASTER_CONFIG_REL}"
if [[ ! -f "${MASTER_CONFIG}" ]]; then
  echo "Master config not found: ${MASTER_CONFIG}" >&2
  exit 1
fi

# Clean runtime environment
unset LD_PRELOAD 2>/dev/null || true
unset PYTHONPATH 2>/dev/null || true
unset PYTHONHOME 2>/dev/null || true
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
hash -r

export SIMSOMA_HOME="${PROJECT_ROOT}"
export SIMSOMA_CODE_DIR="${SIMSOMA_HOME}/simSOMA_corefunc"
export SIMSOMA_INPUT_DIR="${SIMSOMA_HOME}/simSOMA_inputs"
export SIMSOMA_CONFIG_DIR="${SIMSOMA_HOME}/simSOMA_configs"
export SIMSOMA_OUTPUT_DIR="${SIMSOMA_HOME}/simSOMA_output"
mkdir -p "${SIMSOMA_OUTPUT_DIR}"

choose_python() {
  if [[ -n "${SIMSOMA_PYTHON:-}" ]]; then
    if [[ ! -x "${SIMSOMA_PYTHON}" ]]; then
      echo "SIMSOMA_PYTHON is not executable: ${SIMSOMA_PYTHON}" >&2
      exit 1
    fi
    echo "${SIMSOMA_PYTHON}"
    return
  fi

  if [[ -x "${PROJECT_ROOT}/.venv/bin/python" ]]; then
    echo "${PROJECT_ROOT}/.venv/bin/python"
    return
  fi

  if command -v python >/dev/null 2>&1; then
    command -v python
    return
  fi

  if command -v python3 >/dev/null 2>&1; then
    command -v python3
    return
  fi

  echo "No Python interpreter found." >&2
  echo "Either create the project venv or set SIMSOMA_PYTHON explicitly." >&2
  exit 1
}

PYTHON_CMD="$(choose_python)"
PY_VER="$("${PYTHON_CMD}" -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')"

if [[ "${PY_VER}" != "3.11" && "${ALLOW_OTHER_PYTHON}" != "1" ]]; then
  echo "Refusing to run with Python ${PY_VER}; expected 3.11 for reproducibility." >&2
  echo "Use one of:" >&2
  echo "  PYTHON_BIN=python3.11 bash simSOMA_scripts/01_setup_env.sh" >&2
  echo "  SIMSOMA_PYTHON=/path/to/python3.11 bash simSOMA_scripts/00_pipeline.sh ..." >&2
  echo "  SIMSOMA_ALLOW_OTHER_PYTHON=1 bash simSOMA_scripts/00_pipeline.sh ..." >&2
  exit 1
fi

REQ_FILE="${PROJECT_ROOT}/requirements.lock.txt"
if [[ ! -f "${REQ_FILE}" ]]; then
  REQ_FILE="${PROJECT_ROOT}/requirements.txt"
fi

check_environment() {
  "${PYTHON_CMD}" - <<'PY'
import importlib
import importlib.metadata as md
import os
import pathlib
import re
import sys

project_root = pathlib.Path(os.environ["SIMSOMA_HOME"])
req_path = project_root / ("requirements.lock.txt" if (project_root / "requirements.lock.txt").exists() else "requirements.txt")
if not req_path.exists():
    print(f"Requirements file not found: {req_path}", file=sys.stderr)
    raise SystemExit(1)

errors = []
checked = []

def parse_version_tuple(s: str):
    nums = re.findall(r"\d+", s)
    return tuple(int(x) for x in nums[:6])

# package-name -> import-name mapping
IMPORT_NAME_MAP = {
    "fonttools": "fontTools",
    "pillow": "PIL",
    "python-dateutil": "dateutil",
}

for raw in req_path.read_text().splitlines():
    line = raw.strip()
    if not line or line.startswith("#"):
        continue

    m_eq = re.match(r"^([A-Za-z0-9_.-]+)==([0-9][A-Za-z0-9_.-]*)$", line)
    m_ge = re.match(r"^([A-Za-z0-9_.-]+)>=([0-9][A-Za-z0-9_.-]*)$", line)

    if m_eq:
        package = m_eq.group(1)
        required = m_eq.group(2)
        mode = "eq"
    elif m_ge:
        package = m_ge.group(1)
        required = m_ge.group(2)
        mode = "ge"
    else:
        package = line.split()[0]
        required = None
        mode = None

    package_norm = package.lower()
    module_name = IMPORT_NAME_MAP.get(package_norm, package.replace("-", "_"))

    try:
        importlib.import_module(module_name)
    except Exception as exc:
        errors.append(f"missing import: {module_name} (from package {package}; {exc})")
        continue

    try:
        version = md.version(package)
    except Exception:
        version = None

    if required and version is not None:
        v = parse_version_tuple(version)
        r = parse_version_tuple(required)
        if mode == "eq" and v != r:
            errors.append(f"{package} {version} found, but == {required} required")
            continue
        if mode == "ge" and v < r:
            errors.append(f"{package} {version} found, but >= {required} required")
            continue

    checked.append((package, version))

if errors:
    print("Current Python environment is missing required dependencies:", file=sys.stderr)
    for e in errors:
        print(f"  - {e}", file=sys.stderr)
    print("", file=sys.stderr)
    print("Create or rebuild the project-local environment with:", file=sys.stderr)
    print("  bash simSOMA_scripts/01_setup_env.sh", file=sys.stderr)
    raise SystemExit(1)

print("Dependency check passed:")
for package, version in checked:
    if version is None:
        print(f"  - {package}")
    else:
        print(f"  - {package} {version}")
PY
}

echo "Project root:      ${PROJECT_ROOT}"
echo "Mode:              ${MODE}"
echo "Master config:     ${MASTER_CONFIG}"
echo "Grid splits:       ${N_SPLITS}"
echo "Parallel jobs:     ${N_PARALLEL_JOBS}"
echo "Cleanup splits:    ${CLEANUP_SPLITS}"
echo "Split axis:        ${SPLIT_AXIS}"
echo "Python:            ${PYTHON_CMD}"
"${PYTHON_CMD}" --version

echo "SIMSOMA_HOME=${SIMSOMA_HOME}"
echo "SIMSOMA_CODE_DIR=${SIMSOMA_CODE_DIR}"
echo "SIMSOMA_OUTPUT_DIR=${SIMSOMA_OUTPUT_DIR}"

echo
check_environment

cd "${SIMSOMA_CODE_DIR}"

run_check() {
  echo
  echo "[1/2] Checking topology"
  bash "${PROJECT_ROOT}/simSOMA_scripts/check_topology.sh" "${MASTER_CONFIG}"
  echo
  echo "Topology check completed. Review the generated topology outputs before launching the full run."
}

run_grid() {
  echo
  echo "[2/2] Launching split jobs"
  cmd=(
    "${PYTHON_CMD}" launch_grid_splits.py
    --master-config "${MASTER_CONFIG}"
    --n-splits "${N_SPLITS}"
    --jobs "${N_PARALLEL_JOBS}"
    --split-axis "${SPLIT_AXIS}"
  )

  if [[ "${CLEANUP_SPLITS}" == "yes" ]]; then
    cmd+=(--cleanup-splits)
  elif [[ "${CLEANUP_SPLITS}" != "no" ]]; then
    echo "Invalid CLEANUP_SPLITS value: ${CLEANUP_SPLITS} (use yes or no)" >&2
    exit 1
  fi

  "${cmd[@]}"
}

case "${MODE}" in
  inspect)
    run_check
    ;;
  run)
    run_grid
    ;;
  all)
    run_check
    run_grid
    ;;
esac
