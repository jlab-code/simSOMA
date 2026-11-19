#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

sha256_rel() {
  (cd "${PROJECT_ROOT}" && sha256sum "$@")
}

echo "Project root: ${PROJECT_ROOT}"
echo "Version file: VERSION"
if [[ -f "${PROJECT_ROOT}/VERSION" ]]; then
  echo "Version: $(cat "${PROJECT_ROOT}/VERSION")"
else
  echo "Version: <missing VERSION file>"
fi

echo
echo "Tracked source checksums:"
for f in   simsoma/run_from_config.py   simsoma/pipeline_wrapper.py   simsoma/self_renewal.py   simsoma/pre_branching.py   simsoma/branching.py   simsoma/organ.py   simsoma/summaries.py   simsoma/branch_bias.py   simsoma/topology_io.py   simsoma/inspect_topology.py   simsoma/launch_grid_splits.py   simsoma/plot_topology_json.py
  do
    if [[ -f "${PROJECT_ROOT}/${f}" ]]; then
      sha256_rel "${f}"
    else
      echo "MISSING ${f}"
    fi
  done

echo
echo "Tracked script checksums:"
if compgen -G "${PROJECT_ROOT}/simSOMA_scripts/*.sh" > /dev/null; then
  (cd "${PROJECT_ROOT}" && find simSOMA_scripts -maxdepth 1 -type f -name '*.sh' -print0 | sort -z | xargs -0 sha256sum)
else
  echo "MISSING simSOMA_scripts/*.sh"
fi

echo
echo "Config checksums:"
if [[ -d "${PROJECT_ROOT}/simSOMA_configs" ]]; then
  (cd "${PROJECT_ROOT}" && find simSOMA_configs -type f -name '*.json' -print0 | sort -z | xargs -0 sha256sum)
else
  echo "MISSING simSOMA_configs"
fi

echo
echo "Input topology checksums:"
if [[ -d "${PROJECT_ROOT}/simSOMA_inputs" ]]; then
  (cd "${PROJECT_ROOT}" && find simSOMA_inputs -type f -name '*.json' -print0 | sort -z | xargs -0 sha256sum)
else
  echo "MISSING simSOMA_inputs"
fi
