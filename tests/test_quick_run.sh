#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${PROJECT_ROOT}"

SIMSOMA_SKIP_TOPOLOGY_CONFIRM=1 bash simSOMA_scripts/run_quick_test.sh

test -f simSOMA_output/quick_test_2organs/grid_parameter/aggregated_summaries.csv
test -f simSOMA_output/quick_test_2organs/grid_parameter/parameter_sets.csv
test -f simSOMA_output/topology_check/quick_test_2organs_topology/topology_plot.png
test -f simSOMA_output/topology_check/quick_test_2organs_topology/topology_plot_report.json

echo "simSOMA quick test passed."
