#!/usr/bin/env python3
"""Backward-compatible wrapper; the converter lives in simSOMA_corefunc/topology_csv.py
(also available as `simsoma topology-from-csv`)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "simSOMA_corefunc"))
from topology_csv import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
