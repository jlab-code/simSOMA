#!/usr/bin/env python3
"""Backward-compatible wrapper; the converter is ``simsoma.topology_csv``
(also available as `simsoma topology-from-csv` or `python -m simsoma.topology_csv`)."""
import sys
from pathlib import Path

if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # repository checkout without installation
    from simsoma.topology_csv import main
    raise SystemExit(main())
