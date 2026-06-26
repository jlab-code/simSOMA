"""event_schema.py

Standardized event record schema used by the pipeline wrapper.

This schema is intentionally simple (plain dict-compatible dataclasses).
It exists to enforce stable module I/O, enabling plug-in replacement.

No file I/O here.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Any, Dict, Tuple


Time2D = Tuple[int, int]  # (t_age, t_dev)


@dataclass(frozen=True)
class EventRecord:
    module: str
    event_type: str
    branch_id: str
    target_id: str
    event_time: Time2D
    outputs: Dict[str, Any]

    def to_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        d["event_time"] = [int(self.event_time[0]), int(self.event_time[1])]
        return d


def require_keys(d: Dict[str, Any], keys: Tuple[str, ...], where: str) -> None:
    for k in keys:
        if k not in d:
            raise KeyError(f"Missing key '{k}' in {where}")


def validate_organ_outputs(outputs: Dict[str, Any]) -> None:
    require_keys(outputs, ("sequenced_cells", "allele_counts_by_mutation"), where="organ.outputs")
    ac = outputs["allele_counts_by_mutation"]
    if not isinstance(ac, dict):
        raise TypeError("allele_counts_by_mutation must be a dict")
    sequenced_cells = int(outputs["sequenced_cells"])
    for k, v in ac.items():
        vv = int(v)
        if vv < 1 or vv > sequenced_cells:
            raise ValueError(f"allele_counts_by_mutation[{k}]={vv} outside 1..sequenced_cells={sequenced_cells}")
