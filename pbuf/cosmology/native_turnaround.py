"""Classification helpers; no new background evolution equation is supplied."""
from __future__ import annotations


def turnaround_status(scale_mapping: str, evolution_equation: bool) -> str:
    if scale_mapping not in {"ALREADY_DERIVED", "DERIVABLE_WITHOUT_NEW_ASSUMPTION"}:
        return "BLOCKED"
    return "UNRESOLVED" if not evolution_equation else "UNRESOLVED"
