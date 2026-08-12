"""Cycle status remains downstream of independently derived reversals."""
from __future__ import annotations


def cycle_status(turnaround: str, bounce: str) -> str:
    return "DERIVED" if turnaround == bounce == "DERIVED" else "NOT_DERIVED"
