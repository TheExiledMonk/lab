"""Read-only classification helpers for the DEV229 source gate.

These helpers intentionally do not evolve a state or construct a source.  They
make the gate's consequence explicit: a self-supported trajectory needs a
canonical preparation/release before persistence can be evaluated.
"""
from __future__ import annotations


def derivation_from_release(release_semantics: str) -> str:
    """Return the only allowed top-level result for a release-semantics block."""
    if release_semantics == "NOT_DERIVED":
        return "BLOCKED_SOURCE_RELEASE"
    raise ValueError("DEV229 only classifies the recovered NOT_DERIVED release state")


def localization_from_noncompact(noncompact: bool) -> str:
    """Classify the existing packet support without choosing a threshold."""
    return "NONUNIQUE" if noncompact else "UNASSESSED"
