"""Read-only classification vocabulary for the DEV228 source-state gate."""

from __future__ import annotations


def source_inventory() -> list[dict]:
    """Existing preparations, deliberately not promoted to magnetic bodies."""
    return [
        {"id": "S1", "name": "externally maintained static constraint", "basis": "DEV159 stationary source/contact deformation", "persistence": "SOURCE_MAINTAINED", "localization": "EXACT_GEOMETRIC", "identity": "PREPARATION_PROVENANCE_ONLY", "orientation": "ABSENT", "magnet_eligible": False},
        {"id": "S2", "name": "released propagating excitation", "basis": "DEV182/DEV195/DEV203 packet", "persistence": "PROPAGATING_ONLY", "localization": "NONCOMPACT", "identity": "DYNAMIC_COHERENT", "orientation": "DYNAMIC_ONLY", "magnet_eligible": False},
        {"id": "S3", "name": "dynamically prepared full-state structure", "basis": "DEV196/DEV213 injection semantics", "persistence": "PROPAGATING_ONLY", "localization": "NONCOMPACT", "identity": "PREPARATION_PROVENANCE_ONLY", "orientation": "DYNAMIC_ONLY", "magnet_eligible": False},
        {"id": "S4", "name": "source-maintained deformation", "basis": "DEV211 six-neighbour source-contact equilibrium", "persistence": "SOURCE_MAINTAINED", "localization": "EXACT_GEOMETRIC", "identity": "PREPARATION_PROVENANCE_ONLY", "orientation": "ABSENT", "magnet_eligible": False},
        {"id": "S5", "name": "independent self-supported localized bounded state", "basis": "repository-wide audit", "persistence": "NOT_DERIVED", "localization": "ABSENT", "identity": "ABSENT", "orientation": "ABSENT", "magnet_eligible": False},
    ]


def classifications() -> dict:
    return {
        "NATIVE_SOURCE_PERSISTENCE": "NOT_DERIVED",
        "NATIVE_SOURCE_LOCALIZATION": "NONUNIQUE",
        "NATIVE_SOURCE_IDENTITY": "NONUNIQUE",
        "SOURCE_ORIENTATION_STATE": "DYNAMIC_ONLY",
        "MAGNET_LIKE_SOURCE_STATE_VALIDITY": "NOT_DERIVED",
        "TWO_BODY_SOURCE_COMPOSITION": "BLOCKED_PREPARATION",
        "FINITE_X_NATIVE_COMPOSITION_GATE": "DERIVED",
    }
