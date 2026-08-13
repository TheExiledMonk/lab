"""Read-only DEV213 aggregate-state diagnostics."""
from __future__ import annotations
import numpy as np
from pbuf.excitation.native_vector_pair_dynamics import positive_relations


def state_validity(displacement: np.ndarray, momentum: np.ndarray) -> dict:
    rel = positive_relations(displacement)
    lengths = np.linalg.norm(rel, axis=-1)
    strain = lengths - 1.0
    return {
        "finite": bool(np.isfinite(displacement).all() and np.isfinite(momentum).all()),
        "max_abs_strain": float(np.max(np.abs(strain))),
        "all_bonds_strictly_within_domain": bool(np.all(np.abs(strain) < 1.0)),
        "relation_min_length": float(np.min(lengths)),
        "classification": "VALID" if np.isfinite(displacement).all() and np.isfinite(momentum).all() and np.all(np.abs(strain) < 1.0) else "INVALID",
    }


def support_relation(a: np.ndarray, b: np.ndarray) -> str:
    overlap = np.asarray(a) & np.asarray(b)
    if not np.any(overlap):
        return "DISJOINT"
    if np.array_equal(a, b):
        return "DIRECT_OVERLAP"
    return "PARTIAL_OVERLAP"
