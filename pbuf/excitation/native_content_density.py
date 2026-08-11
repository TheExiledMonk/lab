"""Dev158 audit of local decompositions of the frozen global invariants."""
from __future__ import annotations

import numpy as np

from pbuf.excitation.native_bond_state import positive_gradient
from pbuf.excitation.native_relational_dynamics import N6_COORDINATION


def signed_candidate_density(representation: str, q, memory) -> np.ndarray:
    """Exact summand allocation; intentionally not asserted pointwise positive."""
    q = np.asarray(q, float); m = np.asarray(memory, float); g = positive_gradient(q)
    if representation == "F02":
        return q*q + np.sum(m*m + g*m, axis=-1) / N6_COORDINATION
    if representation == "F03":
        return m*m + np.sum(g*g-g*positive_gradient(m), axis=-1) / N6_COORDINATION
    raise ValueError("representation must be F02 or F03")


def positivity_audit(representation: str, q, memory) -> dict:
    h = signed_candidate_density(representation, q, memory)
    return {"representation": representation, "sum": float(h.sum()),
            "minimum": float(h.min()), "negative_site_count": int(np.count_nonzero(h < 0)),
            "exact_local_summand": True, "pointwise_positive_proven": False,
            "classification": "NOT_DERIVED"}
