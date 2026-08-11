"""Dev158 coefficient-free bridge diagnostics for static and dynamic excursions.

This module is experimental.  It maps existing state variables but deliberately
does not add a static/dynamic coupling or alter either frozen evolution law.
"""
from __future__ import annotations

import numpy as np

from pbuf.excitation.native_bond_state import positive_gradient
from pbuf.excitation.native_relational_dynamics import (
    N6_COORDINATION, f02_invariant, f03_invariant,
)
from pbuf.wl.native_incremental_elastic_energy import bounded_strain_energy


def native_bond_excursion(q: np.ndarray, dx: float = 1.0) -> np.ndarray:
    """Existing positive-axis relational difference, optionally as strain."""
    dx = float(dx)
    if not np.isfinite(dx) or dx <= 0:
        raise ValueError("dx must be finite and positive")
    return positive_gradient(np.asarray(q, float)) / dx


def static_small_excursion(strain, K: float = 1.0,
                           epsilon_max: float = 1.0) -> dict:
    """Compare the exact frozen W with its law-derived quadratic term."""
    e = np.asarray(strain, float)
    exact = bounded_strain_energy(e, K, epsilon_max)
    leading = 0.5 * float(K) * e * e
    remainder = exact - leading
    return {"exact": exact, "quadratic": leading, "remainder": remainder,
            "max_absolute_remainder": float(np.max(np.abs(remainder))),
            "series": "K*epsilon^2/2 + K*epsilon^4/(4*epsilon_max^2) + O(epsilon^6)"}


def dynamic_content(representation: str, q: np.ndarray, memory: np.ndarray) -> float:
    if representation == "F02":
        return f02_invariant(q, memory)
    if representation == "F03":
        return f03_invariant(q, memory)
    raise ValueError("representation must be F02 or F03")


def common_mapping_contract() -> dict:
    return {
        "common_variable": "xi_ab=(q_b-q_a)/Delta_x; Delta_x=1 in Dev156 native cells",
        "static_mapping": "xi_ab=epsilon_ab=(u_b-u_a)/Delta_x",
        "dynamic_mapping": "xi_ab=positive_gradient(q)/Delta_x",
        "status": "DERIVED",
        "normalization_origin": "existing static strain definition and native-cell Dev156 gradient",
        "new_physical_coefficient": False,
        "quadratic_match": "STRUCTURAL",
        "reason_not_exact": "Dev156 exact invariants include auxiliary-memory and time-stagger cross terms",
        "coordination": N6_COORDINATION,
    }
