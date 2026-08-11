"""Isolated reuse audit of bounded static stress as a dynamic restoring force."""
from __future__ import annotations

import numpy as np

from pbuf.excitation.native_bond_state import gradient_adjoint, positive_gradient
from pbuf.excitation.native_relational_dynamics import N6_COORDINATION
from pbuf.wl.native_incremental_elastic_energy import bounded_strain_stress


def bounded_f03_step(q: np.ndarray, retained: np.ndarray, K: float = 1.0,
                     excursion_max: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
    """Kick-drift candidate using the frozen sigma(xi), without promotion."""
    q = np.asarray(q, float); retained = np.asarray(retained, float)
    restoring = bounded_strain_stress(positive_gradient(q), K, excursion_max)
    r1 = retained - gradient_adjoint(restoring) / N6_COORDINATION
    return q + r1, r1


def bounded_f03_inverse(q1: np.ndarray, r1: np.ndarray, K: float = 1.0,
                        excursion_max: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
    """Exact algebraic inverse while the reconstructed excursion is admissible."""
    q0 = np.asarray(q1, float) - np.asarray(r1, float)
    restoring = bounded_strain_stress(positive_gradient(q0), K, excursion_max)
    return q0, np.asarray(r1, float) + gradient_adjoint(restoring) / N6_COORDINATION


def bounded_response_run(q: np.ndarray, steps: int = 12,
                         excursion_max: float = 1.0) -> dict:
    q0 = np.asarray(q, float).copy(); r0 = np.zeros_like(q0)
    qn, rn = q0.copy(), r0.copy()
    peak = float(np.max(np.abs(positive_gradient(qn))))
    stable = True
    try:
        for _ in range(int(steps)):
            qn, rn = bounded_f03_step(qn, rn, excursion_max=excursion_max)
            peak = max(peak, float(np.max(np.abs(positive_gradient(qn)))))
        qr, rr = qn.copy(), rn.copy()
        for _ in range(int(steps)):
            qr, rr = bounded_f03_inverse(qr, rr, excursion_max=excursion_max)
        reverse_error = max(float(np.max(np.abs(qr-q0))), float(np.max(np.abs(rr-r0))))
    except ValueError:
        stable = False; reverse_error = None
    return {"steps": int(steps), "stable_within_bound": stable,
            "maximum_excursion": peak, "exact_reverse_error": reverse_error,
            "reversible": bool(reverse_error is not None and reverse_error < 1e-11),
            "bound_enforced_by_update": False,
            "conserved_nonlinear_functional_derived": False}
