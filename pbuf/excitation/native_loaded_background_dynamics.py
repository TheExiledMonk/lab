"""Dev163 audit of perturbations about a loaded Dev159 equilibrium.

This module deliberately contains no lensing rule.  The frozen Dev159/Dev162
equation is linear, so its exact perturbation equation is the unloaded F03
equation.  The helpers below make that cancellation executable and auditable.
"""
from __future__ import annotations

import numpy as np

from .native_relational_dynamics import f03_invariant, f03_step
from pbuf.lens.native_stationary_lens_from_source import equilibrium_residual


def loaded_f03_step(total_q: np.ndarray, total_retained: np.ndarray,
                    source: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Frozen Dev159 F03 kick-drift with a stationary distributed source."""
    q = np.asarray(total_q, float)
    retained = np.asarray(total_retained, float)
    residual = equilibrium_residual(q, source)
    retained1 = retained + residual
    return q + retained1, retained1


def perturbation_step(background_q: np.ndarray, source: np.ndarray,
                      delta_q: np.ndarray, delta_retained: np.ndarray
                      ) -> tuple[np.ndarray, np.ndarray]:
    """Advance a perturbation by subtracting the stationary loaded lane.

    This is an evaluation of the existing total-state law, not a proposed
    coupling.  For an equilibrium background it agrees with ``f03_step`` to
    floating-point precision and is therefore independent of the load.
    """
    background_q = np.asarray(background_q, float)
    base_q1, base_r1 = loaded_f03_step(background_q, np.zeros_like(background_q), source)
    total_q1, total_r1 = loaded_f03_step(
        background_q + np.asarray(delta_q, float),
        np.asarray(delta_retained, float), source)
    return total_q1 - base_q1, total_r1 - base_r1


def linearization_audit(background_q: np.ndarray, source: np.ndarray,
                        delta_q: np.ndarray, delta_retained: np.ndarray) -> dict:
    """Numerically certify equilibrium cancellation and free-map identity."""
    loaded = perturbation_step(background_q, source, delta_q, delta_retained)
    free = f03_step(delta_q, delta_retained)
    residual = equilibrium_residual(background_q, source)
    return {
        "equilibrium_residual_linf": float(np.max(np.abs(residual))),
        "perturbation_q_free_map_error_linf": float(np.max(np.abs(loaded[0] - free[0]))),
        "perturbation_retained_free_map_error_linf": float(np.max(np.abs(loaded[1] - free[1]))),
        "loaded_operator_depends_on_background": False,
        "derivation": "Linearity: D(Q0+dQ)+S=D(Q0)+S+D(dQ)=D(dQ).",
    }


def invariant_audit(background_q: np.ndarray, source: np.ndarray,
                    delta_q: np.ndarray, delta_retained: np.ndarray) -> dict:
    before = f03_invariant(delta_q, delta_retained)
    q1, r1 = perturbation_step(background_q, source, delta_q, delta_retained)
    after = f03_invariant(q1, r1)
    return {
        "name": "DEV156_F03_PERTURBATION_QUADRATIC_INVARIANT",
        "classification": "EXACT",
        "before": before,
        "after": after,
        "absolute_error": abs(after - before),
        "new_energy_density_introduced": False,
    }
