"""Exact finite-step native bond-stress response for DEV204.

This is a read-only decomposition of the frozen DEV167 force law.  It does
not define a field, add a degree of freedom, or approximate a time derivative.
"""
from __future__ import annotations

import numpy as np


def finite_step_force_response(strain: np.ndarray, next_strain: np.ndarray,
                               unit: np.ndarray, next_unit: np.ndarray) -> dict[str, np.ndarray]:
    """Split ``F1-F0`` exactly into magnitude, direction, and cross terms.

    With ``F=sigma(e) u``, ``ds=sigma(e1)-sigma(e0)`` and ``du=u1-u0``:
    ``dF = ds*u + sigma(e)*du + ds*du``.  The final term is retained because
    the production trajectory is finite-step and nonlinear.
    """
    e0, e1 = np.asarray(strain), np.asarray(next_strain)
    u0, u1 = np.asarray(unit), np.asarray(next_unit)
    s0 = e0 / (1.0 - e0 * e0)
    s1 = e1 / (1.0 - e1 * e1)
    delta_sigma = s1 - s0
    delta_unit = u1 - u0
    magnitude = delta_sigma[..., None] * u0
    orientation = s0[..., None] * delta_unit
    cross = delta_sigma[..., None] * delta_unit
    total = magnitude + orientation + cross
    return {"delta_force": total, "magnitude": magnitude,
            "orientation": orientation, "cross": cross,
            "delta_sigma": delta_sigma, "delta_unit": delta_unit}
