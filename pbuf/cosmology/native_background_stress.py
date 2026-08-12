"""Frozen DEV167 constitutive response under a uniform scale deformation."""
from __future__ import annotations

import numpy as np
from .native_background_strain import extension


def stress_from_extension(epsilon: float) -> float:
    if not np.isfinite(epsilon) or abs(epsilon) >= 1:
        raise ValueError("requires finite |epsilon| < 1")
    return float(epsilon / (1.0 - epsilon * epsilon))


def stress_derivative(epsilon: float) -> float:
    if not np.isfinite(epsilon) or abs(epsilon) >= 1:
        raise ValueError("requires finite |epsilon| < 1")
    return float((1.0 + epsilon * epsilon) / (1.0 - epsilon * epsilon) ** 2)


def homogeneous_stress(scale: float) -> float:
    return stress_from_extension(extension(scale))


def homogeneous_potential(scale: float) -> float:
    e = extension(scale)
    if not np.isfinite(e) or abs(e) >= 1:
        raise ValueError("requires finite |epsilon| < 1")
    return float(-0.5 * np.log1p(-e * e))


def restoring_generalized_force(scale: float) -> float:
    """Negative derivative of one-bond stored potential with uniform scale."""
    return -homogeneous_stress(scale)
