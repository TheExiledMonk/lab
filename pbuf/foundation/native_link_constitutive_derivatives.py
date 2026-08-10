"""Coefficient-free longitudinal descriptors for the frozen bounded-strain law."""
from __future__ import annotations

import numpy as np


def energy(strain, K=1.0, epsilon_max=1.0):
    e = np.asarray(strain, dtype=float)
    return -0.5 * K * epsilon_max**2 * np.log1p(-(e / epsilon_max) ** 2)


def stress(strain, K=1.0, epsilon_max=1.0):
    e = np.asarray(strain, dtype=float)
    return K * e / (1.0 - (e / epsilon_max) ** 2)


def tangent_stiffness(strain, K=1.0, epsilon_max=1.0):
    """Exact d sigma / d epsilon."""
    e = np.asarray(strain, dtype=float)
    q = (e / epsilon_max) ** 2
    return K * (1.0 + q) / (1.0 - q) ** 2


def constitutive_curvature(strain, K=1.0, epsilon_max=1.0):
    """Exact d2 W / d epsilon2; identical to tangent_stiffness."""
    return tangent_stiffness(strain, K, epsilon_max)


def link_stretch_ratio(separation, reference_separation):
    return np.asarray(separation, dtype=float) / np.asarray(reference_separation, dtype=float)


def strain_from_separation(separation, reference_separation):
    return link_stretch_ratio(separation, reference_separation) - 1.0


def validate_analytic_derivatives(strains=None, K=1.0, epsilon_max=1.0, step=1e-6):
    if strains is None:
        strains = np.linspace(-0.9, 0.9, 181)
    e = np.asarray(strains, dtype=float)
    fd_sigma = (stress(e + step, K, epsilon_max) - stress(e - step, K, epsilon_max)) / (2 * step)
    fd_W2 = (energy(e + step, K, epsilon_max) - 2 * energy(e, K, epsilon_max) + energy(e - step, K, epsilon_max)) / step**2
    exact = tangent_stiffness(e, K, epsilon_max)
    return {
        "tangent_max_abs_error": float(np.max(np.abs(fd_sigma - exact))),
        "curvature_max_abs_error": float(np.max(np.abs(fd_W2 - exact))),
        "tangent_valid": bool(np.allclose(fd_sigma, exact, rtol=2e-8, atol=2e-8)),
        "curvature_valid": bool(np.allclose(fd_W2, exact, rtol=5e-4, atol=5e-4)),
        "equivalent": bool(np.array_equal(tangent_stiffness(e, K, epsilon_max), constitutive_curvature(e, K, epsilon_max))),
    }


def descriptors(strain, reference_separation=1.0):
    e = np.asarray(strain, dtype=float)
    r0 = float(reference_separation)
    r = r0 * (1.0 + e)
    return {"strain": e, "separation": r, "link_stretch_ratio": link_stretch_ratio(r, r0),
            "stress": stress(e), "energy_density": energy(e),
            "tangent_stiffness": tangent_stiffness(e), "constitutive_curvature": constitutive_curvature(e)}
