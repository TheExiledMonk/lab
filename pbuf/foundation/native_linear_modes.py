"""Exact tangent spectrum of the unloaded DEV167 N6 pair dynamics."""
from __future__ import annotations

import numpy as np


def allowed_wavevectors(shape: tuple[int, int, int]) -> tuple[np.ndarray, np.ndarray]:
    """Return integer FFT indices and their exact periodic wavevectors."""
    axes = [np.fft.fftfreq(n, d=1.0) * (2.0 * np.pi) for n in shape]
    inds = [np.rint(np.fft.fftfreq(n) * n).astype(int) for n in shape]
    grid_i = np.stack(np.meshgrid(*inds, indexing="ij"), axis=-1).reshape(-1, 3)
    grid_k = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, 3)
    return grid_i, grid_k


def tangent_symbol(k: np.ndarray) -> np.ndarray:
    """DEV167 unloaded force symbol L(k)=diag(-4 sin²(k_i/2)).

    This follows by differentiating the central bond force at unit extension:
    only its bond-parallel component survives at first order.
    """
    k = np.asarray(k, dtype=float)
    return np.diag(-4.0 * np.sin(k / 2.0) ** 2)


def update_matrix(k: np.ndarray, numerical_step: float) -> np.ndarray:
    """Six-dimensional kick-drift tangent update in (U,Q) order."""
    L = tangent_symbol(k)
    eye = np.eye(3)
    h = float(numerical_step)
    return np.block([[eye + h * h * L, h * eye], [h * L, eye]])


def oscillatory_frequency(eigenvalue: complex, numerical_step: float) -> float | None:
    """Principal native frequency only for unit-modulus nonzero eigenvalues."""
    if not np.isclose(abs(eigenvalue), 1.0, rtol=0.0, atol=1e-12):
        return None
    angle = abs(float(np.angle(eigenvalue)))
    return angle / float(numerical_step) if angle > 1e-13 else 0.0


def small_k_continuum_statement() -> dict:
    return {
        "equation": "partial_t^2 u_i = partial_i^2 u_i + O(partial_i^4, h^2)",
        "coefficient_origin": "unit tangent stiffness of the frozen DEV167 bounded stress at epsilon=0",
        "not_imported_continuum_elasticity": True,
    }
