"""Exact DEV167 tangent objects on an arbitrary loaded N6 state.

This module deliberately contains no prestress parameter: every coefficient is
read directly from the supplied native displacement field.
"""
from __future__ import annotations
import numpy as np

from pbuf.excitation.native_vector_pair_dynamics import positive_relations


def stress_derivative(extension: np.ndarray) -> np.ndarray:
    """d[e/(1-e**2)]/de, preserving the frozen DEV167 normalization."""
    e = np.asarray(extension, dtype=float)
    return (1.0 + e * e) / (1.0 - e * e) ** 2


def bond_tangent(displacement: np.ndarray) -> dict[str, np.ndarray]:
    """Positive-bond parallel/geometric stiffness and 3x3 force Jacobian."""
    r = positive_relations(displacement)
    length = np.linalg.norm(r, axis=-1)
    extension = length - 1.0
    unit = r / length[..., None]
    stress = extension / (1.0 - extension * extension)
    k_parallel = stress_derivative(extension)
    k_perp = stress / length
    eye = np.eye(3)
    nn = unit[..., :, None] * unit[..., None, :]
    matrix = k_parallel[..., None, None] * nn + k_perp[..., None, None] * (eye - nn)
    return {"relation": r, "length": length, "extension": extension, "stress": stress,
            "unit": unit, "k_parallel": k_parallel, "k_perp": k_perp, "matrix": matrix}


def tangent_net_force(displacement: np.ndarray, delta_displacement: np.ndarray) -> np.ndarray:
    """Exact Frechet derivative of ``net_force`` at a loaded state."""
    b = bond_tangent(displacement)
    du = np.asarray(delta_displacement, dtype=float)
    out = np.zeros_like(du)
    for axis in range(3):
        drel = np.roll(du, -1, axis=axis) - du
        df = np.einsum("...ij,...j->...i", b["matrix"][..., axis, :, :], drel)
        out += df - np.roll(df, 1, axis=axis)
    return out


def local_node_tangent(displacement: np.ndarray) -> np.ndarray:
    """Diagonal 3x3 displacement block of the exact N6 tangent at every node."""
    b = bond_tangent(displacement)["matrix"]
    out = -np.sum(b, axis=-3)
    for axis in range(3):
        out -= np.roll(b[..., axis, :, :], 1, axis=axis)
    return out


def tangent_step(displacement: np.ndarray, delta_u: np.ndarray, delta_p: np.ndarray,
                 numerical_step: float = 0.04) -> tuple[np.ndarray, np.ndarray]:
    """Derivative of the frozen kick-drift map along a supplied trajectory state."""
    dp = delta_p + numerical_step * tangent_net_force(displacement, delta_u)
    return delta_u + numerical_step * dp, dp
