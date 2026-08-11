"""Coefficient-free Dev156 candidate dynamics (experimental, unpromoted)."""
from __future__ import annotations
import numpy as np
from .native_bond_state import positive_gradient, gradient_adjoint, relational_imbalance

N6_COORDINATION = 6


def f01_step(q: np.ndarray) -> np.ndarray:
    """First-order relational relaxation: the N6 neighbor mean."""
    return np.asarray(q, float) + relational_imbalance(q) / N6_COORDINATION


def f03_step(q: np.ndarray, retained: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Reversible kick-drift map driven only by normalized N6 imbalance."""
    r1 = np.asarray(retained, float) + relational_imbalance(q) / N6_COORDINATION
    return np.asarray(q, float) + r1, r1


def f03_inverse(q1: np.ndarray, r1: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    q0 = np.asarray(q1, float) - np.asarray(r1, float)
    return q0, np.asarray(r1, float) - relational_imbalance(q0) / N6_COORDINATION


def f02_step(q: np.ndarray, bonds: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Bond-storage kick-drift candidate; 1/6 is fixed by N6 coordination."""
    b1 = np.asarray(bonds, float) + positive_gradient(q)
    return np.asarray(q, float) - gradient_adjoint(b1) / N6_COORDINATION, b1


def f02_inverse(q1: np.ndarray, b1: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    q0 = np.asarray(q1, float) + gradient_adjoint(b1) / N6_COORDINATION
    return q0, np.asarray(b1, float) - positive_gradient(q0)


def f02_invariant(q: np.ndarray, bonds: np.ndarray) -> float:
    """Exact node-plus-bond quadratic invariant of the F02 map."""
    q = np.asarray(q, float); b = np.asarray(bonds, float)
    return float(np.sum(q*q) + np.sum(b*b)/N6_COORDINATION
                 + np.sum(positive_gradient(q)*b)/N6_COORDINATION)


def f04_step(q: np.ndarray, retained: np.ndarray, bonds: np.ndarray):
    """Inventory candidate with both memories; no claim of minimality."""
    b1 = np.asarray(bonds, float) + positive_gradient(q)
    r1 = np.asarray(retained, float) - gradient_adjoint(b1) / N6_COORDINATION
    return np.asarray(q, float) + r1, r1, b1


def f04_inverse(q1: np.ndarray, r1: np.ndarray, b1: np.ndarray):
    q0 = np.asarray(q1, float) - np.asarray(r1, float)
    r0 = np.asarray(r1, float) + gradient_adjoint(b1) / N6_COORDINATION
    b0 = np.asarray(b1, float) - positive_gradient(q0)
    return q0, r0, b0


def f03_invariant(q: np.ndarray, retained: np.ndarray) -> float:
    """Exact quadratic invariant derived from the kick-drift update."""
    r = np.asarray(retained, float)
    g = positive_gradient(q)
    # For r'=r-Kq, q'=q+r' with K=G*G/6, this is exactly invariant.
    return float(np.sum(r*r) + np.sum(g*g)/N6_COORDINATION
                 - np.sum(g*positive_gradient(r))/N6_COORDINATION)
