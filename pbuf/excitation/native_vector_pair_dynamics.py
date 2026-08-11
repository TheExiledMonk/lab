"""DEV167 native vector relations and distance-bound central-pair dynamics.

Coordinates are dimensionless lattice relations and ``step`` is a native
progression index.  Nothing in this module assigns an SI length or time.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

POSITIVE_DIRECTIONS = np.eye(3, dtype=np.float64)
N6_OFFSETS = ((1, 0, 0), (-1, 0, 0), (0, 1, 0),
              (0, -1, 0), (0, 0, 1), (0, 0, -1))


@dataclass
class VectorPairState:
    """Integrable N6 relation state represented by node displacements."""

    displacement: np.ndarray
    momentum: np.ndarray
    progression_step: int = 0

    def __post_init__(self) -> None:
        self.displacement = np.asarray(self.displacement, dtype=np.float64)
        self.momentum = np.asarray(self.momentum, dtype=np.float64)
        if self.displacement.ndim != 4 or self.displacement.shape[-1] != 3:
            raise ValueError("displacement must have shape (Nx, Ny, Nz, 3)")
        if min(self.displacement.shape[:3]) < 3:
            raise ValueError("all three periodic N6 dimensions must be present")
        if self.momentum.shape != self.displacement.shape:
            raise ValueError("momentum must match displacement")
        if not np.isfinite(self.displacement).all() or not np.isfinite(self.momentum).all():
            raise ValueError("state must be finite")


def positive_relations(displacement: np.ndarray) -> np.ndarray:
    """Return r[a,+axis]=e_axis+x[a+axis]-x[a], shape (...,3,3)."""
    u = np.asarray(displacement, dtype=np.float64)
    return np.stack([POSITIVE_DIRECTIONS[i] + np.roll(u, -1, axis=i) - u
                     for i in range(3)], axis=-2)


def directed_relations(displacement: np.ndarray) -> np.ndarray:
    """Materialize the six oriented relations in canonical N6 order."""
    rp = positive_relations(displacement)
    out = []
    for axis in range(3):
        plus = rp[..., axis, :]
        minus = -np.roll(plus, 1, axis=axis)
        out.extend((plus, minus))
    return np.stack(out, axis=-2)


def relation_antisymmetry_error(displacement: np.ndarray) -> float:
    r = directed_relations(displacement)
    errors = [np.max(np.abs(r[..., 2*a, :] + np.roll(r[..., 2*a+1, :], -1, axis=a)))
              for a in range(3)]
    return float(max(errors))


def bounded_stress(extension: np.ndarray) -> np.ndarray:
    """Frozen K=epsilon_max=1 constitutive magnitude."""
    e = np.asarray(extension, dtype=np.float64)
    if np.any(np.abs(e) >= 1.0):
        raise ValueError("bounded-strain relation requires |epsilon| < 1")
    return e / (1.0 - e*e)


def pair_forces(displacement: np.ndarray) -> np.ndarray:
    """Force on the tail of each positive bond; stretched bonds pull inward."""
    r = positive_relations(displacement)
    length = np.linalg.norm(r, axis=-1)
    return bounded_stress(length - 1.0)[..., None] * r / length[..., None]


def net_force(displacement: np.ndarray) -> np.ndarray:
    """Sum reciprocal central forces at every node."""
    fp = pair_forces(displacement)
    force = np.zeros_like(displacement, dtype=np.float64)
    for axis in range(3):
        f = fp[..., axis, :]
        force += f - np.roll(f, 1, axis=axis)
    return force


def pair_reciprocity_error(displacement: np.ndarray) -> float:
    """Reciprocity is exact by construction; audit its materialized form."""
    fp = pair_forces(displacement)
    errors = []
    for axis in range(3):
        reverse = -np.roll(fp[..., axis, :], 1, axis=axis)
        errors.append(np.max(np.abs(fp[..., axis, :] + np.roll(reverse, -1, axis=axis))))
    return float(max(errors))


def potential(displacement: np.ndarray) -> float:
    e = np.linalg.norm(positive_relations(displacement), axis=-1) - 1.0
    if np.any(np.abs(e) >= 1.0):
        return float("inf")
    return float(np.sum(-0.5 * np.log1p(-e*e)))


def invariant(displacement: np.ndarray, momentum: np.ndarray) -> float:
    """Native Hamiltonian candidate (not an SI-energy claim)."""
    return float(0.5*np.sum(np.asarray(momentum)**2) + potential(displacement))


def step(state: VectorPairState, numerical_step: float = 0.05,
         external_force: np.ndarray | None = None) -> VectorPairState:
    """Reversible symplectic kick-drift map; numerical_step is convergence-only."""
    f = net_force(state.displacement)
    if external_force is not None:
        f = f + np.asarray(external_force, dtype=np.float64)
    p1 = state.momentum + numerical_step*f
    u1 = state.displacement + numerical_step*p1
    return VectorPairState(u1, p1, state.progression_step + 1)


def inverse_step(state: VectorPairState, numerical_step: float = 0.05,
                 external_force: np.ndarray | None = None) -> VectorPairState:
    u0 = state.displacement - numerical_step*state.momentum
    f = net_force(u0)
    if external_force is not None:
        f = f + np.asarray(external_force, dtype=np.float64)
    p0 = state.momentum - numerical_step*f
    return VectorPairState(u0, p0, state.progression_step - 1)


def pair_power_flux(displacement: np.ndarray, momentum: np.ndarray) -> np.ndarray:
    """Oriented positive-bond interaction power J_ab=-F_ab·(p_a+p_b)/2.

    The reverse orientation is exactly its negative.  Its node divergence,
    together with local kinetic/potential partition, gives energy continuity.
    """
    fp = pair_forces(displacement)
    p = np.asarray(momentum, dtype=np.float64)
    return np.stack([-np.sum(fp[..., a, :]*(p + np.roll(p, -1, axis=a))/2.0, axis=-1)
                     for a in range(3)], axis=-1)


def source_contact_force(shape: tuple[int, int, int], center: tuple[int, int, int],
                         magnitude: float = 0.02) -> np.ndarray:
    """Repulsive one-cell contact at the six N6 neighbors of a source node."""
    out = np.zeros(shape + (3,), dtype=np.float64)
    for offset in N6_OFFSETS:
        idx = tuple((center[i] + offset[i]) % shape[i] for i in range(3))
        out[idx] += magnitude*np.asarray(offset, dtype=np.float64)
    return out


def relax_source_equilibrium(shape: tuple[int, int, int], center: tuple[int, int, int],
                             magnitude: float = 0.02, tolerance: float = 2e-9,
                             max_iterations: int = 20000) -> tuple[np.ndarray, dict]:
    """Conjugate-gradient static solve; damping is not used in propagation."""
    from scipy.optimize import minimize
    ext = source_contact_force(shape, center, magnitude)
    def unpack(v):
        u = v.reshape(shape + (3,))
        return u - np.mean(u, axis=(0, 1, 2), keepdims=True)
    def objective(v):
        u = unpack(v)
        return potential(u) - float(np.sum(ext*u))
    def jac(v):
        u = unpack(v)
        g = -(net_force(u) + ext)
        g -= np.mean(g, axis=(0, 1, 2), keepdims=True)
        return g.ravel()
    result = minimize(objective, np.zeros(np.prod(shape)*3), jac=jac,
                      method="L-BFGS-B", options={"gtol": tolerance, "maxiter": max_iterations})
    u = unpack(result.x)
    residual = net_force(u) + ext
    return u, {"success": bool(result.success), "iterations": int(result.nit),
               "max_force_residual": float(np.max(np.abs(residual))),
               "message": str(result.message)}
