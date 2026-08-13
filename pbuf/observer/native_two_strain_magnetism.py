"""Static two-structure diagnostics built exclusively from frozen DEV167 bonds."""
from __future__ import annotations

import numpy as np

from pbuf.excitation.native_vector_pair_dynamics import net_force, potential
from pbuf.observer.native_interaction_energy import bond_state, interaction_residual, node_force_from_positive_bonds


def source_contact(shape: tuple[int, int, int], centers: tuple[tuple[int, int, int], ...], magnitude: float) -> np.ndarray:
    """Exact translated copies of DEV167's six-neighbour source contact."""
    from pbuf.excitation.native_vector_pair_dynamics import source_contact_force
    return sum((source_contact_force(shape, c, magnitude) for c in centers), np.zeros(shape + (3,)))


def relax_contacts(shape: tuple[int, int, int], centers: tuple[tuple[int, int, int], ...],
                   magnitude: float, tolerance: float = 2e-9) -> tuple[np.ndarray, dict]:
    """Static equilibrium under existing DEV167 contact loading; no new profile."""
    from scipy.optimize import minimize
    ext = source_contact(shape, centers, magnitude)
    def unpack(v):
        u = v.reshape(shape + (3,))
        return u - np.mean(u, axis=(0, 1, 2), keepdims=True)
    def objective(v):
        u = unpack(v)
        return potential(u) - float(np.sum(ext * u))
    def jac(v):
        u = unpack(v)
        g = -(net_force(u) + ext)
        g -= np.mean(g, axis=(0, 1, 2), keepdims=True)
        return g.ravel()
    result = minimize(objective, np.zeros(np.prod(shape) * 3), jac=jac, method="L-BFGS-B",
                      options={"gtol": tolerance, "maxiter": 20000})
    u = unpack(result.x)
    return u, {"success": bool(result.success), "iterations": int(result.nit),
               "max_force_residual": float(np.max(np.abs(net_force(u) + ext))),
               "message": str(result.message)}


def support_mask(shape: tuple[int, int, int], center: tuple[int, int, int]) -> np.ndarray:
    """Predeclared DEV167 source-contact support: its exact six neighbour nodes."""
    m = np.zeros(shape, dtype=bool)
    for off in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)):
        m[tuple((center[i] + off[i]) % shape[i] for i in range(3))] = True
    return m


def state_fields(displacement: np.ndarray) -> dict[str, np.ndarray | float]:
    bonds = bond_state(displacement)
    return {**bonds, "node_force": node_force_from_positive_bonds(bonds["force"]),
            "potential": float(np.sum(bonds["energy"]))}


def sum_support(force: np.ndarray, mask: np.ndarray) -> np.ndarray:
    return np.sum(force[mask], axis=0)


def torque_support(force: np.ndarray, mask: np.ndarray, center: tuple[int, int, int]) -> np.ndarray:
    coords = np.argwhere(mask).astype(np.float64)
    r = coords - np.asarray(center, dtype=np.float64)
    # Supports stay clear of periodic boundaries in the fixed DEV211 geometry.
    return np.sum(np.cross(r, force[mask]), axis=0)


def interaction_fields(ab: dict, a: dict, b: dict, quiet: dict) -> dict:
    return {key: interaction_residual(ab[key], a[key], b[key], quiet[key])
            for key in ("node_force", "force", "energy", "extension", "stress", "relation")}
