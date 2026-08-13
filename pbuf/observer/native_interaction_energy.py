"""DEV211 read-only native bond energy and interaction helpers."""
from __future__ import annotations

import numpy as np

from pbuf.excitation.native_vector_pair_dynamics import bounded_stress, positive_relations


def bond_state(displacement: np.ndarray) -> dict[str, np.ndarray]:
    """Positive N6 bond geometry, DEV167 force, and DEV167 potential."""
    relation = positive_relations(displacement)
    length = np.linalg.norm(relation, axis=-1)
    extension = length - 1.0
    stress = bounded_stress(extension)
    force = stress[..., None] * relation / length[..., None]
    energy = -0.5 * np.log1p(-extension * extension)
    return {"relation": relation, "length": length, "extension": extension,
            "stress": stress, "force": force, "energy": energy}


def total_potential(displacement: np.ndarray) -> float:
    return float(np.sum(bond_state(displacement)["energy"]))


def node_force_from_positive_bonds(force: np.ndarray) -> np.ndarray:
    """Sum reciprocal positive-bond forces into the force on every node."""
    out = np.zeros(force.shape[:-2] + (3,), dtype=np.float64)
    for axis in range(3):
        out += force[..., axis, :] - np.roll(force[..., axis, :], 1, axis=axis)
    return out


def interaction_residual(ab: np.ndarray, a: np.ndarray, b: np.ndarray,
                         quiet: np.ndarray) -> np.ndarray:
    return ab - a - b + quiet
