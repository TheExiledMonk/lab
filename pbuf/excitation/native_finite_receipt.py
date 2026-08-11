"""Finite receipt construction for the frozen DEV167 vector-pair dynamics.

All coordinates are lattice indices and all times are progression steps.  The
module adds no propagation law: it only records node and positive-bond plane
crossings from :mod:`native_vector_pair_dynamics` states.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .native_vector_pair_dynamics import pair_forces, pair_power_flux, positive_relations


@dataclass(frozen=True)
class NativeReceivedState:
    source_positions: np.ndarray
    received_positions: np.ndarray
    directions: np.ndarray
    weights: np.ndarray
    progression_steps: np.ndarray
    native_cell_ids: np.ndarray
    local_displacement: np.ndarray
    local_momentum: np.ndarray
    local_flux: np.ndarray
    local_content_candidates: np.ndarray
    representation: str

    def __post_init__(self) -> None:
        n = len(np.asarray(self.weights))
        for name in ("source_positions", "received_positions", "directions",
                     "local_displacement", "local_momentum", "local_flux"):
            if np.asarray(getattr(self, name)).shape != (n, 3):
                raise ValueError(f"{name} must have shape (N,3)")
        if np.asarray(self.progression_steps).shape != (n,):
            raise ValueError("progression_steps must have shape (N,)")
        if np.asarray(self.local_content_candidates).shape != (n, 4):
            raise ValueError("local_content_candidates must contain W01--W04")
        if np.any(np.asarray(self.weights) < 0) or not np.isfinite(self.weights).all():
            raise ValueError("receipt weights must be finite and nonnegative")

    def arrays(self) -> dict[str, np.ndarray]:
        return {k: np.asarray(getattr(self, k)) for k in (
            "source_positions", "received_positions", "directions", "weights",
            "progression_steps", "native_cell_ids", "local_displacement",
            "local_momentum", "local_flux", "local_content_candidates")}


def local_content_candidates(displacement: np.ndarray, momentum: np.ndarray) -> np.ndarray:
    """Return W01 kinetic, W02 bond share, W03 sum, and W04 flux magnitude."""
    kinetic = .5*np.sum(np.asarray(momentum, float)**2, axis=-1)
    relations = positive_relations(displacement)
    extension = np.linalg.norm(relations, axis=-1)-1.0
    bond = -.5*np.log1p(-extension*extension)
    potential_share = np.zeros_like(kinetic)
    for axis in range(3):
        potential_share += .5*(bond[..., axis] + np.roll(bond[..., axis], 1, axis=axis))
    flux = flux_vectors(displacement, momentum)
    return np.stack((kinetic, potential_share, kinetic+potential_share,
                     np.linalg.norm(flux, axis=-1)), axis=-1)


def flux_vectors(displacement: np.ndarray, momentum: np.ndarray) -> np.ndarray:
    """Vectorize the three positive-bond DEV167 signed pair fluxes."""
    j = pair_power_flux(displacement, momentum)
    directions = positive_relations(displacement)
    directions /= np.linalg.norm(directions, axis=-1, keepdims=True)
    return np.sum(j[..., None]*directions, axis=-2)


def unit_directions(vectors: np.ndarray) -> np.ndarray:
    vectors = np.asarray(vectors, float)
    norms = np.linalg.norm(vectors, axis=-1, keepdims=True)
    return np.divide(vectors, norms, out=np.zeros_like(vectors), where=norms > 0)


def crossing_bond_flux(displacement: np.ndarray, momentum: np.ndarray, plane_x: int) -> np.ndarray:
    """Signed +x flux through the predeclared x=plane_x lattice face."""
    if not 1 <= plane_x < displacement.shape[0]:
        raise ValueError("receipt plane must be an interior positive-x bond face")
    return pair_power_flux(displacement, momentum)[plane_x-1, ..., 0]


def plane_node_snapshot(displacement: np.ndarray, momentum: np.ndarray, plane_x: int):
    """Return complete local state on the node layer immediately after the face."""
    content = local_content_candidates(displacement, momentum)[plane_x]
    flux = flux_vectors(displacement, momentum)[plane_x]
    return displacement[plane_x], momentum[plane_x], flux, content
