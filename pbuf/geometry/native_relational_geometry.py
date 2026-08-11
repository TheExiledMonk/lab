"""Coefficient-free information audit for scalar states on a periodic N6 graph.

This module deliberately separates an oriented *scalar difference* from a
spatial bond vector.  The former follows exactly from a node scalar; the latter
requires a constitutive geometry law which Dev156/159/162 do not contain.
"""
from __future__ import annotations

import numpy as np

from pbuf.excitation.native_bond_state import N6_OFFSETS, relational_differences


def directional_scalar_information(q: np.ndarray) -> np.ndarray:
    """Return the exact (+z,-z,+y,-y,+x,-x) neighbor differences."""
    return relational_differences(np.asarray(q, dtype=float))


def scalar_asymmetry(q: np.ndarray) -> np.ndarray:
    """Three axis-wise scalar asymmetries, in array order (z,y,x)."""
    d = directional_scalar_information(q)
    return np.stack((d[..., 0] - d[..., 1], d[..., 2] - d[..., 3],
                     d[..., 4] - d[..., 5]), axis=-1)


def information_contract(q: np.ndarray) -> dict:
    q = np.asarray(q, dtype=float)
    d = directional_scalar_information(q)
    return {
        "node_state_rank": "SCALAR",
        "node_state_shape": list(q.shape),
        "n6_offsets_array_order_zyx": [list(x) for x in N6_OFFSETS],
        "directed_neighbor_differences_exactly_reconstructible": True,
        "directional_scalar_information_nonzero": bool(np.any(d != 0)),
        "directed_difference_peak": float(np.max(np.abs(d))),
        "spatial_vector_components_stored": False,
        "bond_length_semantics_present": False,
        "bond_orientation_deformation_semantics_present": False,
        "scalar_difference_is_not_assumed_to_be_length": True,
    }


def undeformed_bond_vectors(shape: tuple[int, int, int]) -> np.ndarray:
    """Cartesian N6 control vectors at every node, in native cell units."""
    offsets = np.asarray(N6_OFFSETS, dtype=float)
    return np.broadcast_to(offsets, tuple(shape) + offsets.shape).copy()


def reciprocity_error_for_scalar_edges(q: np.ndarray) -> float:
    """Reciprocity of oriented scalar edge differences (not spatial vectors)."""
    d = directional_scalar_information(q)
    errors = []
    for plus, minus, axis in ((0, 1, 0), (2, 3, 1), (4, 5, 2)):
        errors.append(np.max(np.abs(d[..., plus] + np.roll(d[..., minus], -1, axis=axis))))
    return float(max(errors))
