"""Coefficient-free local N6 field representation for frozen DEV167 states."""
from __future__ import annotations
import numpy as np
from pbuf.excitation.native_vector_pair_dynamics import directed_relations, bounded_stress


def n6_field(displacement: np.ndarray, momentum: np.ndarray) -> dict[str, np.ndarray]:
    """Return the ordered (+x,-x,+y,-y,+z,-z) local dynamical state."""
    relation = directed_relations(displacement)
    length = np.linalg.norm(relation, axis=-1)
    strain = length - 1.0
    unit = relation / length[..., None]
    force = bounded_stress(strain)[..., None] * unit
    return {"strain": strain, "relation": relation, "unit": unit,
            "force": force, "momentum": np.asarray(momentum, dtype=np.float64)}


def symmetric_antisymmetric(force: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Exact opposite-bond sums/differences, axis order x,y,z."""
    f = np.asarray(force)
    return f[..., 0::2, :] + f[..., 1::2, :], f[..., 0::2, :] - f[..., 1::2, :]
