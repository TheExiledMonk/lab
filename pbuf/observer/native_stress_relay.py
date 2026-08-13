"""Read-only N6 receiver diagnostics used by DEV209.

The functions materialise existing DEV167 relations, forces, and power; they
do not alter the state or introduce a receiver object.
"""
from __future__ import annotations

import numpy as np

from pbuf.excitation.native_vector_pair_dynamics import directed_relations, pair_forces, pair_power_flux
from pbuf.observer.native_stress_response import finite_step_force_response


def directed_forces(displacement: np.ndarray) -> np.ndarray:
    """DEV167 positive pair forces in the canonical directed-N6 ordering."""
    positive = pair_forces(displacement)
    rows = []
    for axis in range(3):
        plus = positive[..., axis, :]
        rows.extend((plus, -np.roll(plus, 1, axis=axis)))
    return np.stack(rows, axis=-2)


def receiver_state(displacement: np.ndarray, momentum: np.ndarray, mask: np.ndarray) -> dict[str, np.ndarray]:
    """Return native state on a fixed bookkeeping mask, never a selected mask."""
    relation = directed_relations(displacement)[mask]
    force = directed_forces(displacement)[mask]
    return {"relation": relation, "stress": force, "momentum": momentum[mask],
            "power_flux": pair_power_flux(displacement, momentum)[mask]}


def finite_components(previous_relation: np.ndarray, next_relation: np.ndarray) -> dict[str, np.ndarray]:
    """The unchanged DEV204 split for a directed N6 bond sequence."""
    l0 = np.linalg.norm(previous_relation, axis=-1)
    l1 = np.linalg.norm(next_relation, axis=-1)
    return finite_step_force_response(l0 - 1.0, l1 - 1.0,
                                      previous_relation / l0[..., None],
                                      next_relation / l1[..., None])
