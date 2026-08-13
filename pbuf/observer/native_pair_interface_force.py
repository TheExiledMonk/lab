"""Exact N6 interface-bond inventory and reciprocal transfer observer."""
from __future__ import annotations

import numpy as np
from pbuf.excitation.native_vector_pair_dynamics import pair_forces


def interface_bonds(omega_a: np.ndarray, omega_b: np.ndarray) -> dict[str, np.ndarray]:
    """Canonical positive-orientation inventory of all direct A/B N6 edges."""
    rows = []
    for axis in range(3):
        forward_b = np.roll(omega_b, -1, axis=axis)
        forward_a = np.roll(omega_a, -1, axis=axis)
        for node in np.argwhere(omega_a & forward_b):
            rows.append((node, (node + np.eye(3, dtype=int)[axis]) % omega_a.shape, axis, 1))
        # Canonical orientation remains A -> B even when physical positive bond is B -> A.
        for node in np.argwhere(omega_b & forward_a):
            rows.append(((node + np.eye(3, dtype=int)[axis]) % omega_a.shape, node, axis, -1))
    if not rows:
        return {k: np.empty((0, 3), dtype=int) if k in ('node_a', 'node_b') else np.empty(0, dtype=int)
                for k in ('node_a', 'node_b', 'axis', 'orientation')}
    return {'node_a': np.asarray([r[0] for r in rows], dtype=int),
            'node_b': np.asarray([r[1] for r in rows], dtype=int),
            'axis': np.asarray([r[2] for r in rows], dtype=int),
            'orientation': np.asarray([r[3] for r in rows], dtype=int)}


def transfer(displacement: np.ndarray, inventory: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    """Return F_A<-B and its independently materialized reciprocal F_B<-A."""
    fp = pair_forces(displacement)
    ab = np.zeros(3); ba = np.zeros(3)
    for a, b, axis, orientation in zip(inventory['node_a'], inventory['node_b'], inventory['axis'], inventory['orientation']):
        # pair_forces is force on positive-bond tail.  Multiply to orient A -> B.
        f_a_to_b = fp[tuple(a)][axis] if orientation == 1 else -fp[tuple(b)][axis]
        ab += f_a_to_b
        ba -= f_a_to_b
    return ab, ba


def bond_transfer(displacement: np.ndarray, inventory: dict[str, np.ndarray]) -> np.ndarray:
    """Return the force on A from B for every stored canonical A->B bond.

    The order is deliberately the order in the persisted DEV217 inventory.  It
    is therefore suitable for audits which must prove that no interface bond
    was selected, reordered, or reconstructed from a later force result.
    """
    fp = pair_forces(displacement)
    rows = []
    for a, b, axis, orientation in zip(inventory['node_a'], inventory['node_b'],
                                       inventory['axis'], inventory['orientation']):
        rows.append(fp[tuple(a)][axis] if orientation == 1 else -fp[tuple(b)][axis])
    return np.asarray(rows, dtype=float).reshape((-1, 3))
