"""Frozen DEV226 N6 nearest-neighbour tensor-contraction audit."""
from __future__ import annotations

import numpy as np

from pbuf.analysis.native_handedness_representation import axial_dual, frobenius_neighbor_relation


def unique_n6_bonds(shape: tuple[int, int, int]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """One periodic positive-direction bond per node and primitive axis."""
    nodes = np.indices(shape).reshape(3, -1).T
    pairs, axes = [], []
    for axis in range(3):
        neighbour = nodes.copy()
        neighbour[:, axis] = (neighbour[:, axis] + 1) % shape[axis]
        pairs.append(np.stack((nodes, neighbour), axis=1))
        axes.append(np.full(len(nodes), axis, dtype=np.int8))
    return np.concatenate(pairs), np.concatenate(axes), np.arange(3, dtype=np.int8)


def contract_bonds(tensor: np.ndarray, pairs: np.ndarray) -> np.ndarray:
    """Return raw C_ab=A(a):A(b), with no normalization or tolerance."""
    a = tensor[:, pairs[:, 0, 0], pairs[:, 0, 1], pairs[:, 0, 2]]
    b = tensor[:, pairs[:, 1, 0], pairs[:, 1, 1], pairs[:, 1, 2]]
    return frobenius_neighbor_relation(a, b)


def sign_counts(values: np.ndarray, axes: np.ndarray) -> dict:
    """Exact machine sign classes and axis-resolved counts."""
    negative, positive, zero = values < 0, values > 0, values == 0
    counts = {"opposed": negative.sum(1), "aligned": positive.sum(1), "zero": zero.sum(1)}
    by_axis = {}
    for index, name in enumerate("xyz"):
        mask = axes == index
        by_axis[name] = {key: value[:, mask].sum(1) for key, value in
                         (("opposed", negative), ("aligned", positive), ("zero", zero))}
        by_axis[name]["total"] = np.full(values.shape[0], int(mask.sum()), dtype=int)
    counts["total"] = np.full(values.shape[0], values.shape[1], dtype=int)
    counts["by_axis"] = by_axis
    return counts


def zero_causes(tensor: np.ndarray, pairs: np.ndarray, values: np.ndarray) -> dict:
    a = tensor[:, pairs[:, 0, 0], pairs[:, 0, 1], pairs[:, 0, 2]]
    b = tensor[:, pairs[:, 1, 0], pairs[:, 1, 1], pairs[:, 1, 2]]
    null_a = np.all(a == 0, axis=(-2, -1)); null_b = np.all(b == 0, axis=(-2, -1))
    z = values == 0
    return {"BOTH_NULL": (z & null_a & null_b).sum(1),
            "NULL_ENDPOINT": (z & (null_a ^ null_b)).sum(1),
            "ORTHOGONAL_NONZERO": (z & ~null_a & ~null_b).sum(1)}


def axial_equivalence_error(tensor: np.ndarray, pairs: np.ndarray, values: np.ndarray) -> float:
    a = tensor[:, pairs[:, 0, 0], pairs[:, 0, 1], pairs[:, 0, 2]]
    b = tensor[:, pairs[:, 1, 0], pairs[:, 1, 1], pairs[:, 1, 2]]
    return float(np.max(np.abs(values - 2 * np.einsum("...i,...i->...", axial_dual(a), axial_dual(b)))))
