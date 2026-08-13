"""Coefficient-free spatial winding observer for an existing N6 state.

This module deliberately creates no field: it reads DEV167 positive bond
forces, normalizes nonzero vectors exactly, and contracts them with fixed
lattice-edge tangents on an already supplied closed loop.
"""
from __future__ import annotations

import numpy as np

from pbuf.excitation.native_vector_pair_dynamics import pair_forces


def square_yz_loop(shape: tuple[int, int, int], center: tuple[int, int, int], radius: int):
    """The +x-oriented yz square, traversed +y, +z, -y, -z.

    Edges are ``(tail node, positive-axis, sign)``.  The routine rejects a
    periodic wrap rather than silently treating numerical-torus topology as a
    local loop.
    """
    x, cy, cz = center
    ny, nz = shape[1:]
    if radius <= 1 or 2 * radius >= min(ny, nz):
        raise ValueError("radius must be a contractible multi-cell loop")
    edges = []
    y0, y1, z0, z1 = cy - radius, cy + radius, cz - radius, cz + radius
    for y in range(y0, y1): edges.append(((x, y, z0), 1, 1))
    for z in range(z0, z1): edges.append(((x, y1, z), 2, 1))
    for y in range(y1, y0, -1): edges.append(((x, y, z1), 1, -1))
    for z in range(z1, z0, -1): edges.append(((x, y0, z), 2, -1))
    return tuple(edges)


def circulation(displacement: np.ndarray, edges):
    """Return every normalized-direction contribution and the exact sum."""
    positive = pair_forces(displacement)
    directions, tangents, raw, contributions, raw_contributions = [], [], [], [], []
    for tail, axis, sign in edges:
        tail = tuple(int(v) for v in tail)
        if sign == 1:
            force = positive[tail + (axis,)]
        else:
            head = list(tail); head[axis] -= 1
            force = -positive[tuple(head) + (axis,)]
        magnitude = np.linalg.norm(force)
        tangent = np.zeros(3); tangent[axis] = sign
        if magnitude == 0.0:
            return {"defined": False, "directions": np.asarray(directions), "tangents": np.asarray(tangents),
                    "raw": np.asarray(raw), "contributions": np.asarray(contributions),
                    "raw_contributions": np.asarray(raw_contributions)}
        direction = force / magnitude
        directions.append(direction); tangents.append(tangent); raw.append(force)
        contributions.append(float(direction @ tangent)); raw_contributions.append(float(force @ tangent))
    return {"defined": True, "directions": np.asarray(directions), "tangents": np.asarray(tangents),
            "raw": np.asarray(raw), "contributions": np.asarray(contributions),
            "raw_contributions": np.asarray(raw_contributions),
            "circulation": float(np.sum(contributions)), "raw_circulation": float(np.sum(raw_contributions))}


def reflect_x(displacement: np.ndarray, center_x: int) -> np.ndarray:
    """Exact lattice reflection about an integer x plane, including polar vectors."""
    n = displacement.shape[0]
    indices = (2 * center_x - np.arange(n)) % n
    result = displacement[indices].copy()
    result[..., 0] *= -1.0
    return result
