"""Passive reconstruction of the frozen G3D direction normalization (Dev147)."""
from __future__ import annotations

import numpy as np


def decompose_update(previous_direction, response, path_step):
    """Return the complete pre-normalization state without changing propagation."""
    n = np.asarray(previous_direction, dtype=np.float64)
    r = np.asarray(response, dtype=np.float64)
    if n.shape[-1] != 3 or r.shape != n.shape:
        raise ValueError("direction and response must have matching (..., 3) shapes")
    delta = float(path_step) * r
    raw = n + delta
    magnitude = np.linalg.norm(raw, axis=-1)
    if np.any(magnitude <= 0) or not np.all(np.isfinite(magnitude)):
        raise ValueError("raw direction must have finite positive magnitude")
    normalized = raw / magnitude[..., None]
    longitudinal_scalar = np.sum(delta * n, axis=-1)
    parallel = longitudinal_scalar[..., None] * n
    transverse = delta - parallel
    return {
        "raw_vector": raw,
        "raw_magnitude": magnitude,
        "normalized_vector": normalized,
        "normalization_factor": magnitude,
        "delta_vector": delta,
        "delta_magnitude": np.linalg.norm(delta, axis=-1),
        "longitudinal_scalar": longitudinal_scalar,
        "parallel_vector": parallel,
        "parallel_magnitude": np.linalg.norm(parallel, axis=-1),
        "transverse_vector": transverse,
        "transverse_magnitude": np.linalg.norm(transverse, axis=-1),
    }


def audit_history(initial_direction, responses, path_step):
    """Replay the exact direction recurrence and retain every normally lost value."""
    responses = np.asarray(responses, dtype=np.float64)
    n = np.asarray(initial_direction, dtype=np.float64)
    if responses.ndim < 2 or responses.shape[-1] != 3:
        raise ValueError("responses must have shape (steps, ..., 3)")
    rows = []
    for response in responses:
        row = decompose_update(n, response, path_step)
        rows.append(row)
        n = row["normalized_vector"]
    return {key: np.stack([row[key] for row in rows]) for key in rows[0]} if rows else {}


def normalization_classification():
    return {
        "operation": "v_raw -> v_raw / |v_raw|",
        "classification": "GEOMETRIC_DIRECTION_NORMALIZATION",
        "loss_class": "N02 removes arbitrary numerical scale",
        "physical_c_derivation": False,
        "reason": "the magnitude depends on path_step and response representation and is not a persistent state",
    }

