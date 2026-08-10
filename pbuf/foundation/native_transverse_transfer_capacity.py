"""Composition and accounting for loaded-link response and Dev152 transport."""
from __future__ import annotations
import numpy as np


def orthogonal_transport(angle):
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, -s], [s, c]], dtype=float)


def compose(response, rotation, state, ordering="capacity_then_rotation"):
    C, R, x = np.asarray(response), np.asarray(rotation), np.asarray(state)
    if ordering == "capacity_then_rotation": return R @ (C @ x)
    if ordering == "rotation_then_capacity": return C @ (R @ x)
    raise ValueError(ordering)


def norm_accounting(before, after):
    a, b = np.asarray(before), np.asarray(after)
    n0, n1 = float(np.sum(a*a)), float(np.sum(b*b))
    return {"before": n0, "after": n1, "drift": n1-n0, "conserved": bool(np.isclose(n0, n1, atol=1e-12))}


def ordering_audit(response, rotation):
    commutator = np.asarray(response) @ np.asarray(rotation) - np.asarray(rotation) @ np.asarray(response)
    return {"commutator_norm": float(np.linalg.norm(commutator)),
            "classification": "COMMUTING_EQUIVALENT" if np.linalg.norm(commutator) < 1e-12 else "ORDERING_PHYSICALLY_DISTINCT"}
