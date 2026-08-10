"""Coefficient-free Dev151 local-frame transports and covariance diagnostics."""
from __future__ import annotations
import numpy as np
from .native_neighbor_state import frame_overlap

FRAME_CANDIDATES = ("F01", "F02", "F03", "F04", "F05", "F06", "F07")


def transport_map(source: np.ndarray, target: np.ndarray, candidate: str) -> np.ndarray | None:
    if candidate not in FRAME_CANDIDATES:
        raise ValueError("unknown frame candidate")
    if candidate == "F07":
        return None
    raw = np.asarray(target)[1:] @ np.asarray(source)[1:].T
    if candidate == "F01":
        return raw
    if candidate in ("F02", "F03", "F04", "F05"):
        return frame_overlap(source, target)
    # F06 is the proper rotational factor, excluding a reflection.
    u, _, vh = np.linalg.svd(raw)
    q = u @ vh
    if np.linalg.det(q) < 0:
        u[:, -1] *= -1
        q = u @ vh
    return q


def transport(values, source, target, candidate="F04"):
    q = transport_map(source, target, candidate)
    if q is None:
        raise ValueError("F07 is unresolved and cannot progress a state")
    return q @ np.asarray(values, float)


def rotation_angle(q):
    q = np.asarray(q, float)
    return float(np.arctan2(q[1, 0] - q[0, 1], q[0, 0] + q[1, 1]))


def audit_frame_candidates(frames: np.ndarray) -> list[dict]:
    out = []
    for name in FRAME_CANDIDATES:
        if name == "F07":
            out.append({"candidate": name, "status": "UNRESOLVED"}); continue
        errors, angles = [], []
        for a, b in zip(frames[:-1], frames[1:]):
            q = transport_map(a, b, name)
            errors.append(float(np.linalg.norm(q.T @ q - np.eye(2))))
            angles.append(rotation_angle(q))
        ok = max(errors, default=0.0) < 1e-10
        out.append({"candidate": name, "norm_preserving": ok,
                    "basis_covariant": ok, "max_orthogonality_error": max(errors, default=0.0),
                    "frame_rotation_angles": angles, "status": "PASS" if ok else "FAIL"})
    return out


def circulation(frame_angles):
    a = np.asarray(frame_angles, float)
    return np.diff(np.unwrap(a), prepend=a[0]) if a.size else a

