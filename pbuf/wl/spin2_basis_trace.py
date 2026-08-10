"""Target-blind coordinate-frame and spin-2 basis tracing utilities.

The convention is Q=[[q1,q2],[q2,-q1]] and a physical orthogonal map A acts
as Q' = A Q A.T.  No target arrays are accepted by this module.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Iterable

import numpy as np

TOL = 1e-12
ORIENTATIONS_DEG = (0.0, 22.5, 45.0, 67.5, 90.0, 112.5, 135.0, 157.5)


@dataclass(frozen=True)
class CoordinateFrame:
    name: str
    x_axis_definition: str
    y_axis_definition: str
    handedness: str
    origin: str
    axis_signs: tuple[str, str]
    axis_order: tuple[str, str]
    units: str
    source_file: str
    source_function: str
    status: str = "RESOLVED"

    def manifest(self) -> dict:
        row = asdict(self)
        row["axis_signs"] = list(self.axis_signs)
        row["axis_order"] = list(self.axis_order)
        return row


def rotation_matrix(degrees: float) -> np.ndarray:
    phi = math.radians(degrees)
    return np.array([[math.cos(phi), -math.sin(phi)],
                     [math.sin(phi), math.cos(phi)]], dtype=np.float64)


def tensor_from_components(q1: float, q2: float) -> np.ndarray:
    return np.array([[q1, q2], [q2, -q1]], dtype=np.float64)


def components_from_tensor(q: np.ndarray) -> np.ndarray:
    q = np.asarray(q, dtype=np.float64)
    return np.array([(q[0, 0] - q[1, 1]) / 2.0,
                     (q[0, 1] + q[1, 0]) / 2.0])


def spin2_matrix(a: np.ndarray) -> np.ndarray:
    """Derive the component action, including reflections, from A Q A.T."""
    a = np.asarray(a, dtype=np.float64)
    if a.shape != (2, 2) or not np.allclose(a.T @ a, np.eye(2), atol=TOL):
        raise ValueError("basis transform must be a 2x2 orthogonal matrix")
    return np.column_stack([
        components_from_tensor(a @ tensor_from_components(*basis) @ a.T)
        for basis in ((1.0, 0.0), (0.0, 1.0))
    ])


def apply_spin2(q1, q2, a: np.ndarray):
    s = spin2_matrix(a)
    return s[0, 0] * q1 + s[0, 1] * q2, s[1, 0] * q1 + s[1, 1] * q2


def synthetic_state(angle_deg: float) -> np.ndarray:
    angle = math.radians(2.0 * angle_deg)
    return np.array([math.cos(angle), math.sin(angle)], dtype=np.float64)


def classify_transform(a: np.ndarray) -> dict:
    a = np.asarray(a, float)
    det = float(np.linalg.det(a))
    reflection = det < 0
    proper_angle = None if reflection else math.degrees(math.atan2(a[1, 0], a[0, 0]))
    identity = np.allclose(a, np.eye(2), atol=TOL)
    swap = bool(abs(a[0, 1]) > 1 - TOL and abs(a[1, 0]) > 1 - TOL)
    if identity: label = "IDENTITY"
    elif reflection and swap: label = "AXIS_SWAP"
    elif reflection: label = "REFLECTION"
    elif swap: label = "PHYSICAL_ROTATION"
    else: label = "PHYSICAL_ROTATION"
    return {"classification": label, "determinant": det,
            "rotation_angle_if_proper": proper_angle,
            "reflection_present": reflection, "axis_swap": swap,
            "x_flip": bool(np.allclose(a, np.diag([-1, 1]), atol=TOL)),
            "y_flip": bool(np.allclose(a, np.diag([1, -1]), atol=TOL)),
            "equivalent_physical_rotation_deg": proper_angle,
            "spin2_phase_deg": None if proper_angle is None else 2 * proper_angle}


def transform_record(source: str, target: str, a: np.ndarray, reference: str) -> dict:
    a = np.asarray(a, float)
    return {"source_frame": source, "target_frame": target,
            "matrix_2x2": a.tolist(), "spin2_matrix_2x2": spin2_matrix(a).tolist(),
            **classify_transform(a), "source_code_reference": reference,
            "basis_round_trip_max_error": float(np.max(np.abs(np.linalg.inv(a) @ a - np.eye(2)))),
            "spin2_round_trip_max_error": float(np.max(np.abs(np.linalg.inv(spin2_matrix(a)) @ spin2_matrix(a) - np.eye(2))))}


def compose(transforms: Iterable[np.ndarray]) -> np.ndarray:
    result = np.eye(2)
    for transform in transforms:
        result = np.asarray(transform, float) @ result
    return result


def trace_orientations(frame_transforms: list[tuple[str, np.ndarray]]) -> list[dict]:
    rows = []
    for angle in ORIENTATIONS_DEG:
        q = synthetic_state(angle)
        stages = [{"frame": frame_transforms[0][0], "physical_angle_deg": angle,
                   "q1": float(q[0]), "q2": float(q[1]),
                   "tensor_matrix": tensor_from_components(*q).tolist()}]
        cumulative = np.eye(2)
        for name, a in frame_transforms[1:]:
            cumulative = np.asarray(a) @ cumulative
            qt = spin2_matrix(cumulative) @ q
            stages.append({"frame": name, "physical_angle_deg": angle,
                           "q1": float(qt[0]), "q2": float(qt[1]),
                           "tensor_matrix": tensor_from_components(*qt).tolist()})
        rows.append({"physical_angle_deg": angle, "stages": stages})
    return rows
