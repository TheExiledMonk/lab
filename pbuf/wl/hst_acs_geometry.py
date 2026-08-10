"""Continuous, radiometry-free ACS detector geometry primitives (Dev134)."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping
import numpy as np

CONTRACT_VERSION = "PBUF_HST_ACS_WFC_GEOMETRY_V1"
MAPPING_STATUSES = ("ACTIVE_CHIP", "INTER_CHIP_GAP", "OUTSIDE_DETECTOR", "BOUNDARY", "INVALID")


@dataclass(frozen=True)
class ACSDetectorGeometricEvent:
    event_uid: str
    exposure_uid: str
    instrument_frame_position: tuple[float, float, float]
    instrument_frame_direction: tuple[float, float, float]
    ideal_detector_u: float
    ideal_detector_v: float
    physical_detector_x: float
    physical_detector_y: float
    chip_id: str | None
    active_detector_status: str
    gap_status: str
    distortion_status: str
    reverse_metadata: Mapping[str, Any]
    dev133_provenance: Any
    dev132_provenance: Any
    bundle_provenance: Any
    latent_reverse_provenance: Any


@dataclass(frozen=True)
class AffineDetectorTransform:
    """Frozen source-derived ideal-to-physical continuous transform."""
    matrix: np.ndarray
    offset: np.ndarray
    exposure_uid: str = "SYNTHETIC"

    def __post_init__(self):
        m, b = np.asarray(self.matrix, np.float64), np.asarray(self.offset, np.float64)
        if m.shape != (2, 2) or b.shape != (2,) or not np.all(np.isfinite(m)) or not np.all(np.isfinite(b)):
            raise ValueError("ACS_TRANSFORM_INVALID")
        object.__setattr__(self, "matrix", m); object.__setattr__(self, "offset", b)

    def forward(self, uv: np.ndarray) -> np.ndarray:
        return np.asarray(uv, np.float64) @ self.matrix.T + self.offset

    def reverse(self, xy: np.ndarray) -> np.ndarray:
        if abs(np.linalg.det(self.matrix)) <= np.finfo(float).eps:
            raise ValueError("ACS_TRANSFORM_NOT_REVERSIBLE")
        return (np.asarray(xy, np.float64) - self.offset) @ np.linalg.inv(self.matrix).T

    @property
    def reverse_classification(self) -> str:
        return "EXACTLY_REVERSIBLE" if abs(abs(np.linalg.det(self.matrix)) - 1) <= 1e-14 else "NUMERICALLY_REVERSIBLE"


@dataclass(frozen=True)
class RectangularChip:
    chip_id: str
    x0: float
    x1: float
    y0: float
    y1: float


def classify_detector_points(xy: np.ndarray, chips: tuple[RectangularChip, ...], *, tolerance: float = 0.) -> tuple[np.ndarray, np.ndarray]:
    """Classify continuous coordinates with lower-inclusive/upper-exclusive interiors."""
    p = np.asarray(xy, np.float64); status = np.full(len(p), "OUTSIDE_DETECTOR", dtype="U20")
    chip_id = np.full(len(p), "", dtype="U16"); finite = np.all(np.isfinite(p), axis=1); status[~finite] = "INVALID"
    for c in chips:
        boundary = finite & (((np.abs(p[:,0]-c.x0)<=tolerance)|(np.abs(p[:,0]-c.x1)<=tolerance)) & (p[:,1]>=c.y0)&(p[:,1]<=c.y1) |
                            ((np.abs(p[:,1]-c.y0)<=tolerance)|(np.abs(p[:,1]-c.y1)<=tolerance)) & (p[:,0]>=c.x0)&(p[:,0]<=c.x1))
        inside = finite & (p[:,0]>=c.x0) & (p[:,0]<c.x1) & (p[:,1]>=c.y0) & (p[:,1]<c.y1) & ~boundary
        status[inside] = "ACTIVE_CHIP"; chip_id[inside] = c.chip_id
        status[boundary] = "BOUNDARY"; chip_id[boundary] = c.chip_id
    if len(chips) == 2:
        a, b = sorted(chips, key=lambda c: c.y0)
        gap = finite & (p[:,0] >= max(a.x0,b.x0)) & (p[:,0] < min(a.x1,b.x1)) & (p[:,1] >= a.y1) & (p[:,1] < b.y0)
        status[gap] = "INTER_CHIP_GAP"
    return status, chip_id


def basis_audit(ex: np.ndarray, ey: np.ndarray, normal: np.ndarray, tol: float = 1e-12) -> dict[str, Any]:
    b = np.stack((ex, ey, normal)).astype(float); gram = b @ b.T
    return {"finite": bool(np.all(np.isfinite(b))), "orthonormal": bool(np.allclose(gram, np.eye(3), rtol=0, atol=tol)),
            "handedness": float(np.dot(np.cross(b[0], b[1]), b[2])), "valid": bool(np.all(np.isfinite(b)) and np.allclose(gram,np.eye(3),rtol=0,atol=tol))}

