"""Provenance-strict PBUF receiver to tangent-plane bridge (Dev135)."""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Mapping
import hashlib
import json
import numpy as np

BRIDGE_VERSION = "PBUF_RECEIVER_TO_HST_COORDINATE_BRIDGE_V1"
SCALE_CLASSES = ("EXPLICIT_PHYSICAL_SCALE", "EXPLICIT_ANGULAR_SCALE",
                 "DERIVABLE_PHYSICAL_SCALE", "DERIVABLE_ANGULAR_SCALE")


def canonical_sha256(value: Any) -> str:
    blob = json.dumps(value, sort_keys=True, separators=(",", ":"),
                      allow_nan=False).encode()
    return hashlib.sha256(blob).hexdigest()


@dataclass(frozen=True)
class PBUFReceiverSkyBridge:
    bridge_version: str
    a0_origin_definition: str
    a0_basis_u: tuple[float, float, float]
    a0_basis_v: tuple[float, float, float]
    a0_normal: tuple[float, float, float]
    scale_type: str
    scale_value: float
    scale_units: str
    scale_provenance: Mapping[str, Any]
    sky_origin: tuple[float, float]
    sky_basis: str
    orientation_transform: tuple[tuple[float, float], tuple[float, float]]
    reflection_status: bool
    gauge_status: Mapping[str, str]
    derivation_chain: tuple[str, ...]
    source_manifest_sha256: str
    reverse_classification: str = "NUMERICALLY_REVERSIBLE"
    scale_uncertainty: float | None = None
    uncertainty_source: str | None = None

    def __post_init__(self) -> None:
        if self.bridge_version != BRIDGE_VERSION:
            raise ValueError("PBUF_SKY_BRIDGE_VERSION_INVALID")
        if self.scale_type not in SCALE_CLASSES or not np.isfinite(self.scale_value) or self.scale_value <= 0:
            raise ValueError("PBUF_RECEIVER_PHYSICAL_SCALE_NOT_ESTABLISHED")
        r = np.asarray(self.orientation_transform, dtype=np.float64)
        if r.shape != (2, 2) or not np.allclose(r.T @ r, np.eye(2), atol=1e-12, rtol=0):
            raise ValueError("PBUF_RECEIVER_SKY_ATTITUDE_NOT_ESTABLISHED")
        if bool(np.linalg.det(r) < 0) != self.reflection_status:
            raise ValueError("A0_SKY_REFLECTION_METADATA_MISMATCH")

    def forward(self, a0: np.ndarray) -> np.ndarray:
        p = np.asarray(a0, dtype=np.float64)
        return np.asarray(self.sky_origin) + self.scale_value * (p @ np.asarray(self.orientation_transform).T)

    def reverse(self, sky: np.ndarray) -> np.ndarray:
        p = np.asarray(sky, dtype=np.float64) - np.asarray(self.sky_origin)
        return (p / self.scale_value) @ np.asarray(self.orientation_transform)

    def manifest(self) -> dict[str, Any]:
        return asdict(self)

    @property
    def transform_sha256(self) -> str:
        return canonical_sha256(self.manifest())


@dataclass(frozen=True)
class HSTReferenceGeometricEvent:
    event_uid: str
    a0_u: float
    a0_v: float
    a0_direction: tuple[float, float, float]
    hst_ref_u: float
    hst_ref_v: float
    hst_ref_direction: tuple[float, float, float]
    scale_provenance: Mapping[str, Any]
    origin_provenance: Mapping[str, Any]
    attitude_provenance: Mapping[str, Any]
    bridge_status: str
    reverse_metadata: Mapping[str, Any]
    full_upstream_provenance_reference: Mapping[str, Any]


def roundtrip_audit(bridge: PBUFReceiverSkyBridge, points: np.ndarray) -> dict[str, Any]:
    p = np.asarray(points, dtype=np.float64)
    err = np.linalg.norm(bridge.reverse(bridge.forward(p)) - p, axis=-1)
    return {"max_error": float(np.max(err)), "rms": float(np.sqrt(np.mean(err**2))),
            "p95": float(np.percentile(err, 95)), "p99": float(np.percentile(err, 99)),
            "orientation_error": float(np.max(np.abs(np.asarray(bridge.orientation_transform).T @
                                                       np.asarray(bridge.orientation_transform)-np.eye(2)))),
            "identity_mismatches": int(np.count_nonzero(err > 1e-12))}
