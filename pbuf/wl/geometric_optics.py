"""Target-blind geometric optical transport (Dev133).

Only positions, unit directions, and provenance are handled here.  Interaction
records are immutable and append-only; no radiometric or wave state is created.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
import hashlib
from typing import Any, Mapping, Sequence

import numpy as np

from pbuf.wl.arrival_formation import BASIS_TOL, EPS_T, ReceiverPlane, intersect_rays
from pbuf.wl.optical_interaction_state import canonical_sha256

CONTRACT_VERSION = "PBUF_GEOMETRIC_OPTICAL_INTERACTION_V1"
INTERACTIONS = ("FREE_PROPAGATION", "APERTURE_TEST", "IDEAL_DIRECTION_TRANSFORM", "PLANE_INTERSECTION")
REVERSE_CLASSES = ("EXACTLY_REVERSIBLE", "REVERSIBLE_WITH_METADATA", "LOCALLY_REVERSIBLE",
                   "MULTIVALUED_INVERSE", "INFORMATION_LOSSY", "NONREVERSIBLE", "UNKNOWN")
APERTURE_STATUSES = ("PASS_APERTURE", "BLOCKED_BY_APERTURE", "ON_APERTURE_BOUNDARY",
                     "INVALID_APERTURE_INTERSECTION")


def _v3(value: Any, name: str) -> np.ndarray:
    value = np.asarray(value, dtype=np.float64)
    if value.shape != (3,) or not np.all(np.isfinite(value)):
        raise ValueError(f"{name} must be a finite 3-vector")
    return value


@dataclass(frozen=True)
class OpticalSurface:
    surface_id: str
    surface_type: str
    origin: np.ndarray
    basis_u: np.ndarray
    basis_v: np.ndarray
    normal: np.ndarray
    aperture_type: str = "NONE"
    aperture_parameters: Mapping[str, float] = field(default_factory=dict)
    interaction_type: str = "PLANE_INTERSECTION"
    interaction_parameters: Mapping[str, Any] = field(default_factory=dict)
    reverse_classification: str = "EXACTLY_REVERSIBLE"
    normal_convention: str = "basis_u_cross_basis_v_equals_normal"

    def __post_init__(self):
        for name in ("origin", "basis_u", "basis_v", "normal"):
            object.__setattr__(self, name, _v3(getattr(self, name), name))
        basis = np.stack((self.basis_u, self.basis_v, self.normal))
        if not np.allclose(basis @ basis.T, np.eye(3), rtol=0, atol=BASIS_TOL):
            raise ValueError("OPTICAL_SURFACE_BASIS_INVALID")
        if not np.isclose(np.dot(np.cross(self.basis_u, self.basis_v), self.normal), 1., rtol=0, atol=BASIS_TOL):
            raise ValueError("OPTICAL_SURFACE_HANDEDNESS_INVALID")
        if self.interaction_type not in INTERACTIONS:
            raise ValueError("unsupported optical interaction")
        if self.reverse_classification not in REVERSE_CLASSES:
            raise ValueError("invalid reverse classification")
        if self.aperture_type not in ("NONE", "CIRCULAR"):
            raise ValueError("unsupported aperture")
        if self.aperture_type == "CIRCULAR" and not float(self.aperture_parameters.get("radius", 0)) > 0:
            raise ValueError("circular aperture requires positive radius")
        if self.interaction_type == "IDEAL_DIRECTION_TRANSFORM" and not float(self.interaction_parameters.get("focal_distance", 0)) > 0:
            raise ValueError("ideal transform requires positive focal_distance")

    def receiver_plane(self) -> ReceiverPlane:
        return ReceiverPlane(self.origin, self.basis_u, self.basis_v, self.normal)

    def manifest(self) -> dict[str, Any]:
        d = asdict(self)
        for name in ("origin", "basis_u", "basis_v", "normal"):
            d[name] = getattr(self, name).tolist()
        d["aperture_parameters"] = dict(self.aperture_parameters)
        d["interaction_parameters"] = dict(self.interaction_parameters)
        return d


@dataclass(frozen=True)
class GeometricOpticalEvent:
    event_uid: str
    optical_record_uid: str
    surface_id: str
    interaction_index: int
    input_position: tuple[float, float, float]
    input_direction: tuple[float, float, float]
    intersection_position: tuple[float, float, float] | None
    surface_coordinates: tuple[float, float] | None
    incoming_direction: tuple[float, float, float]
    outgoing_direction: tuple[float, float, float] | None
    interaction_status: str
    aperture_status: str
    parent_event_identity: str
    reverse_metadata: Mapping[str, Any]
    upstream_provenance_reference: Any
    latent_state_reference: Any
    reverse_classification: str


def optical_record_uid(event_uid: Any, surface_id: str, interaction_index: int) -> str:
    payload = f"{event_uid}\0{surface_id}\0{int(interaction_index)}".encode()
    return hashlib.sha256(payload).hexdigest()


def aperture_test(u: np.ndarray, v: np.ndarray, radius: float, *, tolerance: float = EPS_T) -> np.ndarray:
    u, v = np.asarray(u, float), np.asarray(v, float)
    r = np.hypot(u, v)
    status = np.full(r.shape, "INVALID_APERTURE_INTERSECTION", dtype="U32")
    finite = np.isfinite(r)
    boundary = finite & np.isclose(r, radius, rtol=0, atol=tolerance * max(1., radius))
    status[finite & (r < radius) & ~boundary] = "PASS_APERTURE"
    status[finite & (r > radius) & ~boundary] = "BLOCKED_BY_APERTURE"
    status[boundary] = "ON_APERTURE_BOUNDARY"  # boundary is included in transmission
    return status


def ideal_direction_transform(points: np.ndarray, surface: OpticalSurface) -> np.ndarray:
    points = np.asarray(points, float)
    q = surface.origin + float(surface.interaction_parameters["focal_distance"]) * surface.normal
    delta = q - points
    norm = np.linalg.norm(delta, axis=1)
    if np.any(~np.isfinite(norm) | (norm <= 0)):
        raise ValueError("OPTICAL_DIRECTION_NORMALIZATION_FAILURE")
    return delta / norm[:, None]


def propagate_to_surface(position: np.ndarray, direction: np.ndarray, surface: OpticalSurface):
    """Reuse Dev130's frozen plane-intersection implementation and semantics."""
    p, t, d, denominator, status = intersect_rays(position, direction, surface.receiver_plane())
    delta = p - surface.origin
    u, v = delta @ surface.basis_u, delta @ surface.basis_v
    aperture = np.full(len(t), "PASS_APERTURE", dtype="U32")
    if surface.aperture_type == "CIRCULAR":
        aperture = aperture_test(u, v, float(surface.aperture_parameters["radius"]))
    invalid = ~np.isin(status, ("FORWARD_INTERSECTION", "ON_SURFACE"))
    aperture[invalid] = "INVALID_APERTURE_INTERSECTION"
    outgoing = d.copy()
    if surface.interaction_type == "IDEAL_DIRECTION_TRANSFORM":
        usable = ~invalid & np.isin(aperture, ("PASS_APERTURE", "ON_APERTURE_BOUNDARY"))
        outgoing[:] = np.nan
        outgoing[usable] = ideal_direction_transform(p[usable], surface)
    transmitted = ~invalid & np.isin(aperture, ("PASS_APERTURE", "ON_APERTURE_BOUNDARY"))
    outgoing[~transmitted] = np.nan
    norms = np.linalg.norm(outgoing[transmitted], axis=1)
    if np.any(np.abs(norms - 1.) > 512 * np.finfo(float).eps):
        raise ValueError("OPTICAL_DIRECTION_NORMALIZATION_FAILURE")
    return {"intersection_position": p, "intersection_t": t, "incoming_direction": d,
            "outgoing_direction": outgoing, "surface_u": u, "surface_v": v,
            "intersection_status": status, "aperture_status": aperture,
            "transmitted": transmitted, "denominator": denominator}


def make_interaction_records(event_uids: Sequence[Any], position: np.ndarray, direction: np.ndarray,
                             surface: OpticalSurface, interaction_index: int, *,
                             parent_event_uids: Sequence[Any] | None = None,
                             upstream_references: Sequence[Any] | None = None,
                             latent_references: Sequence[Any] | None = None):
    result = propagate_to_surface(position, direction, surface)
    n = len(result["intersection_t"])
    parent_event_uids = event_uids if parent_event_uids is None else parent_event_uids
    upstream_references = event_uids if upstream_references is None else upstream_references
    latent_references = event_uids if latent_references is None else latent_references
    records = []
    for i in range(n):
        point = result["intersection_position"][i]
        out = result["outgoing_direction"][i]
        blocked = result["aperture_status"][i] == "BLOCKED_BY_APERTURE"
        records.append(GeometricOpticalEvent(
            str(event_uids[i]), optical_record_uid(event_uids[i], surface.surface_id, interaction_index),
            surface.surface_id, int(interaction_index), tuple(map(float, position[i])), tuple(map(float, direction[i])),
            tuple(map(float, point)) if np.all(np.isfinite(point)) else None,
            (float(result["surface_u"][i]), float(result["surface_v"][i])) if np.all(np.isfinite(point)) else None,
            tuple(map(float, result["incoming_direction"][i])), tuple(map(float, out)) if np.all(np.isfinite(out)) else None,
            str(result["intersection_status"][i]), str(result["aperture_status"][i]), str(parent_event_uids[i]),
            {"intersection_t": float(result["intersection_t"][i]), "stored_incoming_direction": tuple(map(float, result["incoming_direction"][i])),
             "focal_distance": surface.interaction_parameters.get("focal_distance"), "upstream_state_retained": True},
            upstream_references[i], latent_references[i], "INFORMATION_LOSSY" if blocked else surface.reverse_classification))
    return tuple(records), result


def system_manifest(surfaces: Sequence[OpticalSurface]) -> dict[str, Any]:
    ids = [s.surface_id for s in surfaces]
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate surface_id")
    return {"contract": CONTRACT_VERSION, "surface_order": ids, "surfaces": [s.manifest() for s in surfaces],
            "target_access": False, "hst_pixel_access": False, "physical_intensity_formation": False,
            "spectral_optics": False, "phase_optics": False, "diffraction": False}


def system_sha256(surfaces: Sequence[OpticalSurface]) -> str:
    return canonical_sha256(system_manifest(surfaces))

