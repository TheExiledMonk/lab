"""Information-preserving arrival-to-instrument contract (Dev131).

No response law, rasterization, target data, or detector pixels live here.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
import hashlib
import json
from typing import Any, Mapping

import numpy as np

CONTRACT_VERSION = "PBUF_INSTRUMENT_INTERACTION_CONTRACT_V1"
REVERSIBILITY_CLASSES = ("FORWARD_REQUIRED", "REVERSE_REQUIRED", "BIDIRECTIONAL_REQUIRED",
                         "DERIVED_RECOMPUTABLE", "ARCHIVAL_LATENT", "PROVEN_REDUNDANT")
ELIGIBILITY_CLASSES = ("PHYSICAL_INSTRUMENT_INPUT_CANDIDATE", "LATENT_TRANSPORT_STATE",
                       "REVERSE_ONLY_PROVENANCE", "DERIVED_DIAGNOSTIC", "NOT_YET_CLASSIFIED")
AVAILABILITY_CLASSES = ("AVAILABLE", "PARTIALLY_AVAILABLE", "NOT_TRACKED",
                        "NOT_DEFINED_IN_CURRENT_MODEL")


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def canonical_sha256(value: Any) -> str:
    return hashlib.sha256(canonical_json(value).encode()).hexdigest()


def event_uid(ray_index: int, receiver_row_index: int, arrival_index: int) -> str:
    """Stable immutable-identity digest; deliberately not a random UUID."""
    raw = f"PBUF-EVENT-V1:{int(ray_index)}:{int(receiver_row_index)}:{int(arrival_index)}"
    return hashlib.sha256(raw.encode()).hexdigest()


def event_uids(ray_index, receiver_row_index, arrival_index=None) -> np.ndarray:
    ray = np.asarray(ray_index, dtype=np.int64)
    rec = np.asarray(receiver_row_index, dtype=np.int64)
    arr = np.arange(len(ray), dtype=np.int64) if arrival_index is None else np.asarray(arrival_index, dtype=np.int64)
    if not (ray.shape == rec.shape == arr.shape):
        raise ValueError("identity arrays must have equal shapes")
    return np.asarray([event_uid(a, b, c) for a, b, c in zip(ray, rec, arr)], dtype="U64")


@dataclass(frozen=True)
class InstrumentInteractionEvent:
    identity: Mapping[str, Any]
    surface_position: Mapping[str, Any]
    arrival_direction: Mapping[str, Any]
    incidence_state: Mapping[str, Any]
    transport_provenance: Mapping[str, Any]
    receiver_latent_state: Mapping[str, Any]
    bundle_provenance: Mapping[str, Any]
    interaction_eligibility: Mapping[str, str]
    reversibility_metadata: Mapping[str, Any]


@dataclass(frozen=True)
class InstrumentInteractionContract:
    version: str = CONTRACT_VERSION
    event_input_schema: Mapping[str, Any] = field(default_factory=dict)
    provenance_schema: Mapping[str, Any] = field(default_factory=dict)
    latent_state_schema: Mapping[str, Any] = field(default_factory=dict)
    neighbor_bundle_schema: Mapping[str, Any] = field(default_factory=dict)
    reverse_index_schema: Mapping[str, Any] = field(default_factory=dict)
    availability_flags: Mapping[str, Any] = field(default_factory=dict)
    cardinality: str = "many-to-many capable; current multiplicity audited separately"
    target_access: bool = False
    hst_pixel_access: bool = False

    def serialize(self) -> dict:
        return asdict(self)

    @property
    def sha256(self) -> str:
        return canonical_sha256(self.serialize())


def input_availability(available_fields) -> list[dict]:
    fields = set(available_fields)
    definitions = [
        ("I0", "surface location", ("arrival_u", "arrival_v"), "AVAILABLE"),
        ("I1", "arrival direction", ("arrival_dir_u", "arrival_dir_v", "arrival_dir_n"), "AVAILABLE"),
        ("I2", "incidence geometry", ("receiver_incidence_cosine", "receiver_incidence_angle"), "AVAILABLE"),
        ("I3", "event weight / carried signal", (), "NOT_DEFINED_IN_CURRENT_MODEL"),
        ("I4", "wavelength/frequency state", (), "NOT_DEFINED_IN_CURRENT_MODEL"),
        ("I5", "timing/path state", ("path_length", "path_excess"), "PARTIALLY_AVAILABLE"),
        ("I6", "polarization state", (), "NOT_DEFINED_IN_CURRENT_MODEL"),
        ("I7", "phase/coherence state", (), "NOT_DEFINED_IN_CURRENT_MODEL"),
        ("I8", "local bundle geometry", ("final_J11", "final_J12", "final_J21", "final_J22"), "AVAILABLE"),
        ("I9", "transport provenance", ("ray_index", "receiver_row_index", "launch_grid_index"), "AVAILABLE"),
    ]
    out = []
    for code, category, required, default in definitions:
        present = [x for x in required if x in fields]
        status = default
        if required and not present: status = "NOT_TRACKED"
        elif required and len(present) < len(required): status = "PARTIALLY_AVAILABLE"
        out.append({"category": code, "name": category, "availability": status,
                    "tracked_fields": present, "missing_fields": [x for x in required if x not in fields]})
    return out


def default_contract(availability=None) -> InstrumentInteractionContract:
    required = lambda classification, fields: {"classification": classification, "fields": list(fields)}
    return InstrumentInteractionContract(
        event_input_schema={
            "surface_position": required("BIDIRECTIONAL_REQUIRED", ("arrival_u", "arrival_v")),
            "arrival_direction": required("BIDIRECTIONAL_REQUIRED", ("arrival_dir_u", "arrival_dir_v", "arrival_dir_n")),
            "incidence_state": required("FORWARD_REQUIRED", ("receiver_incidence_cosine", "receiver_incidence_angle", "intersection_t")),
        },
        provenance_schema={
            "identity": required("BIDIRECTIONAL_REQUIRED", ("event_uid", "ray_index", "receiver_row_index", "arrival_index")),
            "launch": required("REVERSE_REQUIRED", ("launch_u", "launch_v", "launch_grid_index")),
            "trajectory": required("ARCHIVAL_LATENT", ("immutable Dev128/Dev129 references",)),
        },
        latent_state_schema={
            "path": required("ARCHIVAL_LATENT", ("path_length", "path_excess", "trajectory history reference")),
            "receiver_depth": required("REVERSE_REQUIRED", ("global_receive_position",)),
        },
        neighbor_bundle_schema={
            "launch_neighbors": {"scales": [1, 2, 4, 8, 16, 32], "storage": "regenerable from frozen launch_grid_index"},
            "arrival_neighbors": {"k": [4, 8, 16, 32, 64], "tie_break": "canonical arrival index"},
            "bundle_transport": {"classification": "BIDIRECTIONAL_REQUIRED", "storage": "immutable Dev129 reference"},
        },
        reverse_index_schema={"event_to_launch": "int64[N]", "launch_to_events": "CSR offsets + event indices",
                              "bundle_to_events": "regenerable membership", "receiver_to_event": "CSR capable"},
        availability_flags={x["category"]: x["availability"] for x in (availability or [])},
    )
