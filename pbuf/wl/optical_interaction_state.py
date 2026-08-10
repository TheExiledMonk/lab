"""Physical optical interaction state contract (Dev132).

This module classifies and preserves event state.  It deliberately performs no
optical response, rasterization, detector formation, or source reconstruction.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from enum import Enum
import hashlib
import json
from typing import Any, Mapping

import numpy as np

CONTRACT_VERSION = "PBUF_OPTICAL_INTERACTION_STATE_V1"


class Availability(str, Enum):
    DIRECTLY_AVAILABLE = "DIRECTLY_AVAILABLE"
    DERIVABLE_FROM_FROZEN_STATE = "DERIVABLE_FROM_FROZEN_STATE"
    PARTIALLY_DERIVABLE = "PARTIALLY_DERIVABLE"
    NOT_TRACKED = "NOT_TRACKED"
    NOT_DEFINED_IN_CURRENT_MODEL = "NOT_DEFINED_IN_CURRENT_MODEL"
    REQUIRES_NEW_PHYSICS_OR_SOURCE_INPUT = "REQUIRES_NEW_PHYSICS_OR_SOURCE_INPUT"


RELEVANCE = ("REQUIRED", "POSSIBLY_RELEVANT", "NOT_CURRENTLY_REQUIRED", "UNKNOWN")
DERIVATION_CLASSES = ("EXACT_ALGEBRAIC", "EXACT_GEOMETRIC", "NUMERICAL_SUMMARY",
                      "REQUIRES_ASSUMPTION", "NOT_DERIVABLE")
FLAG_NAMES = ("has_physical_weight", "has_spectral_state", "has_arrival_time",
              "has_phase", "has_polarization", "has_bundle_state",
              "has_reverse_provenance")


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def canonical_sha256(value: Any) -> str:
    return hashlib.sha256(canonical_json(value).encode()).hexdigest()


@dataclass(frozen=True)
class OpticalInteractionState:
    geometry: Mapping[str, Any] = field(default_factory=dict)
    carried_signal: Mapping[str, Any] = field(default_factory=dict)
    spectral: Mapping[str, Any] = field(default_factory=dict)
    temporal: Mapping[str, Any] = field(default_factory=dict)
    phase_coherence: Mapping[str, Any] = field(default_factory=dict)
    polarization: Mapping[str, Any] = field(default_factory=dict)
    bundle_state: Mapping[str, Any] = field(default_factory=dict)
    latent_reverse_state: Mapping[str, Any] = field(default_factory=dict)
    availability_metadata: Mapping[str, bool] = field(default_factory=dict)
    provenance: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        missing = set(FLAG_NAMES) - set(self.availability_metadata)
        if missing:
            raise ValueError(f"missing availability flags: {sorted(missing)}")

    def serialize(self) -> dict[str, Any]:
        return asdict(self)


def state_from_event(event: Mapping[str, Any]) -> OpticalInteractionState:
    """Losslessly classify a synthetic or production event without defaults."""
    present = lambda *names: any(n in event and event[n] is not None for n in names)
    flags = {
        "has_physical_weight": present("physical_weight", "photon_count", "packet_energy"),
        "has_spectral_state": present("wavelength", "frequency", "photon_energy", "spectral_bin"),
        "has_arrival_time": present("arrival_time"),
        "has_phase": present("phase", "complex_amplitude", "wavefront_phase"),
        "has_polarization": present("stokes", "jones", "polarization_vector", "helicity"),
        "has_bundle_state": present("final_J11", "bundle_jacobian", "bundle_area_ratio"),
        "has_reverse_provenance": present("event_uid") and present("ray_index") and present("receiver_row_index") and present("launch_grid_index"),
    }
    pick = lambda names: {k: event[k] for k in names if k in event and event[k] is not None}
    return OpticalInteractionState(
        geometry=pick(("arrival_u","arrival_v","arrival_dir_u","arrival_dir_v","arrival_dir_n","incidence_cosine","incidence_angle","intersection_t")),
        carried_signal=pick(("physical_weight","photon_count","packet_energy","numerical_sampling_weight")),
        spectral=pick(("wavelength","frequency","photon_energy","spectral_bin","filter_band_identity")),
        temporal=pick(("path_length","path_excess","number_of_steps","arrival_time","relative_delay","emission_time")),
        phase_coherence=pick(("phase","phase_offset","complex_amplitude","coherence_length","wavefront_phase")),
        polarization=pick(("stokes","jones","polarization_vector","helicity")),
        bundle_state={k:v for k,v in event.items() if k.startswith(("final_J","final_H","mean_J","rms_J","minimum_area","maximum_area")) or k in ("bundle_jacobian","bundle_area_ratio","bundle_anisotropy")},
        latent_reverse_state=pick(("path_history_reference","receiver_depth","neighbor_provenance","conditioning_state","reverse_candidate_metadata")),
        availability_metadata=flags,
        provenance=pick(("event_uid","arrival_index","ray_index","receiver_row_index","launch_grid_index","launch_u","launch_v")),
    )


def availability_masks(events: list[Mapping[str, Any]]) -> dict[str, np.ndarray]:
    states = [state_from_event(e) for e in events]
    return {name: np.asarray([s.availability_metadata[name] for s in states], dtype=bool)
            for name in FLAG_NAMES}


def validate_derivation_graph(graph: Mapping[str, Any]) -> None:
    edges = graph.get("edges", [])
    adjacency: dict[str, list[str]] = {}
    for edge in edges:
        adjacency.setdefault(edge["source"], []).append(edge["derived"])
    visiting: set[str] = set(); done: set[str] = set()
    def visit(node: str):
        if node in visiting: raise ValueError("OPTICAL_STATE_DERIVATION_CYCLE")
        if node in done: return
        visiting.add(node)
        for nxt in adjacency.get(node, []): visit(nxt)
        visiting.remove(node); done.add(node)
    for node in set(adjacency): visit(node)


def contract_schema() -> dict[str, Any]:
    return {
        "version": CONTRACT_VERSION,
        "groups": ["geometry","carried_signal","spectral","temporal","phase_coherence",
                   "polarization","bundle_state","latent_reverse_state","availability_metadata","provenance"],
        "availability_flags": list(FLAG_NAMES),
        "missing_value_semantics": "explicit event-wise boolean mask; numeric NaN is not semantic",
        "primitive_first": True,
        "target_access": False,
        "hst_pixel_access": False,
        "detector_pixels_generated": False,
    }
