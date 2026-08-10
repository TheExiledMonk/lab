"""Neutral, source-supplied scalar state carried beside a frozen PBUF ray.

This module deliberately assigns no energy, momentum, frequency, wavelength,
phase, or amplitude meaning to ``q_scalar``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np


Q0_SUITE = (0.25, 0.5, 1.0, 2.0, 4.0, 8.0)


@dataclass(frozen=True)
class ZeroMassScalarSource:
    q_scalar: float
    source_uid: str = "source-0"

    def __post_init__(self) -> None:
        if not np.isfinite(self.q_scalar):
            raise ValueError("q_scalar must be finite")


@dataclass(frozen=True)
class ScalarStep:
    path_step: int
    path_position: float
    q_before: float
    q_after: float
    R_step: float
    delta_log_q: float
    local_driver: float | None
    medium_state_before: dict[str, Any] = field(default_factory=dict)
    medium_state_after: dict[str, Any] = field(default_factory=dict)


@dataclass
class ZeroMassScalarState:
    q_emit: float
    q_scalar: float | None = None
    candidate_id: str = "U01"
    history: list[ScalarStep] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.q_emit = float(self.q_emit)
        self.q_scalar = self.q_emit if self.q_scalar is None else float(self.q_scalar)
        if not np.isfinite(self.q_emit) or not np.isfinite(self.q_scalar):
            raise ValueError("scalar state must be finite")

    @property
    def q_receive(self) -> float:
        return float(self.q_scalar)

    @property
    def q_ratio(self) -> float:
        if self.q_emit == 0.0:
            raise ZeroDivisionError("relative transport is undefined for q_emit=0")
        return self.q_receive / self.q_emit

    def received_diagnostic(self, history_reference: str | None = None) -> dict[str, Any]:
        return {"q_emit": self.q_emit, "q_receive": self.q_receive,
                "q_ratio": self.q_ratio, "q_history_reference": history_reference,
                "q_transport_candidate_id": self.candidate_id}


def source_scalar_ontology_contract() -> dict[str, Any]:
    return {"contract": "PBUF_ZERO_MASS_SOURCE_SCALAR_V1",
            "source_scalar_required": True,
            "source_scalar_generated_by_medium": False,
            "source_scalar_initial_condition_allowed": True,
            "absolute_normalization_required": False,
            "relative_transport_primary": True,
            "identified_physical_quantity": None,
            "identification_status": "UNRESOLVED"}


def scalar_state_schema() -> dict[str, Any]:
    return {"state": ["position", "direction", "q_scalar"],
            "history_index": ["path_step", "native_path_position", "propagation_progression"],
            "forbidden_time_index": True,
            "q_semantics": "NEUTRAL"}
