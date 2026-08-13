"""DEV213 preparation algebra for existing DEV182 native packet increments.

This module introduces neither a packet law nor an evolution law.  It makes
the already-audited DEV196 operation explicit so its same-progression-step
order and provenance can be inspected.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from pbuf.excitation.native_vector_pair_dynamics import VectorPairState


@dataclass(frozen=True)
class NativePreparation:
    structure_id: str
    preparation_operation: str
    internal_state_id: str
    placement_operation: str
    source_dev: str
    amplitude_source: str
    orientation_source: str
    displacement_increment: np.ndarray
    momentum_increment: np.ndarray

    def provenance(self) -> dict:
        return {key: getattr(self, key) for key in (
            "structure_id", "preparation_operation", "internal_state_id",
            "placement_operation", "source_dev", "amplitude_source",
            "orientation_source",
        )}


def inject(state: VectorPairState, preparation: NativePreparation) -> VectorPairState:
    """The unchanged DEV196 valid-state operation I_deltaX(X)=X+deltaX."""
    return VectorPairState(
        state.displacement + preparation.displacement_increment,
        state.momentum + preparation.momentum_increment,
        state.progression_step,
    )


def reverse_internal_state(preparation: NativePreparation, state_id: str) -> NativePreparation:
    """Apply DEV212 R_p after geometry and placement have already been fixed."""
    return NativePreparation(
        **{**preparation.__dict__, "internal_state_id": state_id,
           "momentum_increment": -preparation.momentum_increment.copy()}
    )


def exact_support(displacement: np.ndarray, momentum: np.ndarray) -> np.ndarray:
    """Exact nonzero state support; deliberately no magnitude cutoff."""
    return np.any(np.asarray(displacement) != 0.0, axis=-1) | np.any(np.asarray(momentum) != 0.0, axis=-1)
