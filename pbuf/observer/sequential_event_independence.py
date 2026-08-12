"""Pure diagnostics shared by the DEV196 sequential-event audit."""
from __future__ import annotations

import numpy as np


def inject(state, displacement_packet: np.ndarray, momentum_packet: np.ndarray):
    """Apply the DEV182 additive packet operation to a valid native state."""
    from pbuf.excitation.native_vector_pair_dynamics import VectorPairState
    return VectorPairState(state.displacement + displacement_packet,
                           state.momentum + momentum_packet,
                           state.progression_step)


def support_mask(displacement: np.ndarray, momentum: np.ndarray) -> np.ndarray:
    """Exact native occupancy, matching the DEV195 serialized convention."""
    return np.any(np.abs(displacement) > 0, axis=-1) | np.any(np.abs(momentum) > 0, axis=-1)


def support_relation(a: np.ndarray, b: np.ndarray, recurrence_possible: bool) -> str:
    if recurrence_possible:
        return 'PERIODIC_RECURRENCE_CONTAMINATION'
    n = int(np.count_nonzero(a & b))
    if n == 0:
        return 'DISJOINT'
    if n == min(int(a.sum()), int(b.sum())):
        return 'DIRECT_OVERLAP'
    return 'PARTIAL_OVERLAP'


def component_summary(a: np.ndarray, overlap: np.ndarray, recurrence: np.ndarray) -> dict:
    a = np.asarray(a)
    mag = np.abs(a)
    flat = int(np.argmax(mag))
    index = tuple(int(x) for x in np.unravel_index(flat, a.shape))
    # State arrays are (time,x,y,z,component); scalar flux uses no component.
    time = index[0]
    location = index[1:4] if a.ndim == 5 else None
    return {
        'max_absolute_residual': float(mag.flat[flat]),
        'l2_diagnostic': float(np.linalg.norm(a)),
        'maximum_index': index,
        'location': location,
        'time_after_second_launch': time,
        'support_overlap_at_maximum': bool(overlap[(time,) + location]) if location is not None else None,
        'boundary_recurrence_possible_at_maximum': bool(recurrence[time]),
    }
