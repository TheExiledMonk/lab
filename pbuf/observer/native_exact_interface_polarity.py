"""DEV218 exact force helpers over the persisted DEV217 observer."""
from __future__ import annotations

import numpy as np


ATOL = 1e-12


def action_reaction_class(residual: float) -> str:
    return 'EXACT' if residual == 0 else 'ROUND_OFF' if residual <= ATOL else 'VIOLATED'


def radial_class(value: float) -> str:
    if abs(value) <= ATOL:
        return 'ZERO'
    return 'ATTRACTION' if value > 0 else 'REPULSION'


def temporal_class(values: np.ndarray) -> str:
    values = np.asarray(values)
    signs = np.sign(values[np.abs(values) > ATOL])
    if not len(signs):
        return 'ZERO'
    if np.any(signs != signs[0]):
        return 'SIGN_REVERSING'
    # A nonconstant native trace still has a stable radial polarity.
    return 'STEADY_ATTRACTION' if signs[0] > 0 else 'STEADY_REPULSION'


def symmetry_class(a: np.ndarray, b: np.ndarray) -> str:
    error = float(np.max(np.abs(np.asarray(a) - np.asarray(b))))
    return action_reaction_class(error)
