"""DEV212 gate-preserving empty interaction-matrix representation."""
from __future__ import annotations
import numpy as np


def blocked_matrices() -> dict[str, np.ndarray]:
    """No pair composition is implied when the independently justified gate is closed."""
    return {"state_labels": np.asarray([], dtype="U1"), "radial_force": np.empty((0, 0)),
            "torque": np.empty((0, 0, 3)), "interaction_residual": np.empty((0, 0, 3))}
