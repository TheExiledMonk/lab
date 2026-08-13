"""DEV212 signed relational-rotation artifact helpers (no new angular variable)."""
from __future__ import annotations
import numpy as np
from pbuf.observer.native_internal_state_inventory import axial_dual


def sign_inventory(antisymmetric_trace: np.ndarray) -> dict:
    q = axial_dual(antisymmetric_trace)
    return {"antisymmetric_nonzero": bool(np.any(antisymmetric_trace != 0)),
            "component_has_both_signs": bool(np.any(q > 0) and np.any(q < 0)),
            "axial_dual_l2": float(np.linalg.norm(q)),
            "axial_dual_min": float(np.min(q)), "axial_dual_max": float(np.max(q))}
