"""Dev164 embedding gate: describe what the frozen scalar state can support."""
from __future__ import annotations

import numpy as np


def embedding_derivability(q: np.ndarray) -> dict:
    q = np.asarray(q, dtype=float)
    if q.ndim != 3:
        raise ValueError("static state must be a scalar 3D lattice")
    return {
        "candidate": "G00_SCALAR_STATE_ONLY",
        "scalar_node_excursion_available": True,
        "directed_scalar_differences_available": True,
        "deformed_bond_lengths_derivable": False,
        "deformed_bond_directions_derivable": False,
        "global_node_embedding_derivable": False,
        "reason": ("The frozen laws define scalar node excursion and scalar edge "
                   "differences, but no law maps either quantity to spatial bond "
                   "lengths, angles, or vector node displacement."),
        "candidate_displacement_u_times_normalized_asymmetry_accepted": False,
        "rejection": ("Node excursion is not established as displacement magnitude; "
                      "normalizing its scalar gradient would add a geometry law."),
    }
