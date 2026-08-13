"""DEV213 pre-evolution invariant summaries; no force or torque diagnostic."""
from __future__ import annotations
import numpy as np
from pbuf.excitation.native_vector_pair_dynamics import invariant, pair_power_flux, potential


def summary(state) -> dict:
    kinetic = float(0.5 * np.sum(state.momentum ** 2))
    return {
        "total_energy": invariant(state.displacement, state.momentum),
        "potential_energy": potential(state.displacement),
        "kinetic_energy": kinetic,
        "total_momentum": np.sum(state.momentum, axis=(0, 1, 2)),
        "pair_power_flux": np.sum(pair_power_flux(state.displacement, state.momentum), axis=(0, 1, 2)),
    }
