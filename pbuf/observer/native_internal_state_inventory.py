"""Read-only full-state diagnostics used by the DEV212 multistate audit."""
from __future__ import annotations

import numpy as np

from pbuf.excitation.native_vector_pair_dynamics import (
    VectorPairState, invariant, pair_power_flux, potential, step,
)
from pbuf.observer.relational_wave_content import native_symmetric_antisymmetric
from pbuf.excitation.native_vector_pair_dynamics import directed_relations


def reverse_momentum(state: VectorPairState) -> VectorPairState:
    """The exact state involution R_p(u,p)=(u,-p)."""
    return VectorPairState(state.displacement.copy(), -state.momentum.copy(), state.progression_step)


def axial_dual(antisymmetric: np.ndarray) -> np.ndarray:
    """Coefficient-free axial dual Q_i=epsilon_ijk A_jk/2, as a diagnostic only."""
    return np.stack((antisymmetric[..., 1, 2] - antisymmetric[..., 2, 1],
                     antisymmetric[..., 2, 0] - antisymmetric[..., 0, 2],
                     antisymmetric[..., 0, 1] - antisymmetric[..., 1, 0]), axis=-1) / 2


def state_summary(state: VectorPairState) -> dict:
    rel = directed_relations(state.displacement)
    _, anti = native_symmetric_antisymmetric(rel)
    flux = pair_power_flux(state.displacement, state.momentum)
    return {
        "total_energy": invariant(state.displacement, state.momentum),
        "potential_energy": potential(state.displacement),
        "kinetic_energy": float(0.5 * np.sum(state.momentum**2)),
        "displacement_l2": float(np.linalg.norm(state.displacement)),
        "momentum_l2": float(np.linalg.norm(state.momentum)),
        "flux_sum": np.sum(flux, axis=(0, 1, 2)),
        "antisymmetric_l2": float(np.linalg.norm(anti)),
        "axial_dual_sum": np.sum(axial_dual(anti), axis=(0, 1, 2)),
    }


def persistence_trace(state: VectorPairState, updates: int, numerical_step: float) -> dict[str, np.ndarray]:
    """Unchanged DEV167 free evolution; no source, damping, or new diagnostic law."""
    cur = state
    displacement, momentum, energy, flux = [], [], [], []
    for _ in range(updates + 1):
        s = state_summary(cur)
        displacement.append(s["displacement_l2"])
        momentum.append(s["momentum_l2"])
        energy.append(s["total_energy"])
        flux.append(s["flux_sum"])
        cur = step(cur, numerical_step)
    return {"displacement_l2": np.asarray(displacement), "momentum_l2": np.asarray(momentum),
            "total_energy": np.asarray(energy), "flux_sum": np.asarray(flux)}
