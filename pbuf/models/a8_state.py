"""Frozen A8/T1 state — slow + fast scalar fields, evolved under the
frozen transport update. Extracted from a8_three_dimensional_projection_lab001.

The A8 state is consumed by:
  * a8_pair_amplitude.compute_a8_pair_amplitudes
  * a8_pair_amplitude.compute_longitudinal_axis
"""
from __future__ import annotations
import numpy as np

from ..core.conventions import EPS_FLOAT, EPS_ZERO

__all__ = [
    "A8_INIT_DT", "A8_INIT_OMEGA", "A8_INIT_K", "A8_INIT_STEPS",
    "A8_INIT_INJECTION_NOISE", "A8_INIT_COUP_F2S", "A8_INIT_COUP_S2F",
    "A8_INIT_FAST_TIMESCALE", "A8_INIT_SLOW_TIMESCALE",
    "A8_INIT_CLIP",
    "neighbours6_face_reflective_3d",
    "evolve_a8_transport_3d",
    "build_a8_state_3d",
]


# These constants reproduce the frozen A8/T1 update used by the
# previous correction lab. They are referenced from a8_pair_amplitude
# and frozen here so the state module owns them.
A8_INIT_DT = 0.03
A8_INIT_OMEGA = 1.0
A8_INIT_K = 1.0
A8_INIT_STEPS = 160
A8_INIT_INJECTION_NOISE = 0.02
A8_INIT_COUP_F2S = 0.5
A8_INIT_COUP_S2F = 0.5
A8_INIT_FAST_TIMESCALE = 1.0
A8_INIT_SLOW_TIMESCALE = 0.1
A8_INIT_CLIP = 5.0


def neighbours6_face_reflective_3d(u: np.ndarray) -> np.ndarray:
    """Average of 6 face-connected voxels with reflective padding."""
    p = np.pad(u, ((1, 1), (1, 1), (1, 1)), mode="reflect")
    n_xm = p[1:-1, 1:-1, :-2]
    n_xp = p[1:-1, 1:-1, 2:]
    n_ym = p[1:-1, :-2, 1:-1]
    n_yp = p[1:-1, 2:, 1:-1]
    n_zm = p[:-2, 1:-1, 1:-1]
    n_zp = p[2:, 1:-1, 1:-1]
    return (n_xm + n_xp + n_ym + n_yp + n_zm + n_zp) / 6.0


def evolve_a8_transport_3d(u_slow, u_fast, stencil="N6", boundary="reflective"):
    """Run the frozen T1 update for STEPS timesteps.

    Returns (u_slow_final, u_fast_final, history).
    """
    history = []
    history.append(0.5 * (u_slow + u_fast))

    if stencil != "N6" or boundary != "reflective":
        raise ValueError(f"unsupported stencil/boundary: {stencil}/{boundary}")
    n6 = neighbours6_face_reflective_3d

    for _ in range(A8_INIT_STEPS):
        n_slow = n6(u_slow)
        n_fast = n6(u_fast)
        d_fast = A8_INIT_DT * A8_INIT_OMEGA * A8_INIT_K * (
            (n_fast - u_fast)
            + A8_INIT_COUP_S2F * (u_slow - u_fast))
        d_slow = A8_INIT_DT * A8_INIT_SLOW_TIMESCALE * (
            (n_slow - u_slow)
            + A8_INIT_COUP_F2S * (u_fast - u_slow))
        u_fast = np.clip(u_fast + d_fast, -A8_INIT_CLIP, A8_INIT_CLIP)
        u_slow = np.clip(u_slow + d_slow, -A8_INIT_CLIP, A8_INIT_CLIP)
        history.append(0.5 * (u_slow + u_fast))

    return u_slow, u_fast, history


def build_a8_state_3d(rho_3d, strength, seed=12345):
    """Initialise and evolve the A8/T1 state for a 3D mass-density field.

    Returns dict with keys rho_3d, u_slow, u_fast, c_state (the combined
    state at the final timestep).
    """
    rng = np.random.RandomState(seed)
    eq = strength * rho_3d
    u_slow = eq.copy()
    u_fast = eq.copy() + A8_INIT_INJECTION_NOISE * strength * rng.randn(*rho_3d.shape)
    u_slow, u_fast, history = evolve_a8_transport_3d(u_slow, u_fast)
    c_state = history[-1]
    return {"rho_3d": rho_3d.copy(),
            "u_slow": u_slow.copy(),
            "u_fast": u_fast.copy(),
            "c_state": c_state.copy()}


if __name__ == "__main__":
    nz, ny, nx = 5, 8, 8
    rho = np.exp(-((np.arange(nz)[:, None, None] - 2) ** 2 +
                   (np.arange(ny)[None, :, None] - 4) ** 2 +
                   (np.arange(nx)[None, None, :] - 4) ** 2) / 4.0)
    state = build_a8_state_3d(rho, strength=0.18, seed=12345)
    assert state["u_slow"].shape == rho.shape
    assert state["u_fast"].shape == rho.shape
    assert state["c_state"].shape == rho.shape
    # Non-trivial: at least one nonzero voxel in u_fast - u_slow.
    diff = state["u_fast"] - state["u_slow"]
    assert np.max(np.abs(diff)) > 0.0
    print(f"M06 A8 state built: shape={rho.shape}, "
          f"max|u_f-u_s|={float(np.max(np.abs(diff))):.4e}")