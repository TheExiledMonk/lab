"""External quantum bridge mathematics for a zero-rest-mass mode.

This module deliberately contains no PBUF constitutive physics.  SI-valued
absolute transformations and normalization-free ratio transformations are kept
separate so a successful algebra check cannot be mistaken for a native state.
"""
from __future__ import annotations

import numpy as np

BRIDGE_ROLE = "EXTERNAL_ESTABLISHED_QUANTUM_RELATION"
HBAR_SI = 1.054_571_817e-34
C_SI = 299_792_458.0
H_SI = 2.0 * np.pi * HBAR_SI


def _positive(value, name: str):
    a = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(a)) or np.any(a <= 0):
        raise ValueError(f"{name} must be finite and positive")
    return a


def energy_from_momentum(momentum, c: float = C_SI):
    return _positive(momentum, "momentum") * _positive(c, "c")


def momentum_from_energy(energy, c: float = C_SI):
    return _positive(energy, "energy") / _positive(c, "c")


def wave_number_from_momentum(momentum, hbar: float = HBAR_SI):
    return _positive(momentum, "momentum") / _positive(hbar, "hbar")


def wave_number_from_energy(energy, hbar: float = HBAR_SI, c: float = C_SI):
    return _positive(energy, "energy") / (_positive(hbar, "hbar") * _positive(c, "c"))


def wavelength_from_energy(energy, h: float = H_SI, c: float = C_SI):
    return _positive(h, "h") * _positive(c, "c") / _positive(energy, "energy")


def ratio_bridge(energy_ratio):
    """Map E/E0 to p/p0, k/k0, lambda/lambda0 and redshift."""
    r = _positive(energy_ratio, "energy_ratio")
    return {"energy_ratio": r, "momentum_ratio": r.copy(), "k_ratio": r.copy(),
            "wavelength_ratio": 1.0 / r, "one_plus_z": 1.0 / r,
            "redshift": 1.0 / r - 1.0, "bridge_role": BRIDGE_ROLE}


def strain_mode_ratio_bridge(mode_energy_ratio, *, proxy_established: bool):
    """Consume a gated relative strain-mode-energy proxy, never infer one."""
    if not proxy_established:
        raise ValueError("an established mode-energy proxy is required")
    out = ratio_bridge(mode_energy_ratio)
    out["proxy_status"] = "SUPPORTED_PROXY"
    return out


def bridge_contract(**availability):
    base = {"contract": "PBUF_QUANTUM_ZERO_MASS_BRIDGE_V1",
            "p_equals_hbar_k_used": True, "E_equals_pc_used": True,
            "E_equals_hbar_c_k_used": True, "bridge_role": BRIDGE_ROLE,
            "absolute_k_established": False, "relative_k_established": False,
            "absolute_lambda_established": False, "relative_lambda_established": False,
            "phase_required": False, "wavefunction_required": False,
            "native_time_required": False}
    base.update(availability)
    return base
