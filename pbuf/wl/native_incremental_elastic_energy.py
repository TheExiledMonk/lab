"""Passive bounded-strain incremental-energy diagnostics (Dev143).

All quantities are native unless a caller supplies a physical ``cell_volume``.
Nothing in this module identifies medium elastic energy with photon energy.
"""
from __future__ import annotations

import numpy as np


def _finite(a, name: str):
    out = np.asarray(a, dtype=float)
    if np.any(~np.isfinite(out)):
        raise ValueError(f"{name} must be finite")
    return out


def bounded_strain_energy(strain, K: float = 1.0, epsilon_max: float = 1.0):
    """Return W(epsilon)=-K emax^2/2 log(1-(epsilon/emax)^2)."""
    e = _finite(strain, "strain")
    K, em = float(K), float(epsilon_max)
    if not np.isfinite(K) or K <= 0 or not np.isfinite(em) or em <= 0:
        raise ValueError("K and epsilon_max must be finite and positive")
    if np.any(np.abs(e) >= em):
        raise ValueError("bounded strain requires abs(strain) < epsilon_max")
    return -0.5 * K * em**2 * np.log1p(-(e / em) ** 2)


def bounded_strain_stress(strain, K: float = 1.0, epsilon_max: float = 1.0):
    e = _finite(strain, "strain")
    if np.any(np.abs(e) >= epsilon_max):
        raise ValueError("bounded strain requires abs(strain) < epsilon_max")
    return K * e / (1.0 - (e / epsilon_max) ** 2)


def bounded_strain_tangent(strain, K: float = 1.0, epsilon_max: float = 1.0):
    """Exact W'' used by the Taylor audit."""
    e = _finite(strain, "strain"); x = (e / epsilon_max) ** 2
    if np.any(x >= 1):
        raise ValueError("bounded strain requires abs(strain) < epsilon_max")
    return K * (1.0 + x) / (1.0 - x) ** 2


def incremental_elastic_energy(background_strain, perturbation_strain, K: float = 1.0,
                               epsilon_max: float = 1.0):
    """Exact Delta W = W(background + perturbation) - W(background)."""
    bg, de = np.broadcast_arrays(_finite(background_strain, "background_strain"),
                                 _finite(perturbation_strain, "perturbation_strain"))
    return bounded_strain_energy(bg + de, K, epsilon_max) - bounded_strain_energy(bg, K, epsilon_max)


delta_W = incremental_elastic_energy


def taylor_increment(background_strain, perturbation_strain, K: float = 1.0,
                     epsilon_max: float = 1.0, order: int = 2):
    bg, de = np.broadcast_arrays(_finite(background_strain, "background_strain"),
                                 _finite(perturbation_strain, "perturbation_strain"))
    out = bounded_strain_stress(bg, K, epsilon_max) * de
    if order >= 2:
        out = out + 0.5 * bounded_strain_tangent(bg, K, epsilon_max) * de**2
    if order not in (1, 2):
        raise ValueError("only first- and second-order audits are supported")
    return out


def integrate_packet(delta_w, mask=None, cell_volume=1.0):
    """Return signed, positive, and absolute packet integrals."""
    dw = _finite(delta_w, "delta_w")
    selected = np.ones(dw.shape, bool) if mask is None else np.asarray(mask, bool)
    if selected.shape != dw.shape:
        raise ValueError("mask shape must match delta_w")
    vol = np.broadcast_to(_finite(cell_volume, "cell_volume"), dw.shape)[selected]
    if np.any(vol <= 0):
        raise ValueError("cell_volume must be positive")
    values = dw[selected] * vol
    return {"signed": float(values.sum()), "positive": float(np.maximum(values, 0).sum()),
            "absolute": float(np.abs(values).sum()), "cell_count": int(values.size),
            "native_cell_volume_normalized": bool(np.asarray(cell_volume).ndim == 0 and float(cell_volume) == 1.0)}


def excitation_energy(background_strain, perturbation_strain, K=1.0, epsilon_max=1.0):
    """Convex Bregman excess: Delta W - sigma(background)*delta epsilon."""
    return (incremental_elastic_energy(background_strain, perturbation_strain, K, epsilon_max)
            - bounded_strain_stress(background_strain, K, epsilon_max) * np.asarray(perturbation_strain))


def positivity_audit(background_strain, perturbation_strain, K=1.0, epsilon_max=1.0):
    dw = incremental_elastic_energy(background_strain, perturbation_strain, K, epsilon_max)
    ex = excitation_energy(background_strain, perturbation_strain, K, epsilon_max)
    return {"signed_increment_can_be_negative": bool(np.any(dw < 0)),
            "positive_excitation_nonnegative": bool(np.all(ex >= -1e-12)),
            "minimum_signed": float(np.min(dw)), "minimum_excitation": float(np.min(ex))}
