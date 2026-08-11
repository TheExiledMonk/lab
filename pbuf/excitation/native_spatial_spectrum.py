"""Canonical three-dimensional lattice spectrum diagnostics for Dev157."""
from __future__ import annotations

import numpy as np


def native_wavevectors(shape: tuple[int, int, int]) -> tuple[np.ndarray, ...]:
    """Return FFT-ordered native wave-number axes in radians per lattice cell."""
    if len(shape) != 3 or any(int(n) < 2 for n in shape):
        raise ValueError("shape must contain three lattice extents >= 2")
    return tuple(2.0 * np.pi * np.fft.fftfreq(int(n)) for n in shape)


def spectrum3d(field: np.ndarray) -> dict[str, np.ndarray | float]:
    """Unitary 3-D Fourier transform, power, and canonical |k| grid."""
    a = np.asarray(field, dtype=float)
    if a.ndim != 3 or not np.isfinite(a).all():
        raise ValueError("field must be a finite three-dimensional array")
    transform = np.fft.fftn(a, norm="ortho")
    axes = native_wavevectors(a.shape)
    mesh = np.meshgrid(*axes, indexing="ij")
    kmagnitude = np.sqrt(sum(component * component for component in mesh))
    return {
        "transform": transform,
        "power": np.abs(transform) ** 2,
        "k_axes": axes,
        "k_magnitude": kmagnitude,
        "dc_power": float(np.abs(transform[(0, 0, 0)]) ** 2),
    }


def reconstruct(transform: np.ndarray) -> np.ndarray:
    """Inverse of the unitary transform (real fields only)."""
    return np.fft.ifftn(np.asarray(transform), norm="ortho").real


def radial_spectrum(power: np.ndarray, k_magnitude: np.ndarray, bins: int | None = None) -> dict:
    """Shell-sum power without silently excluding the DC mode."""
    p = np.asarray(power, float)
    k = np.asarray(k_magnitude, float)
    if p.shape != k.shape or p.ndim != 3 or np.any(p < 0):
        raise ValueError("power and k_magnitude must be matching nonnegative 3-D arrays")
    bins = int(bins or max(p.shape) // 2 + 1)
    edges = np.linspace(0.0, float(k.max()) + np.finfo(float).eps, bins + 1)
    which = np.minimum(np.digitize(k.ravel(), edges) - 1, bins - 1)
    shell_power = np.bincount(which, weights=p.ravel(), minlength=bins)
    counts = np.bincount(which, minlength=bins)
    centers = 0.5 * (edges[:-1] + edges[1:])
    non_dc = k > 0
    total_non_dc = float(p[non_dc].sum())
    centroid = (float((p[non_dc] * k[non_dc]).sum()) / total_non_dc
                if total_non_dc else None)
    width = (float(np.sqrt((p[non_dc] * (k[non_dc] - centroid) ** 2).sum()
                           / total_non_dc)) if total_non_dc else None)
    eligible = np.flatnonzero(centers > 0)
    peak_index = int(eligible[np.argmax(shell_power[eligible])]) if len(eligible) else 0
    return {
        "k_centers": centers.tolist(), "shell_power": shell_power.tolist(),
        "mode_counts": counts.tolist(), "dc_power": float(p[k == 0].sum()),
        "dc_explicitly_included": True, "non_dc_peak_k": float(centers[peak_index]),
        "non_dc_spectral_centroid": centroid, "non_dc_spectral_width": width,
    }


def native_wavelength(k: float) -> float | None:
    """Native wavelength in lattice cells; DC has no finite wavelength."""
    return None if float(k) == 0.0 else float(2.0 * np.pi / abs(k))
