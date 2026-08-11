"""Progression-frequency diagnostics for the frozen Dev156 maps."""
from __future__ import annotations

import numpy as np

from .native_relational_dynamics import f02_step, f03_step


def mode(shape: tuple[int, int, int], indices: tuple[int, int, int], amplitude=1.0,
         phase=0.0) -> np.ndarray:
    """Real lattice Fourier mode specified by integer FFT indices."""
    if len(shape) != 3 or len(indices) != 3:
        raise ValueError("shape and indices must be three-dimensional")
    grids = np.meshgrid(*(np.arange(n) for n in shape), indexing="ij")
    argument = sum(2.0 * np.pi * int(m) * x / int(n)
                   for x, m, n in zip(grids, indices, shape)) + phase
    return float(amplitude) * np.cos(argument)


def wavevector(shape: tuple[int, int, int], indices: tuple[int, int, int]) -> np.ndarray:
    signed = [m if m <= n // 2 else m - n for m, n in zip(indices, shape)]
    return 2.0 * np.pi * np.asarray(signed, float) / np.asarray(shape, float)


def analytic_progression_frequency(kvector: np.ndarray) -> float:
    """Exact Dev156 lattice dispersion in radians per progression step."""
    k = np.asarray(kvector, float)
    cosine = 1.0 - np.sum(np.sin(k / 2.0) ** 2) / 3.0
    return float(np.arccos(np.clip(cosine, -1.0, 1.0)))


def _coefficient(field: np.ndarray, reference: np.ndarray) -> float:
    return float(np.sum(field * reference) / np.sum(reference * reference))


def measure_mode_frequency(representation: str, shape: tuple[int, int, int],
                           indices: tuple[int, int, int], amplitude=1e-3,
                           steps=24) -> dict:
    """Measure frequency from the exact second-order modal recurrence."""
    reference = mode(shape, indices, amplitude=1.0)
    q = amplitude * reference
    coefficients = [_coefficient(q, reference)]
    if representation == "F02":
        auxiliary = np.zeros(q.shape + (3,))
        step = f02_step
    elif representation == "F03":
        auxiliary = np.zeros_like(q)
        step = f03_step
    else:
        raise ValueError("representation must be F02 or F03")
    for _ in range(steps):
        q, auxiliary = step(q, auxiliary)
        coefficients.append(_coefficient(q, reference))
    c = np.asarray(coefficients)
    denominator = 2.0 * c[1:-1]
    valid = np.abs(denominator) > max(abs(amplitude), 1.0) * 1e-14
    cosine_samples = (c[2:] + c[:-2])[valid] / denominator[valid]
    measured = float(np.arccos(np.clip(np.median(cosine_samples), -1.0, 1.0)))
    kv = wavevector(shape, indices)
    exact = analytic_progression_frequency(kv)
    return {
        "mode_indices": list(map(int, indices)), "k_vector": kv.tolist(),
        "k_magnitude": float(np.linalg.norm(kv)),
        "native_wavelength": None if np.linalg.norm(kv) == 0 else float(2*np.pi/np.linalg.norm(kv)),
        "progression_frequency": measured, "analytic_progression_frequency": exact,
        "absolute_frequency_error": abs(measured - exact),
        "phase_progression": None if np.linalg.norm(kv) == 0 else measured / float(np.linalg.norm(kv)),
        "amplitude": float(amplitude), "steps": int(steps),
    }


def radial_group_progression(k: float, direction=(1.0, 0.0, 0.0)) -> float:
    """Analytic radial derivative dOmega/dk along a specified direction."""
    d = np.asarray(direction, float); d /= np.linalg.norm(d)
    kv = float(k) * d
    omega = analytic_progression_frequency(kv)
    if omega == 0.0:
        return float(1.0 / np.sqrt(6.0))
    derivative_cos = -np.sum(np.sin(kv) * d) / 6.0
    return float(-derivative_cos / np.sin(omega))
