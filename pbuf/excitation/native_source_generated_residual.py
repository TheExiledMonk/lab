"""Spectrum and frozen-dispersion diagnostics for source-generated states."""
from __future__ import annotations

import numpy as np


def analytic_omega_grid(shape: tuple[int, int, int]) -> np.ndarray:
    axes = [2*np.pi*np.fft.fftfreq(n) for n in shape]
    k = np.meshgrid(*axes, indexing="ij")
    rhs = 1.0 - sum(np.sin(x/2.0)**2 for x in k)/3.0
    return np.arccos(np.clip(rhs, -1.0, 1.0))


def generated_spectrum(states: np.ndarray) -> dict:
    states = np.asarray(states, float)
    if states.ndim != 4:
        raise ValueError("states must have shape (time,Nx,Ny,Nz)")
    spatial = np.fft.fftn(states, axes=(1,2,3))
    power = np.mean(np.abs(spatial)**2, axis=0)
    power[(0,0,0)] = 0
    dominant = tuple(int(x) for x in np.unravel_index(np.argmax(power), power.shape))
    temporal = np.fft.fft(spatial, axis=0)
    joint = np.abs(temporal)**2
    frequency_index = int(np.argmax(joint[:, dominant[0], dominant[1], dominant[2]]))
    omega = float(abs(2*np.pi*np.fft.fftfreq(states.shape[0])[frequency_index]))
    exact = float(analytic_omega_grid(states.shape[1:])[dominant])
    total = float(np.sum(power))
    participation = float(total**2/np.sum(power**2)) if total and np.sum(power**2) else 0.0
    axes = [2*np.pi*np.fft.fftfreq(n) for n in states.shape[1:]]
    grid = np.meshgrid(*axes, indexing="ij")
    k2 = sum(x*x for x in grid)
    rms_k = float(np.sqrt(np.sum(k2*power)/total)) if total else 0.0
    directional = [float(np.sum((grid[i]*grid[i])*power)/np.sum(k2*power))
                   if np.sum(k2*power) else 0.0 for i in range(3)]
    return {"dominant_mode_indices": list(dominant), "measured_progression_frequency": omega,
            "dev157_progression_frequency": exact, "absolute_frequency_error": abs(omega-exact),
            "spectral_participation_modes": participation, "total_spectral_power": total,
            "rms_native_wave_number": rms_k,
            "rms_native_wavelength": float(2*np.pi/rms_k) if rms_k else None,
            "directional_k2_fractions": directional,
            "spectral_width_measure":"RMS_NATIVE_WAVE_NUMBER"}


def dispersion_match(states: np.ndarray) -> dict:
    """Compare all energetic space-time bins with the exact Dev157 branches."""
    states = np.asarray(states, float)
    transformed = np.fft.fftn(states, axes=(0,1,2,3))
    power = np.abs(transformed)**2
    power[:,0,0,0] = 0
    omega_bins = abs(2*np.pi*np.fft.fftfreq(states.shape[0]))
    exact = analytic_omega_grid(states.shape[1:])
    nearest_error = np.min(abs(omega_bins[:,None,None,None]-exact[None,...]), axis=0)
    spatial_power = power.sum(axis=0)
    weighted_error = float(np.sum(nearest_error*spatial_power)/np.sum(spatial_power)) if spatial_power.sum() else 0.0
    resolution = float(2*np.pi/states.shape[0])
    return {"weighted_branch_error": weighted_error, "frequency_bin_width": resolution,
            "within_one_frequency_bin": bool(weighted_error <= resolution),
            "classification": "TRUE" if weighted_error <= resolution else "PARTIAL"}
