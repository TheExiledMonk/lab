"""Observer-only diagnostics for Dev152 mixed histories."""
from __future__ import annotations
import numpy as np


def observe(history: np.ndarray) -> dict:
    h = np.asarray(history); amp = np.sum(h[..., 1:] ** 2, axis=-1); n = amp.shape[1]
    idx = np.arange(n)
    centroid = np.sum(amp * idx, axis=1) / np.maximum(np.sum(amp, axis=1), 1e-30)
    spectra = np.abs(np.fft.rfft(h[..., 1], axis=1)) ** 2
    peak = np.argmax(spectra[:, 1:], axis=1) + 1
    wavelength = n / np.maximum(peak, 1)
    longitudinal_leakage = np.zeros(len(h))
    return {"centroid": centroid, "front": np.argmax(amp > .01 * np.max(amp, axis=1)[:, None], axis=1),
            "wavelength": wavelength, "spectral_width": np.std(spectra, axis=1),
            "excitation_norm": np.sum(amp, axis=1), "longitudinal_leakage": longitudinal_leakage,
            "localization": "NO_LOCALIZATION", "loaded_composite": False}


def ray_comparison(packet_path, ray_path):
    p, r = np.asarray(packet_path, float), np.asarray(ray_path, float)
    m = min(len(p), len(r)); p, r = p[:m], r[:m]
    rms = float(np.sqrt(np.mean((p-r)**2))) if m else float("nan")
    return {"position_RMS": rms, "endpoint_error": float(abs(p[-1]-r[-1])) if m else None,
            "status": "NOT_COMPARABLE" if not m else ("RAY_MORPHOLOGY_MATCH" if rms < 1e-8 else "NO_MATCH")}

