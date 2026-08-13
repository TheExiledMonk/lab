"""Read-only local full-state descriptors for the frozen DEV167 lattice."""
from __future__ import annotations
import numpy as np
from pbuf.excitation.native_vector_pair_dynamics import directed_relations, pair_forces, pair_power_flux

def directed_forces(u: np.ndarray) -> np.ndarray:
    positive = pair_forces(u)
    out = []
    for axis in range(3):
        out.extend((positive[..., axis, :], -np.roll(positive[..., axis, :], 1, axis=axis)))
    return np.stack(out, axis=-2)

def descriptor(u: np.ndarray, p: np.ndarray, nodes: np.ndarray) -> dict[str, np.ndarray]:
    """Uncompressed native node momentum plus all six directed bond readouts."""
    r = directed_relations(u)[tuple(nodes.T)]
    length = np.linalg.norm(r, axis=-1)
    return {"momentum": p[tuple(nodes.T)], "relation": r, "strain": length - 1.0,
            "force": directed_forces(u)[tuple(nodes.T)],
            "power_flux": directed_power(u, p)[tuple(nodes.T)]}

def directed_power(u: np.ndarray, p: np.ndarray) -> np.ndarray:
    positive = pair_power_flux(u, p)
    out = []
    for axis in range(3):
        out.extend((positive[..., axis], -np.roll(positive[..., axis], 1, axis=axis)))
    return np.stack(out, axis=-1)

def full_features(state: dict[str, np.ndarray]) -> np.ndarray:
    """Lossless concatenation of descriptor levels only for recurrence distances."""
    return np.concatenate((state["momentum"].reshape(len(state["momentum"]), -1),
                           state["relation"].reshape(len(state["momentum"]), -1),
                           state["strain"].reshape(len(state["momentum"]), -1),
                           state["force"].reshape(len(state["momentum"]), -1),
                           state["power_flux"].reshape(len(state["momentum"]), -1)), axis=-1)

def scale_independent_recurrence(features: np.ndarray) -> np.ndarray:
    """RMS state distance / RMS temporal excursion; no selected tolerance."""
    centered = features - features.mean(axis=0, keepdims=True)
    scale = np.sqrt(np.mean(centered * centered))
    if scale == 0.0: return np.full(features.shape[0] - 1, np.nan)
    return np.array([np.sqrt(np.mean((features[k:] - features[:-k])**2)) / scale
                     for k in range(1, len(features))])
