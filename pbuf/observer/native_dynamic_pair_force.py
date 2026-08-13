"""DEV214 force diagnostics for immutable DEV213 preparation provenance."""
from __future__ import annotations

import numpy as np


def provenance_mask(packet_displacement: np.ndarray) -> np.ndarray:
    """The frozen DEV213 transverse packet footprint, extended over its native x envelope.

    This is preparation bookkeeping, not an evolved structure-identification rule.
    """
    return np.any(np.asarray(packet_displacement)[1] != 0.0, axis=-1)[None, ...].repeat(packet_displacement.shape[0], axis=0)


def support_force(node_force: np.ndarray, mask: np.ndarray) -> np.ndarray:
    return np.sum(np.asarray(node_force)[mask], axis=0)


def classify(force_a: np.ndarray, force_b: np.ndarray, rhat: np.ndarray, atol: float) -> tuple[str, float, float]:
    ar, br = float(force_a @ rhat), float(-force_b @ rhat)
    if np.linalg.norm(force_a + force_b) > atol:
        return "UNRESOLVED", ar, br
    if abs(ar) <= atol and abs(br) <= atol:
        return "ZERO", ar, br
    if ar > 0 and br < 0:
        return "ATTRACTION", ar, br
    if ar < 0 and br > 0:
        return "REPULSION", ar, br
    return "MIXED", ar, br


def temporal_character(values: np.ndarray, atol: float) -> str:
    x = np.asarray(values)
    signs = np.sign(x[np.abs(x) > atol])
    if not len(signs): return "ZERO"
    if np.all(signs == signs[0]): return "STEADY_SIGN"
    return "SIGN_REVERSING"
