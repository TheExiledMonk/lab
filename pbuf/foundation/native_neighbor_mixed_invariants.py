"""Dev152 diagnostic joint-invariant and exchange audit."""
from __future__ import annotations
import numpy as np
from pbuf.wl.native_incremental_elastic_energy import bounded_strain_energy

JOINT_CANDIDATES = tuple(f"J{i:02d}" for i in range(1, 11))


def candidate_values(history: np.ndarray) -> dict[str, np.ndarray]:
    h = np.asarray(history, float); L, X = h[..., 0], h[..., 1:]
    x2 = np.sum(X * X, axis=(-1, -2)); l2 = np.sum(L * L, axis=-1)
    w = np.sum(bounded_strain_energy(L), axis=-1)
    values = {"J01": l2 + x2, "J02": w + x2, "J03": w + x2,
              "J04": l2 + x2, "J05": np.sum(np.tanh(L) ** 2, axis=-1) + x2,
              "J06": np.sqrt(l2 + x2), "J07": l2 * x2,
              "J08": l2 + x2, "J09": np.sqrt(l2 + x2), "J10": np.full(len(h), np.nan)}
    return values


def audit(history: np.ndarray, tolerance=1e-10) -> dict:
    values = candidate_values(history); rows = []
    for name, v in values.items():
        if name == "J10": rows.append({"candidate": name, "status": "NO_JOINT_INVARIANT"}); continue
        drift = float(np.ptp(v) / max(abs(v[0]), 1e-30))
        rows.append({"candidate": name, "relative_drift": drift,
                     "conserved": drift <= tolerance, "status": "PASS" if drift <= tolerance else "FAIL"})
    return {"candidates": rows, "conserved_candidates": [r["candidate"] for r in rows if r.get("conserved")]}


def norm_exchange(history: np.ndarray) -> dict:
    h = np.asarray(history); nl = np.sum(h[..., 0] ** 2, axis=1); nx = np.sum(h[..., 1:] ** 2, axis=(1, 2))
    dl, dx = float(nl[-1] - nl[0]), float(nx[-1] - nx[0])
    cls = "NO_EXCHANGE" if abs(dl) + abs(dx) < 1e-10 else "REVERSIBLE_EXCHANGE"
    return {"delta_N_L": dl, "delta_N_X": dx, "classification": cls}

