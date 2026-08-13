"""Algebraic L/T reporting for DEV201 native tangent eigenvectors."""
from __future__ import annotations

import numpy as np


def longitudinal_transverse(vector: np.ndarray, k: np.ndarray) -> dict:
    u = np.asarray(vector, dtype=complex)
    k = np.asarray(k, dtype=float)
    nk = float(np.linalg.norm(k)); nu = float(np.linalg.norm(u))
    if nk == 0.0 or nu == 0.0:
        return {"f_L": None, "f_T": None, "sector": "ZERO_MODE"}
    khat = k / nk
    ul = np.vdot(khat, u)
    fl = float(abs(ul) ** 2 / (nu * nu)); ft = float(1.0 - fl)
    exact_l = np.isclose(ft, 0.0, atol=1e-12)
    exact_t = np.isclose(fl, 0.0, atol=1e-12)
    return {"f_L": fl, "f_T": ft,
            "sector": "EXACT_LONGITUDINAL" if exact_l else "EXACT_TRANSVERSE" if exact_t else "MIXED"}


def transverse_rank(vectors: np.ndarray, k: np.ndarray) -> int:
    """Exact numerical rank of transverse projections of a candidate eigenspace."""
    k = np.asarray(k, float); nk = np.linalg.norm(k)
    if nk == 0.0: return 0
    P = np.eye(3) - np.outer(k, k) / (nk * nk)
    return int(np.linalg.matrix_rank(P @ np.asarray(vectors), tol=1e-11))
