"""Full-3D source-sheet differential diagnostics for DEV177."""
from __future__ import annotations

import numpy as np


def fit_j3(source: np.ndarray, received: np.ndarray, *, condition_limit: float = 1e10) -> dict:
    """Fit x_r = a + J3 [s1,s2] from a genuine rank-two source sheet.

    Source coordinates are the two varying coordinates discovered from the
    source positions; no missing derivative is represented as zero.
    """
    source = np.asarray(source, float); received = np.asarray(received, float)
    valid = np.isfinite(source).all(1) & np.isfinite(received).all(1)
    source, received = source[valid], received[valid]
    unique, inverse = np.unique(source, axis=0, return_inverse=True)
    # Multiple receipt events belonging to one source are reduced to its mean
    # only for this differential fit; individual receipts remain primary data.
    mean_received = np.array([received[inverse == i].mean(0) for i in range(len(unique))])
    varying = np.ptp(unique, axis=0) > 1e-12
    coords = unique[:, varying]
    result = {"N_RECEIPTS": int(len(source)), "N_UNIQUE_SOURCE_POINTS": int(len(unique)),
              "SOURCE_RANK": int(np.linalg.matrix_rank(coords-coords.mean(0))) if len(coords) else 0}
    if coords.shape[1] != 2 or result["SOURCE_RANK"] < 2 or len(unique) < 3:
        result["J3_STATUS"] = "UNDEFINED_INSUFFICIENT_SOURCE_SUPPORT"
        return result
    a = np.column_stack((np.ones(len(coords)), coords-coords.mean(0)))
    condition = float(np.linalg.cond(a[:, 1:]))
    result["CONDITION_NUMBER"] = condition
    if not np.isfinite(condition) or condition > condition_limit:
        result["J3_STATUS"] = "UNDEFINED_INSUFFICIENT_SOURCE_SUPPORT"
        return result
    coefficients, _, _, _ = np.linalg.lstsq(a, mean_received, rcond=None)
    predicted = a @ coefficients
    j3 = coefficients[1:].T
    metric = j3.T @ j3
    result.update({"J3_STATUS": "DEFINED", "J3": j3, "G3": metric,
                   "FIT_RESIDUAL": float(np.sqrt(np.mean((predicted-mean_received)**2))),
                   "trace": float(np.trace(metric)), "determinant": float(np.linalg.det(metric)),
                   "area_element": float(np.sqrt(max(np.linalg.det(metric), 0.0)))})
    vals, vecs = np.linalg.eigh(metric)
    result["principal_stretches"] = np.sqrt(np.maximum(vals, 0))
    result["anisotropy_ratio"] = float(np.sqrt(vals[-1]/vals[0])) if vals[0] > 0 else None
    result["symmetric_traceless"] = metric - np.eye(2)*np.trace(metric)/2
    result["principal_orientation_source_sheet"] = vecs[:, -1]
    return result
