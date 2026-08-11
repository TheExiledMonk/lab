"""Diagnostic information geometry for the DEV177 native receipt audit.

These routines deliberately have no observer, projection, or fitting policy.
They describe finite matrices only; rank and SVD are never physical equations.
"""
from __future__ import annotations

import numpy as np


def _finite_rows(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    return x[np.isfinite(x).all(axis=1)]


def information_geometry(x: np.ndarray, *, rtol: float = 1e-10) -> dict:
    """Return unstandardized covariance/SVD diagnostics without imputing data."""
    x = _finite_rows(x)
    if len(x) < 2:
        return {"status": "INSUFFICIENT_SUPPORT", "n_rows": int(len(x))}
    centered = x - x.mean(axis=0, keepdims=True)
    singular = np.linalg.svd(centered, compute_uv=False)
    cutoff = (singular[0] * rtol) if len(singular) and singular[0] else 0.0
    rank = int(np.count_nonzero(singular > cutoff))
    power = singular * singular
    if power.sum():
        probabilities = power[power > 0] / power.sum()
        effective_rank = float(np.exp(-np.sum(probabilities * np.log(probabilities))))
    else:
        effective_rank = 0.0
    participation = float(power.sum() ** 2 / np.sum(power ** 2)) if np.any(power) else 0.0
    covariance = np.cov(centered, rowvar=False, ddof=1)
    scale = np.sqrt(np.outer(np.diag(covariance), np.diag(covariance)))
    correlation = np.divide(covariance, scale, out=np.full_like(covariance, np.nan), where=scale > 0)
    return {"status": "DEFINED", "n_rows": int(len(x)), "raw_dimensions": int(x.shape[1]),
            "numerical_rank": rank, "effective_rank": effective_rank,
            "participation_ratio": participation, "singular_values": singular,
            "covariance": covariance, "correlation": correlation,
            "standardized_for_diagnostic": False, "missing_rows_dropped_not_imputed": True}


def linear_recoverability(target: np.ndarray, predictors: np.ndarray) -> dict:
    """Centred least-squares residual, reported only as a redundancy diagnostic."""
    target = np.asarray(target, float); predictors = np.asarray(predictors, float)
    valid = np.isfinite(target).all(1) & np.isfinite(predictors).all(1)
    target, predictors = target[valid], predictors[valid]
    if len(target) < 2 or predictors.shape[1] == 0:
        return {"status": "INSUFFICIENT_SUPPORT"}
    yc = target - target.mean(0); xc = predictors - predictors.mean(0)
    fitted = xc @ np.linalg.lstsq(xc, yc, rcond=None)[0]
    denom = np.linalg.norm(yc)
    return {"status": "DEFINED", "relative_residual": float(np.linalg.norm(yc-fitted)/denom) if denom else 0.0,
            "explained_linear_fraction": float(1 - np.linalg.norm(yc-fitted)**2 / np.linalg.norm(yc)**2) if denom else 1.0}


def status_from_increment(increment: int, residual: float | None) -> str:
    if increment > 0:
        return "INDEPENDENT_INFORMATION"
    if residual is not None and residual > 0.1:
        return "COMPLEMENTARY_INFORMATION"
    if residual is not None and residual < 0.01:
        return "HIGHLY_REDUNDANT"
    return "PARTIALLY_REDUNDANT"
