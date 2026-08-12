"""Exact finite-kernel moment identities for DEV193.

These functions describe a nonnegative *mathematical mixture* of frozen
single-packet columns.  They deliberately make no claim that the mixture is a
simultaneous physical multi-packet state.
"""
from __future__ import annotations

import numpy as np


def input_moments(launch_coordinates: np.ndarray, content: np.ndarray):
    x, s = np.asarray(launch_coordinates, float), np.asarray(content, float)
    total = float(s.sum())
    if total <= 0:
        return np.full(2, np.nan), np.full((2, 2), np.nan)
    mean = (s[:, None] * x).sum(axis=0) / total
    dx = x - mean
    return mean, (dx.T * s) @ dx / total


def kernel_moment_decomposition(weight: np.ndarray, detector_coordinates: np.ndarray,
                                content: np.ndarray):
    """Return exact output covariance as centroid spread plus response spread.

    The transported launch weights are ``w_l = eta_l s_l / sum(eta*s)``.  This
    is simply finite-sum algebra over the positive kernel, not a probabilistic
    or optical assumption.
    """
    k, q, s = np.asarray(weight, float), np.asarray(detector_coordinates, float), np.asarray(content, float)
    eta = k.sum(axis=0)
    received = eta * s
    total = float(received.sum())
    if total <= 0:
        nan2 = np.full(2, np.nan); nan22 = np.full((2, 2), np.nan)
        return {"throughput": eta, "transported_weights": np.full_like(s, np.nan),
                "output_centroid": nan2, "output_covariance": nan22,
                "centroid_covariance": nan22, "response_covariance": nan22}
    p = np.divide(k, eta[None, :], out=np.zeros_like(k), where=eta[None, :] > 0)
    y = p.T @ q
    response = np.empty((len(s), 2, 2))
    for i in range(len(s)):
        dq = q - y[i]
        response[i] = (dq.T * p[:, i]) @ dq if eta[i] > 0 else np.nan
    w = received / total
    out = k @ s
    out_mean = (out[:, None] * q).sum(axis=0) / total
    dq = q - out_mean
    output_cov = (dq.T * out) @ dq / total
    dy = y - out_mean
    centroid_cov = (dy.T * w) @ dy
    response_cov = np.einsum("l,lij->ij", w, response)
    return {"throughput": eta, "conditional_centroids": y, "conditional_covariances": response,
            "transported_weights": w, "output_centroid": out_mean,
            "output_covariance": output_cov, "centroid_covariance": centroid_cov,
            "response_covariance": response_cov,
            "identity_residual": output_cov - centroid_cov - response_cov}
