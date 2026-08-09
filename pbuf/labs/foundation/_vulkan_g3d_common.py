"""Shared, target-blind preparation and reporting for Dev Doc 109 labs."""

from __future__ import annotations

import hashlib
import numpy as np

from pbuf.labs.foundation import native_full_state_2d_reconstruction_decoder_sweep001 as DEC
from pbuf.wl.channels import decode_full_channel_bank
from pbuf.wl.interface import get_interface_vector
from pbuf.wl.launch import launch_25pct, launch_100pct
from pbuf.wl.los import project_interface_to_los
from pbuf.wl.native_response import build_native_response
from pbuf.wl.received_state import build_received_state
from pbuf.wl.reconstruction import build_reconstruction_candidates
from pbuf.wl.screen import build_detector_screen
from pbuf.wl.source import load_cluster_source

RAY_RTOL, RAY_ATOL = 1e-10, 1e-12
SCREEN_RTOL, SCREEN_ATOL = 1e-10, 1e-12
BANK_RTOL, BANK_ATOL = 1e-9, 1e-11


def prepare(cluster: dict, coverage: str) -> dict:
    source = load_cluster_source(cluster)
    native = build_native_response(source["rho3"])
    los = project_interface_to_los(get_interface_vector(native))
    launch = launch_25pct() if coverage == "25pct" else launch_100pct()
    return {"cluster": cluster, "source": source, "native": native, "los": los, "launch": launch}


def downstream(prepared: dict, propagation: dict) -> dict:
    launch = prepared["launch"]
    screen = build_detector_screen(launch, propagation)
    received = build_received_state(launch, propagation, screen)
    decoded = decode_full_channel_bank(screen, received)
    candidates, meta = build_reconstruction_candidates(decoded["bank"], decoded["family"])
    metrics = DEC._compare_candidates(candidates, DEC._targets_after_decoding(prepared["source"]["data"]))
    return {"screen": screen, "received": received, "bank": decoded["bank"],
            "family": decoded["family"], "candidates": candidates, "meta": meta, "metrics": metrics}


def sha256_f64(value) -> str:
    return hashlib.sha256(np.ascontiguousarray(value, dtype=np.float64).tobytes()).hexdigest()


def report(a, b, rtol, atol) -> dict:
    x, y = np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64)
    finite = np.isfinite(x) & np.isfinite(y)
    delta = x[finite] - y[finite]
    denom = max(float(np.sqrt(np.mean(x[finite]**2))) if np.any(finite) else 0.0, 1e-30)
    pearson = float("nan")
    if np.count_nonzero(finite) >= 2 and np.std(x[finite]) and np.std(y[finite]):
        pearson = float(np.corrcoef(x[finite], y[finite])[0, 1])
    return {"pass": bool(np.allclose(x, y, rtol=rtol, atol=atol, equal_nan=True)),
            "max_abs_error": float(np.max(np.abs(delta))) if delta.size else 0.0,
            "relative_rms_error": float(np.sqrt(np.mean(delta**2))/denom) if delta.size else 0.0,
            "pearson": pearson}


def numeric_dict_report(a: dict, b: dict, rtol: float, atol: float) -> dict:
    return {name: report(a[name], b[name], rtol, atol) for name in a
            if name in b and np.issubdtype(np.asarray(a[name]).dtype, np.number)}


def metric_max_error(a: dict, b: dict) -> float:
    errors = []
    for candidate in a:
        for target in a[candidate]:
            for metric in ("pearson", "spearman"):
                x, y = a[candidate][target][metric], b[candidate][target][metric]
                if np.isnan(x) and np.isnan(y): continue
                errors.append(abs(float(x)-float(y)))
    return max(errors, default=0.0)
