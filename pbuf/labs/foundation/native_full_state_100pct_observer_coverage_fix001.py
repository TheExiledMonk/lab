#!/usr/bin/env python3
"""Launch-aware repair for PBUF full-state 100% observer coverage 001.

PR #106 correctly changed source-plane coverage from 25% to 100%, but the reused
PR #104 receipt helper rebuilt its initial 3D ray positions from the historical
25% launcher.  The 100% lane therefore attempted to combine 285156 terminal rays
with 71289 initial rays and failed before any science comparison.

This repair changes no source physics, native response, M10 field, LOS field,
G3D propagation, decoder inventory, observational use, weighting, fitting, or
rescaling.  It changes only the observer-side receipt bookkeeping so that the
actual x0/y0 launch coordinates for each lane are carried into the 3D receipt
extraction.
"""
from __future__ import annotations

import numpy as np

import pbuf.labs.foundation.native_full_state_100pct_observer_coverage001 as LAB
import pbuf.labs.foundation.native_full_received_state_information_retention001 as RET
import pbuf.labs.foundation.native_observable_extraction_method_sweep001 as EX
import pbuf.labs.foundation.native_full_state_2d_reconstruction_decoder_sweep001 as DEC
import pbuf.labs.foundation.native_accumulated_full_lensing001 as G3D


def _binned_received_3d_with_launch(
    screen: dict,
    snap: dict,
    extent: float,
    bins: int,
    x0: np.ndarray,
    y0: np.ndarray,
) -> dict[str, np.ndarray]:
    e1 = RET._finite(screen["e1"])
    e2 = RET._finite(screen["e2"])
    n = RET._finite(screen["normal"])
    u0 = RET._finite(screen["u0"])
    v0 = RET._finite(screen["v0"])
    uf = RET._finite(screen["uf"])
    vf = RET._finite(screen["vf"])

    x0 = RET._finite(x0)
    y0 = RET._finite(y0)
    p0 = np.column_stack((x0, y0, np.zeros_like(x0))).astype(np.float64)
    pf = np.column_stack((snap["x"], snap["y"], snap["z"])).astype(np.float64)
    vel = np.column_stack((snap["vx"], snap["vy"], snap["vz"])).astype(np.float64)

    expected = x0.size
    lengths = {
        "x0": x0.size,
        "y0": y0.size,
        "u0": u0.size,
        "v0": v0.size,
        "uf": uf.size,
        "vf": vf.size,
        "pf": pf.shape[0],
        "vel": vel.shape[0],
    }
    if any(v != expected for v in lengths.values()):
        raise RuntimeError(f"launch/receipt ray-count mismatch: {lengths}")

    w0 = p0 @ n
    wf = pf @ n
    du = uf - u0
    dv = vf - v0
    dw = wf - w0
    t1 = vel @ e1
    t2 = vel @ e2
    tn = vel @ n

    scalars = {"du": du, "dv": dv, "dw": dw, "t1": t1, "t2": t2, "tn": tn}
    out = {name: RET._empty(bins) for name in RET.PRIMARY_3D_BIN_CHANNELS}
    r, c, valid = RET._bin_indices(u0, v0, extent, bins)
    flat = r * bins + c

    for q in np.unique(flat[valid]):
        idx = np.where(valid & (flat == q))[0]
        rr, cc = divmod(int(q), bins)
        if idx.size == 0:
            continue

        for key in ("du", "dv", "dw", "t1", "t2", "tn"):
            vals = scalars[key][idx]
            out[f"mean_{key}"][rr, cc] = float(np.mean(vals))
            out[f"std_{key}"][rr, cc] = float(np.std(vals))

        if idx.size >= 2:
            D = np.column_stack((du[idx], dv[idx], dw[idx]))
            C = np.cov(D, rowvar=False, ddof=1)
            out["cov_du_dv"][rr, cc] = float(C[0, 1])
            out["cov_du_dw"][rr, cc] = float(C[0, 2])
            out["cov_dv_dw"][rr, cc] = float(C[1, 2])

        if idx.size >= 6:
            X = np.column_stack((u0[idx] - np.mean(u0[idx]), v0[idx] - np.mean(v0[idx])))
            Y = np.column_stack((uf[idx] - np.mean(uf[idx]), vf[idx] - np.mean(vf[idx]), wf[idx] - np.mean(wf[idx])))
            try:
                A, *_ = np.linalg.lstsq(X, Y, rcond=None)
            except np.linalg.LinAlgError:
                continue
            out["j3_e1_u"][rr, cc] = float(A[0, 0])
            out["j3_e1_v"][rr, cc] = float(A[1, 0])
            out["j3_e2_u"][rr, cc] = float(A[0, 1])
            out["j3_e2_v"][rr, cc] = float(A[1, 1])
            out["j3_n_u"][rr, cc] = float(A[0, 2])
            out["j3_n_v"][rr, cc] = float(A[1, 2])

    out["_per_ray_full3d"] = np.column_stack((du, dv, dw, t1, t2, tn))
    out["_per_ray_transverse2d"] = np.column_stack((du, dv, t1, t2))
    return out


def _decode_lane_launch_aware(
    data: dict,
    channel: dict,
    chain: dict,
    x0: np.ndarray,
    y0: np.ndarray,
    label: str,
) -> dict:
    snap = chain["checkpoints"][G3D.CHECKPOINT]
    screen = EX._screen_coordinates(x0, y0, snap)
    extracted = EX._extract_all(screen, LAB.EXTENT, LAB.BINS)
    receipt3d = _binned_received_3d_with_launch(
        screen, snap, LAB.EXTENT, LAB.BINS, x0, y0
    )
    bank, family = RET._decoded_bank(extracted, receipt3d)

    candidates, decoder_meta = DEC._build_candidates(bank, family)
    per_ray = RET._per_ray_geometry(receipt3d)
    stages = RET._stage_metrics(bank, family)

    return {
        "label": label,
        "data": data,
        "channel": channel,
        "chain": chain,
        "screen": screen,
        "bank": bank,
        "family": family,
        "candidates": candidates,
        "decoder_meta": decoder_meta,
        "per_ray_information_geometry": per_ray,
        "stage_information_geometry": stages,
        "support": LAB._source_support(x0, y0),
    }


def main() -> int:
    LAB._decode_lane = _decode_lane_launch_aware
    return LAB.main()


if __name__ == "__main__":
    raise SystemExit(main())
