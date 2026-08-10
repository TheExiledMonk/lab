"""Launch-aware received 3D state, with no hidden launch reconstruction."""

import numpy as np

from pbuf.labs.foundation import native_full_received_state_information_retention001 as RET
from .config import EXTENT, OBS_BINS
from .launch import RayLaunch


def build_received_state(launch: RayLaunch, propagation: dict, screen: dict) -> dict:
    snap = propagation["final_snapshot"]
    arrays = {
        "launch.x0": launch.x0, "launch.y0": launch.y0,
        "screen.u0": screen["u0"], "screen.v0": screen["v0"],
        "screen.uf": screen["uf"], "screen.vf": screen["vf"],
        **{f"final_snapshot.{k}": snap[k] for k in ("x", "y", "z", "vx", "vy", "vz")},
    }
    lengths = {name: int(np.asarray(value).shape[0]) for name, value in arrays.items()}
    if len(set(lengths.values())) != 1:
        raise ValueError(f"received-state ray-count mismatch: {lengths}")

    e1, e2, n = (RET._finite(screen[k]) for k in ("e1", "e2", "normal"))
    x0, y0 = RET._finite(launch.x0), RET._finite(launch.y0)
    u0, v0, uf, vf = (RET._finite(screen[k]) for k in ("u0", "v0", "uf", "vf"))
    p0 = np.column_stack((x0, y0, np.zeros_like(x0))).astype(np.float64)
    pf = np.column_stack((snap["x"], snap["y"], snap["z"])).astype(np.float64)
    vel = np.column_stack((snap["vx"], snap["vy"], snap["vz"])).astype(np.float64)
    du, dv = uf - u0, vf - v0
    scalars = {
        "du": du, "dv": dv, "dw": pf @ n - p0 @ n,
        "t1": vel @ e1, "t2": vel @ e2, "tn": vel @ n,
    }
    out = {name: RET._empty(OBS_BINS) for name in RET.PRIMARY_3D_BIN_CHANNELS}
    r, c, valid = RET._bin_indices(u0, v0, EXTENT, OBS_BINS)
    flat = r * OBS_BINS + c
    for q in np.unique(flat[valid]):
        idx = np.where(valid & (flat == q))[0]
        rr, cc = divmod(int(q), OBS_BINS)
        for key, values in scalars.items():
            out[f"mean_{key}"][rr, cc] = float(np.mean(values[idx]))
            out[f"std_{key}"][rr, cc] = float(np.std(values[idx]))
        if idx.size >= 2:
            cov = np.cov(np.column_stack((du[idx], dv[idx], scalars["dw"][idx])), rowvar=False, ddof=1)
            out["cov_du_dv"][rr, cc], out["cov_du_dw"][rr, cc], out["cov_dv_dw"][rr, cc] = cov[0, 1], cov[0, 2], cov[1, 2]
        if idx.size >= 6:
            X = np.column_stack((u0[idx] - np.mean(u0[idx]), v0[idx] - np.mean(v0[idx])))
            Y = np.column_stack((uf[idx] - np.mean(uf[idx]), vf[idx] - np.mean(vf[idx]), (pf @ n)[idx] - np.mean((pf @ n)[idx])))
            try:
                A, *_ = np.linalg.lstsq(X, Y, rcond=None)
            except np.linalg.LinAlgError:
                continue
            for key, value in zip(("j3_e1_u", "j3_e1_v", "j3_e2_u", "j3_e2_v", "j3_n_u", "j3_n_v"),
                                  (A[0, 0], A[1, 0], A[0, 1], A[1, 1], A[0, 2], A[1, 2])):
                out[key][rr, cc] = float(value)
    out["_per_ray_full3d"] = np.column_stack(tuple(scalars.values()))
    out["_per_ray_transverse2d"] = np.column_stack((du, dv, scalars["t1"], scalars["t2"]))
    return out
