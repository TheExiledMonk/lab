"""Target-blind transport/receiver relations and diagnostic reconstruction."""
from __future__ import annotations

import numpy as np


def effective_rank(x, tol=None):
    x = np.asarray(x, float)
    if x.size == 0 or x.shape[1] == 0:
        return 0
    x = np.nan_to_num(x-x.mean(0)); s = np.linalg.svd(x, compute_uv=False)
    threshold = (max(x.shape)*np.finfo(float).eps*s[0]) if tol is None and len(s) else (tol or 0)
    return int(np.sum(s > threshold))


def spin2_alignment(theta1, theta2):
    return np.cos(2.0*(np.asarray(theta1)-np.asarray(theta2)))


def relation_bank(transport, first, receiver_bank):
    """Restricted relational families; empty if no receiver field is available."""
    if not receiver_bank:
        return {}
    du, dv = transport["delta_u"].ravel(), transport["delta_v"].ravel()
    theta_t = np.arctan2(dv, du)
    out = {}
    prefixes = sorted({k.rsplit("_grad_", 1)[0] for k in receiver_bank if "_grad_" in k})
    for p in prefixes:
        gu, gv = receiver_bank.get(p+"_grad_0"), receiver_bank.get(p+"_grad_1")
        if gu is None or gv is None:
            continue
        gnorm = np.hypot(gu, gv)+np.finfo(float).eps
        out[p+"_transport_parallel"] = (du*gu+dv*gv)/gnorm
        out[p+"_transport_perpendicular"] = (-du*gv+dv*gu)/gnorm
        out[p+"_transport_spin2_alignment"] = spin2_alignment(theta_t, np.arctan2(gv, gu))
    return out


def neutral_textures(shape):
    """Deterministic uniform/grid, dots, and four fixed-angle bar probes."""
    y, x = np.indices(shape); period = max(8, min(shape)//16)
    tex = {"uniform": np.ones(shape),
           "checker_grid": (((x//period)+(y//period)) % 2).astype(float),
           "isotropic_dots": (((x % period)==period//2)&((y % period)==period//2)).astype(float)}
    for angle, key in ((0, "bars_0"), (45, "bars_45"), (90, "bars_90"), (135, "bars_135")):
        t = np.deg2rad(angle); coord=x*np.cos(t)+y*np.sin(t)
        tex[key] = ((np.mod(coord, period)) < max(1, period//4)).astype(float)
    return tex


def rasterize(u, v, weights, resolution=64, bounds=(-8., 8.)):
    lo, hi = bounds; u=np.asarray(u).ravel();v=np.asarray(v).ravel();w=np.asarray(weights).ravel()
    ix=np.floor((u-lo)*resolution/(hi-lo)).astype(int);iy=np.floor((v-lo)*resolution/(hi-lo)).astype(int)
    ok=np.isfinite(u)&np.isfinite(v)&np.isfinite(w)&(ix>=0)&(iy>=0)&(ix<resolution)&(iy<resolution)
    cell=iy[ok]*resolution+ix[ok]
    sums=np.bincount(cell,weights=w[ok],minlength=resolution**2)
    count=np.bincount(cell,minlength=resolution**2)
    return np.divide(sums,count,out=np.zeros_like(sums),where=count>0).reshape(resolution,resolution)


def reconstruct_neutral(transport, resolution=64):
    shape=transport["u0"].shape; textures=neutral_textures(shape); out={}
    for name, image in textures.items():
        out["before_"+name]=rasterize(transport["u0"],transport["v0"],image,resolution)
        out["transport_only_"+name]=rasterize(transport["uf"],transport["vf"],image,resolution)
    return out
