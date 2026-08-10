"""Correspondence-preserving launch-to-receipt transport fields (Dev 121).

Only frozen ray checkpoints are accepted.  All differences are evaluated on the
exact regular launch topology; no received-cell neighbourhood is inferred.
"""
from __future__ import annotations

import numpy as np

from .bundle_transport_geometry import reconstruct_launch_topology, received_state


RAW_KEYS = ("u0", "v0", "uf", "vf", "dx", "dy", "dz", "rx", "ry", "rz")


def build_dual_transport(rays, w_reference=None):
    """Preserve raw T_i and return scalar/vector fields in launch coordinates."""
    topology = reconstruct_launch_topology(rays)
    shape = topology["shape"]
    state = received_state(rays)
    u0 = np.asarray(rays["u0"], float).reshape(shape)
    v0 = np.asarray(rays["v0"], float).reshape(shape)
    uf = state["received_u"].reshape(shape)
    vf = state["received_v"].reshape(shape)
    wf = state["received_w"].reshape(shape)
    if w_reference is None:
        # A constant origin shift contains no transport geometry.  The median is
        # deterministic and robust, and the unshifted received_w is retained.
        w_reference = float(np.median(wf[np.isfinite(wf)]))
    out = {
        "ray_id": topology["ray_id"], "u0": u0, "v0": v0,
        "uf": uf, "vf": vf, "wf": wf,
        "delta_u": uf-u0, "delta_v": vf-v0,
        "delta_w": wf-w_reference,
        "dir_u": state["received_dir_u"].reshape(shape),
        "dir_v": state["received_dir_v"].reshape(shape),
        "dir_w": state["received_dir_w"].reshape(shape),
        "w_reference": np.asarray(w_reference),
        "launch_spacing_u": np.asarray(topology["spacing"][0]),
        "launch_spacing_v": np.asarray(topology["spacing"][1]),
    }
    # Preserve every native raw ray quantity, including world coordinates.
    for key in RAW_KEYS:
        out["raw_"+key] = np.asarray(rays[key]).reshape(shape)
    return out


def _gradient(a, du, dv):
    # axis 1 is u and axis 0 is v in the verified row-major launch topology.
    dv_a, du_a = np.gradient(np.asarray(a, float), dv, du, edge_order=2)
    return du_a, dv_a


def first_order_transport(transport):
    """First derivatives of displacement, depth, and received direction."""
    du = float(transport["launch_spacing_u"])
    dv = float(transport["launch_spacing_v"])
    out = {"ray_id": transport["ray_id"]}
    for name in ("delta_u", "delta_v", "wf", "dir_u", "dir_v", "dir_w"):
        gu, gv = _gradient(transport[name], du, dv)
        out[f"d_u_{name}"] = gu
        out[f"d_v_{name}"] = gv
    return out


def second_order_transport(transport, first=None):
    """uu, uv, and vv derivatives for displacement and received depth."""
    first = first if first is not None else first_order_transport(transport)
    du = float(transport["launch_spacing_u"])
    dv = float(transport["launch_spacing_v"])
    out = {"ray_id": transport["ray_id"]}
    for name in ("delta_u", "delta_v", "wf"):
        duu, duv_a = _gradient(first[f"d_u_{name}"], du, dv)
        dvu_a, dvv = _gradient(first[f"d_v_{name}"], du, dv)
        out[f"d_uu_{name}"] = duu
        out[f"d_uv_{name}"] = .5*(duv_a+dvu_a)
        out[f"d_vv_{name}"] = dvv
    return out


def feature_matrix(bank, names=None):
    names = names or [k for k in bank if k != "ray_id" and np.asarray(bank[k]).ndim == 2]
    return np.column_stack([np.asarray(bank[k]).ravel() for k in names])
