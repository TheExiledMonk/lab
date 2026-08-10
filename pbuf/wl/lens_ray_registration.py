"""Target-blind lens-footprint registration and received-ray redistribution tools."""
from __future__ import annotations

import hashlib
import numpy as np
from scipy import ndimage

COHORT_NAMES = ("RAY_UNCLASSIFIED", "RAY_LENS_CORE", "RAY_LENS_INNER",
                "RAY_LENS_OUTER", "RAY_NEAR_LENS", "RAY_FAR_CONTROL")


def fingerprint(a) -> str:
    return hashlib.sha256(np.ascontiguousarray(a, dtype=np.float64).tobytes()).hexdigest()


def register_source_geometry_to_propagation(source_field, extent):
    """Derive source voxel-centre coordinates in propagation coordinates.

    The canonical proxy and launch screen share the configured square extent;
    no observed image or fitted transform enters this mapping.
    """
    a = np.asarray(source_field, float)
    ny, nx = a.shape
    x = np.linspace(-extent + extent/nx, extent - extent/nx, nx)
    y = np.linspace(-extent + extent/ny, extent - extent/ny, ny)
    return {"matrix": np.eye(3), "offset": np.zeros(3), "x": x, "y": y,
            "rule": "shared canonical source/proxy and propagation extent; voxel centres"}


def register_propagation_geometry_to_observer(e1, e2):
    e1, e2 = np.asarray(e1, float), np.asarray(e2, float)
    n = np.cross(e1, e2); n /= np.linalg.norm(n)
    return {"matrix": np.vstack((e1, e2, n)), "offset": np.zeros(3),
            "rule": "frozen detector orthonormal basis"}


def build_frozen_lens_masks(source_field):
    """Build exclusive deterministic source masks before ray inspection."""
    a = np.asarray(source_field, float)
    positive = a[a > 0]
    if not positive.size:
        raise ValueError("empty source loading")
    p70, p90 = np.percentile(positive, (70, 90))
    explicit = a > 0
    # A resampled positive proxy can have a low nonzero numerical floor over
    # essentially the entire field; that is not an explicit support boundary.
    # Freeze a source-only lower support percentile in that case.
    if explicit.mean() > .8:
        p30 = np.percentile(positive, 30)
        support = a >= p30
        support_rule = ">=p30 positive source loading (no explicit boundary)"
    else:
        p30 = None; support = explicit; support_rule = ">0 explicit support"
    core = a >= p90
    inner_total = a >= p70
    outer_total = support
    near_total = ndimage.binary_dilation(support, iterations=1)
    near = near_total & ~support
    p10 = np.percentile(positive, 10)
    far = (a <= p10) & ~near_total
    if not np.any(far): far = a <= p10
    far_rule = "lowest 10% source-loading region outside the near shell (source-only control)"
    return {"L0_CORE": core, "L1_INNER": inner_total & ~core,
            "L2_OUTER": outer_total & ~inner_total, "L3_NEAR_FIELD": near,
            "L4_FAR_CONTROL": far}, {"rule": f"{support_rule}; >=p90 core; >=p70 inner; one-cell near shell; {far_rule}",
                                      "positive_p10": None if p10 is None else float(p10), "positive_p30": None if p30 is None else float(p30), "positive_p70": float(p70), "positive_p90": float(p90)}


def classify_launch_rays(launch_u, launch_v, masks, extent):
    shape = next(iter(masks.values())).shape; ny, nx = shape
    col = np.floor((np.asarray(launch_u)+extent)/(2*extent)*nx).astype(int)
    row = np.floor((np.asarray(launch_v)+extent)/(2*extent)*ny).astype(int)
    valid = np.isfinite(launch_u) & np.isfinite(launch_v) & (row>=0)&(row<ny)&(col>=0)&(col<nx)
    ids = np.zeros(len(col), np.uint8)
    for cid, name in enumerate(("L0_CORE","L1_INNER","L2_OUTER","L3_NEAR_FIELD","L4_FAR_CONTROL"), 1):
        hit = np.zeros(len(col), bool); hit[valid] = masks[name][row[valid], col[valid]]
        ids[(ids == 0) & hit] = cid
    return ids


def observer_histogram(u, v, bins, extent, weights=None):
    h, _, _ = np.histogram2d(u, v, bins=bins, range=[[-extent,extent],[-extent,extent]], weights=weights)
    return h


def trace_received_positions(rays, cohort_id):
    e1,e2=np.asarray(rays["e1"],float),np.asarray(rays["e2"],float)
    n=np.cross(e1,e2); n/=np.linalg.norm(n)
    xyz=np.column_stack((rays["rx"],rays["ry"],rays["rz"]))
    w=xyz@n
    return {"ray_id":np.arange(len(cohort_id),dtype=np.int64), "cohort_id":cohort_id,
            "launch_u":np.asarray(rays["u0"]), "launch_v":np.asarray(rays["v0"]),
            "received_u":np.asarray(rays["uf"]), "received_v":np.asarray(rays["vf"]), "received_w":w,
            "delta_u":np.asarray(rays["uf"])-rays["u0"], "delta_v":np.asarray(rays["vf"])-rays["v0"],
            "delta_w":w-(np.column_stack((rays["launch_x"],rays["launch_y"],np.zeros(len(w))))@n),
            "dir_u":np.asarray(rays["dx"]), "dir_v":np.asarray(rays["dy"]), "dir_w":np.asarray(rays["dz"]),
            "expected_received_u":np.asarray(rays["u0"]), "expected_received_v":np.asarray(rays["v0"])}


def quantiles(x):
    x=np.asarray(x,float); x=x[np.isfinite(x)]
    return {f"q{q:02d}":float(np.percentile(x,q)) for q in (5,25,50,75,95)} if x.size else {f"q{q:02d}":None for q in (5,25,50,75,95)}


def cohort_statistics(trace):
    out={}; dr=np.hypot(trace["delta_u"],trace["delta_v"])
    for cid,name in enumerate(COHORT_NAMES):
        m=trace["cohort_id"]==cid; n=int(m.sum())
        def med(k): return float(np.median(trace[k][m])) if n else None
        def var(k): return float(np.var(trace[k][m])) if n else None
        out[name]={"ray_count":n,"median_delta_u":med("delta_u"),"median_delta_v":med("delta_v"),
          "median_transverse_displacement":float(np.median(dr[m])) if n else None,
          "rms_transverse_displacement":float(np.sqrt(np.mean(dr[m]**2))) if n else None,
          "received_u_variance":var("received_u"),"received_v_variance":var("received_v"),
          "received_depth_variance":var("received_w"),
          "received_direction_variance":float(np.mean([var(k) or 0 for k in ("dir_u","dir_v","dir_w")])) if n else None,
          "transverse_displacement_quantiles":quantiles(dr[m])}
    return out


def separability(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float); a=a[np.all(np.isfinite(a),1)]; b=b[np.all(np.isfinite(b),1)]
    if min(len(a),len(b))<2: return {"centroid_distance":None,"within_between_ratio":None,"mahalanobis_distance":None}
    ma,mb=a.mean(0),b.mean(0); d=ma-mb; cov=(np.cov(a,rowvar=False)+np.cov(b,rowvar=False))/2
    return {"centroid_distance":float(np.linalg.norm(d)),
            "within_between_ratio":float((np.mean(np.var(a,0))+np.mean(np.var(b,0)))/(np.mean(d*d)+1e-30)),
            "mahalanobis_distance":float(np.sqrt(max(0,d@np.linalg.pinv(np.atleast_2d(cov))@d)))}
