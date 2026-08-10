"""Target-blind multiscale relations for the frozen Dev121 R3 transport bank.

The routines in this module contain no lens-cohort or shear access.  They operate
only on scalar fields and the immutable semantic manifest returned below.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import ndimage, signal

SCALES = (2, 4, 8, 16, 32)
EPS = np.finfo(np.float64).eps


def canonical_manifest():
    """Return the exact column order used for Dev121's rank-29 R3 matrix."""
    rows = []
    def add(name, source, order, component, basis, group):
        rows.append({"index": len(rows), "name": name, "source_field": source,
                     "derivative_order": order, "component_type": component,
                     "coordinate_basis": basis, "group": group})
    for name, component in (("u0", "launch_u"), ("v0", "launch_v"),
                            ("uf", "received_u"), ("vf", "received_v"),
                            ("wf", "depth"), ("dir_u", "direction_u"),
                            ("dir_v", "direction_v"), ("dir_w", "direction_w")):
        add(name, "dual_transport", 0, component, "launch_grid", "G0_raw_transport")
    for name in ("d_u_delta_u", "d_v_delta_u", "d_u_delta_v", "d_v_delta_v"):
        add(name, "first_order_transport", 1, "transverse_derivative", "launch_uv", "G1_first_order_transverse")
    for name in ("d_u_wf", "d_v_wf"):
        add(name, "first_order_transport", 1, "depth_derivative", "launch_uv", "G2_first_order_depth")
    for name in ("d_u_dir_u", "d_v_dir_u", "d_u_dir_v", "d_v_dir_v", "d_u_dir_w", "d_v_dir_w"):
        add(name, "first_order_transport", 1, "direction_derivative", "launch_uv", "G3_first_order_direction")
    for target in ("delta_u", "delta_v"):
        for axis in ("uu", "uv", "vv"):
            add(f"d_{axis}_{target}", "second_order_transport", 2, "transverse_derivative", "launch_uv", "G4_second_order_transverse")
    for axis in ("uu", "uv", "vv"):
        add(f"d_{axis}_wf", "second_order_transport", 2, "depth_derivative", "launch_uv", "G5_second_order_depth")
    assert len(rows) == 29
    return rows


def canonical_json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def sha256_json(value):
    return hashlib.sha256((canonical_json(value) + "\n").encode()).hexdigest()


def load_second_order_bank(cluster_dir):
    """Load and validate the frozen 29-channel Dev121 representation."""
    cluster_dir = Path(cluster_dir)
    files = {}
    for stem in ("dual_transport", "first_order_transport", "second_order_transport"):
        with np.load(cluster_dir / f"{stem}.npz", allow_pickle=False) as z:
            files[stem] = {k: np.asarray(z[k]) for k in z.files}
    fields = {}
    for row in canonical_manifest():
        key, source = row["name"], row["source_field"]
        if key not in files[source]:
            raise ValueError(f"frozen Dev121 channel missing: {source}:{key}")
        a = np.asarray(files[source][key], dtype=np.float64)
        if a.ndim != 2:
            raise ValueError(f"channel {key} is not a launch-grid scalar field")
        fields[key] = a
    if len({a.shape for a in fields.values()}) != 1:
        raise ValueError("Dev121 channel shapes differ")
    return fields, files["dual_transport"], canonical_manifest()


def group_indices(manifest=None):
    manifest = canonical_manifest() if manifest is None else manifest
    return {g: [r["index"] for r in manifest if r["group"] == g]
            for g in dict.fromkeys(r["group"] for r in manifest)}


def ablation_lanes(manifest=None):
    manifest = canonical_manifest() if manifest is None else manifest
    groups = np.array([r["group"] for r in manifest]); order = np.array([r["derivative_order"] for r in manifest])
    component = np.array([r["component_type"] for r in manifest]); all_i = np.arange(len(manifest))
    lanes = {"ALL29": all_i, "NO_RAW": all_i[order != 0], "NO_FIRST_ORDER": all_i[order != 1],
             "NO_SECOND_ORDER": all_i[order != 2], "RAW_ONLY": all_i[order == 0],
             "FIRST_ONLY": all_i[order == 1], "SECOND_ONLY": all_i[order == 2],
             "TRANSVERSE_ONLY": all_i[np.char.find(component.astype(str), "transverse") >= 0],
             "DEPTH_ONLY": all_i[np.char.find(component.astype(str), "depth") >= 0],
             "DIRECTION_ONLY": all_i[np.char.find(component.astype(str), "direction") >= 0]}
    return {k: v.tolist() for k, v in lanes.items() if len(v)}


def support_kernel(radius, mode="square"):
    radius = int(radius); y, x = np.mgrid[-radius:radius+1, -radius:radius+1]
    if mode == "square": mask = np.ones_like(x, bool)
    elif mode == "radial": mask = x*x + y*y <= radius*radius
    else: raise ValueError("mode must be square or radial")
    return mask.astype(np.float64) / mask.sum()


def local_mean_variance(field, radius, mode="square"):
    x = np.asarray(field, dtype=np.float64)
    if mode == "square":
        size=2*int(radius)+1
        mean=ndimage.uniform_filter(x,size=size,mode="constant",cval=0.0)
        mean2=ndimage.uniform_filter(x*x,size=size,mode="constant",cval=0.0)
    else:
        k = support_kernel(radius, mode)
        mean = signal.fftconvolve(x, k, mode="same")
        mean2 = signal.fftconvolve(x*x, k, mode="same")
    valid = np.zeros(x.shape, bool); valid[radius:x.shape[0]-radius, radius:x.shape[1]-radius] = True
    return mean, np.maximum(mean2 - mean*mean, 0.0), valid


def derivatives(field, spacing=(1.0, 1.0)):
    """Float64 gradient and Hessian invariants (axis 1 is u, axis 0 is v)."""
    x = np.asarray(field, dtype=np.float64); dv, du = np.gradient(x, spacing[1], spacing[0], edge_order=2)
    duu = np.gradient(du, spacing[0], axis=1, edge_order=2)
    dvv = np.gradient(dv, spacing[1], axis=0, edge_order=2)
    duv = .5*(np.gradient(du, spacing[1], axis=0, edge_order=2)+np.gradient(dv, spacing[0], axis=1, edge_order=2))
    diff = duu-dvv
    return {"gradient_u":du, "gradient_v":dv, "gradient_magnitude":np.hypot(du,dv),
            "gradient_orientation":np.arctan2(dv,du), "trace":duu+dvv,
            "determinant":duu*dvv-duv*duv, "eigenvalue_difference":np.hypot(diff,2*duv),
            "principal_axis":.5*np.arctan2(2*duv,diff), "hessian_q1":diff, "hessian_q2":2*duv}


def spatial_quadrupole(field, radius, variant="signed", mode="square"):
    """Spatial quadrupole using signed residual or positive energy weights."""
    x=np.asarray(field,np.float64); r=int(radius); y,u=np.mgrid[-r:r+1,-r:r+1]
    mask=support_kernel(r,mode)>0; base=np.ones(mask.shape,float)*mask
    if variant == "signed":
        mean=ndimage.convolve(x,base/base.sum(),mode="constant",cval=0.0); w=x-mean
    elif variant == "energy": w=x*x
    else: raise ValueError("variant must be signed or energy")
    if mode == "square":
        den=ndimage.uniform_filter(np.abs(w),size=2*r+1,mode="constant",cval=0.0)*(2*r+1)**2+EPS
    else: den=ndimage.convolve(np.abs(w),base,mode="constant",cval=0.0)+EPS
    # The center-specific mean cancels from both symmetric moment kernels.
    numerator_field = x if variant == "signed" else w
    if mode == "square":
        one=np.ones(2*r+1); coord=np.arange(-r,r+1,dtype=float); square=coord*coord
        sum_u2=ndimage.correlate1d(ndimage.correlate1d(numerator_field,square,axis=1,mode="constant",cval=0),one,axis=0,mode="constant",cval=0)
        sum_v2=ndimage.correlate1d(ndimage.correlate1d(numerator_field,one,axis=1,mode="constant",cval=0),square,axis=0,mode="constant",cval=0)
        cross=ndimage.correlate1d(ndimage.correlate1d(numerator_field,coord,axis=1,mode="constant",cval=0),coord,axis=0,mode="constant",cval=0)
        q1=(sum_u2-sum_v2)/den;q2=2*cross/den
    else:
        q1=ndimage.convolve(numerator_field,base*(u*u-y*y),mode="constant",cval=0.0)/den
        q2=ndimage.convolve(numerator_field,base*(2*u*y),mode="constant",cval=0.0)/den
    return {"q1":q1,"q2":q2,"q_abs":np.hypot(q1,q2),"q_angle":.5*np.arctan2(q2,q1)}


def cross_quadrupole(a, b, radius, mode="square"):
    ma,_,_=local_mean_variance(a,radius,mode); mb,_,_=local_mean_variance(b,radius,mode)
    return spatial_quadrupole((np.asarray(a)-ma)*(np.asarray(b)-mb),radius,"energy" if False else "signed",mode)


def matrix_diagnostics(x):
    x=np.asarray(x,np.float64); x=np.nan_to_num(x-x.mean(0)); scale=x.std(0)
    # Relational banks mix coordinates, derivatives, variances and moments.  A
    # frozen unit-variance normalization prevents physical units from defining
    # the numerical-rank tolerance; constant columns remain exactly zero.
    x=np.divide(x,scale,out=np.zeros_like(x),where=scale>0);s=np.linalg.svd(x,compute_uv=False)
    if not len(s) or s[0] == 0: return {"effective_rank":0,"variance_95_components":0,"variance_99_components":0,"condition_number":0.0}
    tol=max(x.shape)*np.finfo(float).eps*s[0]; rank=int((s>tol).sum()); v=s*s; c=np.cumsum(v)/v.sum()
    return {"effective_rank":rank,"variance_95_components":int(np.searchsorted(c,.95)+1),
            "variance_99_components":int(np.searchsorted(c,.99)+1),
            "condition_number":float(s[0]/s[rank-1]) if rank else 0.0}


def spin2_rotate(q1, q2, phi):
    c,s=np.cos(2*phi),np.sin(2*phi)
    return q1*c-q2*s, q1*s+q2*c


def scale_persistence(q_a, q_b):
    a1,a2=q_a; b1,b2=q_b; aa=np.hypot(a1,a2); bb=np.hypot(b1,b2); ok=(aa>EPS)&(bb>EPS)
    orientation=float(np.mean((a1[ok]*b1[ok]+a2[ok]*b2[ok])/(aa[ok]*bb[ok]))) if ok.any() else 0.0
    mag=float(np.corrcoef(aa.ravel(),bb.ravel())[0,1]) if aa.std() and bb.std() else 0.0
    sign=float(np.mean((a1[ok]*b1[ok]+a2[ok]*b2[ok])>=0)) if ok.any() else 0.0
    return {"orientation_persistence":orientation,"magnitude_correlation":mag,"sign_consistency":sign}
