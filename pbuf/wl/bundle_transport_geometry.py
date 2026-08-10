"""Target-blind launch-bundle transport geometry (Dev 120).

This module intentionally has no benchmark, lens-mask, or observed-shear imports.
It works from stable per-ray launch/receive identities and exact regular-grid
adjacency.  ``A`` below is a native transport matrix, not a GR lens Jacobian.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import ndimage

SCALES = (1, 2, 4, 8)
PARITY_EPS = 1e-10
DIR_EPS = 1e-15


def structural_sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def reconstruct_launch_topology(rays):
    """Return exact row-major regular topology, refusing inferred kNN topology."""
    n = len(rays["uf"])
    side = int(round(np.sqrt(n)))
    if side * side != n or "launch_x" not in rays or "launch_y" not in rays:
        raise RuntimeError("DEV120_REQUIRED_LAUNCH_RECEIVE_IDENTITY_MISSING")
    x = np.asarray(rays["launch_x"]).reshape(side, side)
    y = np.asarray(rays["launch_y"]).reshape(side, side)
    if not (np.allclose(x, x[0:1], rtol=0, atol=2e-7) and
            np.allclose(y, y[:, 0:1], rtol=0, atol=2e-7)):
        raise RuntimeError("DEV120_REQUIRED_LAUNCH_RECEIVE_IDENTITY_MISSING")
    dx = np.diff(x[0]); dy = np.diff(y[:, 0])
    if np.any(dx <= 0) or np.any(dy <= 0) or not (
            np.allclose(dx, np.median(dx), rtol=2e-5) and
            np.allclose(dy, np.median(dy), rtol=2e-5)):
        raise RuntimeError("DEV120_REQUIRED_LAUNCH_RECEIVE_IDENTITY_MISSING")
    return {"shape": (side, side), "x": x, "y": y,
            "spacing": (float(np.median(dx)), float(np.median(dy))),
            "ray_id": np.arange(n, dtype=np.int64).reshape(side, side)}


def received_state(rays):
    e1, e2 = np.asarray(rays["e1"]), np.asarray(rays["e2"])
    ew = np.cross(e1, e2)
    pos = np.column_stack((rays["rx"], rays["ry"], rays["rz"]))
    direction = np.column_stack((rays["dx"], rays["dy"], rays["dz"]))
    return {
        "received_u": np.asarray(rays["uf"], float),
        "received_v": np.asarray(rays["vf"], float),
        "received_w": pos @ ew,
        "received_dir_u": direction @ e1,
        "received_dir_v": direction @ e2,
        "received_dir_w": direction @ ew,
    }


def _box_sum(a, radius):
    size = 2 * radius + 1
    return ndimage.uniform_filter(np.asarray(a, float), size=size,
                                  mode="constant", cval=0.) * (size * size)


def _local_affine(x, y, outputs, radius):
    """Vectorized deterministic least squares on square launch neighborhoods."""
    one = np.ones_like(x)
    n = _box_sum(one, radius)
    sx, sy = _box_sum(x, radius), _box_sum(y, radius)
    mx, my = sx / n, sy / n
    cxx = _box_sum(x*x, radius) / n - mx*mx
    cyy = _box_sum(y*y, radius) / n - my*my
    cxy = _box_sum(x*y, radius) / n - mx*my
    det = cxx*cyy-cxy*cxy
    cond = np.full_like(x, np.inf)
    tr = cxx+cyy; disc = np.sqrt(np.maximum((cxx-cyy)**2+4*cxy*cxy, 0))
    lo, hi = .5*(tr-disc), .5*(tr+disc)
    np.divide(hi, lo, out=cond, where=lo>0)
    mats = np.zeros(x.shape + (len(outputs), 2))
    means = []
    residual_var = np.zeros_like(x)
    for j, z in enumerate(outputs):
        mz = _box_sum(z, radius)/n; means.append(mz)
        cxz = _box_sum(x*z, radius)/n-mx*mz
        cyz = _box_sum(y*z, radius)/n-my*mz
        mats[..., j, 0] = np.divide(cyy*cxz-cxy*cyz, det,
                                    out=np.zeros_like(x), where=np.abs(det)>1e-30)
        mats[..., j, 1] = np.divide(cxx*cyz-cxy*cxz, det,
                                    out=np.zeros_like(x), where=np.abs(det)>1e-30)
        explained = mats[..., j, 0]*cxz + mats[..., j, 1]*cyz
        residual_var += np.maximum(_box_sum(z*z, radius)/n-mz*mz-explained, 0)
    pred_center = np.stack(means, -1) + np.einsum(
        "...ij,...j->...i", mats, np.stack((x-mx, y-my), -1))
    center_resid = np.linalg.norm(np.stack(outputs, -1)-pred_center, axis=-1)
    med = ndimage.median_filter(center_resid, size=2*radius+1, mode="nearest")
    rank = (lo > np.finfo(float).eps).astype(np.int8) * 2
    valid = (n >= 6) & (rank == 2) & (cond <= 1e10) & np.all(np.isfinite(mats), axis=(-2, -1))
    return mats, n.astype(np.int32), cond, np.sqrt(residual_var), med, rank, valid


def _polar_fields(A):
    # A=R U through batched SVD.  Singular values are principal stretches.
    u, s, vh = np.linalg.svd(A)
    R = u @ vh
    U = np.swapaxes(vh, -1, -2) @ (s[..., :, None] * vh)
    rot = np.arctan2(R[..., 1, 0], R[..., 0, 0])
    q1, q2 = U[..., 0, 0]-U[..., 1, 1], 2*U[..., 0, 1]
    qabs = np.hypot(q1, q2)
    axis = .5*np.arctan2(q2, q1)
    aniso = (s[..., 0]-s[..., 1])/(np.abs(s).sum(-1)+1e-30)
    return R, U, s, rot, q1, q2, qabs, axis, aniso


def _triangle_flips(x, y, uf, vf):
    before = (x[:-1, 1:]-x[:-1, :-1])*(y[1:, :-1]-y[:-1, :-1]) - \
             (y[:-1, 1:]-y[:-1, :-1])*(x[1:, :-1]-x[:-1, :-1])
    after = (uf[:-1, 1:]-uf[:-1, :-1])*(vf[1:, :-1]-vf[:-1, :-1]) - \
            (vf[:-1, 1:]-vf[:-1, :-1])*(uf[1:, :-1]-uf[:-1, :-1])
    out = np.zeros_like(x, dtype=bool)
    out[:-1, :-1] = np.signbit(before) != np.signbit(after)
    return out


def fit_bundle_scale(rays, topology, scale, S_ray, D_ray):
    shape = topology["shape"]
    x = np.asarray(rays["u0"], float).reshape(shape)
    y = np.asarray(rays["v0"], float).reshape(shape)
    state = received_state(rays)
    uf, vf, wf = (state[k].reshape(shape) for k in
                  ("received_u", "received_v", "received_w"))
    B, count, cond, rms3, med, rank, valid = _local_affine(x, y, (uf, vf, wf), scale)
    A = B[..., :2, :]
    _, U, sv, rotation, q1, q2, qabs, qangle, aniso = _polar_fields(A)
    det = np.linalg.det(A)
    parity = np.zeros(shape, np.int8)
    parity[det > PARITY_EPS] = 1; parity[det < -PARITY_EPS] = -1
    theta = np.arctan2(state["received_dir_v"], state["received_dir_u"]).reshape(shape)
    dnorm = np.hypot(state["received_dir_u"], state["received_dir_v"]).reshape(shape)
    dm = (dnorm > DIR_EPS).astype(float)
    den = _box_sum(dm, scale)
    zr = np.divide(_box_sum(np.cos(2*theta)*dm, scale), den,
                   out=np.zeros(shape), where=den>0)
    zi = np.divide(_box_sum(np.sin(2*theta)*dm, scale), den,
                   out=np.zeros(shape), where=den>0)
    # Ordinary vector mean is retained as the explicit control.
    orr = np.divide(_box_sum(np.cos(theta)*dm, scale), den,
                    out=np.zeros(shape), where=den>0)
    ori = np.divide(_box_sum(np.sin(theta)*dm, scale), den,
                    out=np.zeros(shape), where=den>0)
    flips = _triangle_flips(x, y, uf, vf)
    sign_change = ((parity != np.roll(parity, 1, 0)) | (parity != np.roll(parity, 1, 1))) & (parity != 0)
    sign_change[[0], :] = False; sign_change[:, [0]] = False
    near = np.abs(det) <= PARITY_EPS
    fold = sign_change.astype(np.uint8) + flips.astype(np.uint8) + near.astype(np.uint8)
    S = np.asarray(S_ray).reshape(shape); D = np.asarray(D_ray).reshape(shape)
    gs = np.gradient(S); gd = np.gradient(D)
    ts, td = np.arctan2(gs[0], gs[1]), np.arctan2(gd[0], gd[1])
    asq = np.cos(2*(ts-qangle)); adq = np.cos(2*(td-qangle))
    sd_align = np.cos(2*(ts-td))
    gw_u, gw_v = B[..., 2, 0], B[..., 2, 1]
    gw_angle = np.arctan2(gw_v, gw_u)
    return {
        "scale": np.full(shape, scale, np.int8), "ray_id": topology["ray_id"],
        "launch_center_u": x, "launch_center_v": y,
        "received_center_u": uf, "received_center_v": vf, "received_center_w": wf,
        "bundle_matrix_2d": A, "bundle_matrix_3d": B,
        "ray_count": count, "condition_number": cond, "fit_residual": rms3,
        "fit_residual_median": med, "valid_rank": rank, "valid": valid,
        "rotation_angle": rotation, "stretch_tensor": U,
        "stretch_lambda1": sv[..., 0], "stretch_lambda2": sv[..., 1],
        "isotropic_stretch": sv.mean(-1), "anisotropic_stretch": aniso,
        "bundle_q1": q1, "bundle_q2": q2, "bundle_q_abs": qabs,
        "bundle_q_angle": qangle, "bundle_det": det, "bundle_abs_det": np.abs(det),
        "parity_class": parity, "triangle_orientation_flip": flips,
        "fold_sign_change": sign_change, "fold_score": fold,
        "depth_gradient_u": gw_u, "depth_gradient_v": gw_v,
        "depth_gradient_abs": np.hypot(gw_u, gw_v),
        "depth_gradient_orientation": gw_angle,
        "depth_q_alignment": np.cos(2*(qangle-gw_angle)),
        "spin2_direction_real": zr, "spin2_direction_imag": zi,
        "spin2_direction_abs": np.hypot(zr, zi),
        "spin2_direction_angle": .5*np.arctan2(zi, zr),
        "ordinary_direction_real": orr, "ordinary_direction_imag": ori,
        "ordinary_direction_abs": np.hypot(orr, ori),
        "S": S, "D": D, "S_gradient_Q_alignment": asq,
        "D_gradient_Q_alignment": adq, "S_D_spin2_alignment": sd_align,
        "Q_S_weighted_q1": q1*S, "Q_S_weighted_q2": q2*S,
        "Q_D_weighted_q1": q1*D, "Q_D_weighted_q2": q2*D,
        "Q_SD_relational_q1": q1*sd_align, "Q_SD_relational_q2": q2*sd_align,
    }


def _cell_ids(u, v, resolution, bounds=(-8., 8.)):
    lo, hi = bounds
    ix = np.floor((u-lo)*resolution/(hi-lo)).astype(int)
    iy = np.floor((v-lo)*resolution/(hi-lo)).astype(int)
    ok = np.isfinite(u)&np.isfinite(v)&(ix>=0)&(iy>=0)&(ix<resolution)&(iy<resolution)
    return ix*resolution+iy, ok


def _component_statistics(cell, ok, shape, resolution):
    """Connected launch components per received cell using exact 4-adjacency."""
    n = cell.size; parent = np.arange(n, dtype=np.int64)
    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]; a = parent[a]
        return a
    def union(a, b):
        ra, rb = find(int(a)), find(int(b))
        if ra != rb: parent[rb] = ra
    ids = np.arange(n).reshape(shape); cg = cell.reshape(shape); vg = ok.reshape(shape)
    for a, b in ((ids[:, :-1], ids[:, 1:]), (ids[:-1, :], ids[1:, :])):
        eq = ok[a.ravel()] & ok[b.ravel()] & (cell[a.ravel()] == cell[b.ravel()])
        for aa, bb in zip(a.ravel()[eq], b.ravel()[eq]): union(aa, bb)
    roots = np.array([find(i) for i in range(n)])
    pairs = np.column_stack((cell[ok], roots[ok]))
    unique, counts = np.unique(pairs, axis=0, return_counts=True)
    comp = np.bincount(unique[:, 0], minlength=resolution**2)
    entropy = np.zeros(resolution**2)
    for c in np.unique(unique[:, 0]):
        w = counts[unique[:, 0] == c].astype(float); p = w/w.sum()
        entropy[c] = -np.sum(p*np.log(p+1e-30))/np.log(len(p)) if len(p)>1 else 0.
    return comp.reshape(resolution, resolution), entropy.reshape(resolution, resolution), roots


def deposit_bundle(bundle, resolution, bounds=(-8., 8.)):
    u=bundle["received_center_u"].ravel(); v=bundle["received_center_v"].ravel()
    cell, ok = _cell_ids(u, v, resolution, bounds)
    occ=np.bincount(cell[ok], minlength=resolution**2).reshape(resolution,resolution)
    def mean(name):
        a=np.asarray(bundle[name]).ravel(); good=ok&np.isfinite(a)
        s=np.bincount(cell[good],weights=a[good],minlength=resolution**2)
        c=np.bincount(cell[good],minlength=resolution**2)
        return np.divide(s,c,out=np.zeros_like(s),where=c>0).reshape(resolution,resolution)
    components, entropy, roots = _component_statistics(cell, ok, bundle["scale"].shape, resolution)
    parity_pos=mean("parity_class")
    # Mixed parity requires both signs, not a near-zero mean proxy.
    pos=np.bincount(cell[ok&(bundle["parity_class"].ravel()>0)],minlength=resolution**2)
    neg=np.bincount(cell[ok&(bundle["parity_class"].ravel()<0)],minlength=resolution**2)
    mixed=((pos>0)&(neg>0)).reshape(resolution,resolution)
    return {"occupancy":occ,"launch_component_count":components,
            "launch_origin_entropy":entropy,"bundle_q1":mean("bundle_q1"),
            "bundle_q2":mean("bundle_q2"),"bundle_q_abs":mean("bundle_q_abs"),
            "D":mean("D"),"S":mean("S"),"fold_density":mean("fold_score"),
            "parity_mean":parity_pos,"parity_mixed":mixed,
            "depth_gradient_abs":mean("depth_gradient_abs"),
            "ray_component_label":roots}


def received_overlap_diagnostic(bundle, deposited, resolution=64, bounds=(-8., 8.)):
    """Retain overlap and separation for disconnected launch components."""
    u=bundle["received_center_u"].ravel();v=bundle["received_center_v"].ravel()
    lu=bundle["launch_center_u"].ravel();lv=bundle["launch_center_v"].ravel()
    cell,ok=_cell_ids(u,v,resolution,bounds);roots=deposited["ray_component_label"];rows=[]
    for c in np.flatnonzero(deposited["launch_component_count"].ravel()>1):
        ids=np.flatnonzero(ok&(cell==c)); labels=np.unique(roots[ids]);cent=[]
        for label in labels:
            q=ids[roots[ids]==label];cent.append((lu[q].mean(),lv[q].mean(),u[q].mean(),v[q].mean()))
        for i in range(len(cent)):
            for j in range(i+1,len(cent)):
                a,b=cent[i],cent[j];rows.append((c,np.hypot(a[0]-b[0],a[1]-b[1]),1.,np.hypot(a[2]-b[2],a[3]-b[3])))
    arr=np.asarray(rows,float).reshape(-1,4)
    return {"received_cell_id":arr[:,0].astype(np.int32),"launch_separation":arr[:,1],
            "received_overlap_fraction":arr[:,2],"received_centroid_distance":arr[:,3]}


def downsample_mean(a, factor):
    a=np.asarray(a); n=a.shape[0]//factor
    return a[:n*factor,:n*factor].reshape(n,factor,n,factor).mean((1,3))


def pearson(a,b):
    a=np.asarray(a).ravel(); b=np.asarray(b).ravel(); m=np.isfinite(a)&np.isfinite(b)
    return float(np.corrcoef(a[m],b[m])[0,1]) if m.sum()>2 and np.std(a[m])*np.std(b[m])>0 else 0.


def resolution_metrics(maps):
    out={}
    for hi in (128,256):
        a,b=maps[64],maps[hi]; f=hi//64
        q1,q2=downsample_mean(b["bundle_q1"],f),downsample_mean(b["bundle_q2"],f)
        mag=np.hypot(q1,q2); base=np.hypot(a["bundle_q1"],a["bundle_q2"])
        orient=float(np.mean(np.cos(np.arctan2(q2,q1)-np.arctan2(a["bundle_q2"],a["bundle_q1"]))))
        out[f"64_vs_{hi}"]={"Q_map_Pearson":pearson(base,mag),"Q_orientation_agreement":orient,
          "fold_map_overlap":float(np.mean((a["fold_density"]>0)==(downsample_mean(b["fold_density"],f)>0))),
          "parity_boundary_overlap":float(np.mean(a["parity_mixed"]== (downsample_mean(b["parity_mixed"],f)>0))),
          "D_map_preservation":pearson(a["D"],downsample_mean(b["D"],f)),
          "launch_origin_entropy_preservation":pearson(a["launch_origin_entropy"],downsample_mean(b["launch_origin_entropy"],f))}
    m=out["64_vs_256"]
    if m["Q_map_Pearson"]>=.95 and m["Q_orientation_agreement"]>=.90: cls="OBSERVER_RESOLUTION_STABLE"
    elif m["Q_map_Pearson"]>=.80 and m["Q_orientation_agreement"]>=.70: cls="OBSERVER_RESOLUTION_MODERATE_SENSITIVITY"
    else: cls="OBSERVER_RESOLUTION_LIMITED"
    out["classification"]=cls
    return out


def multiscale_metrics(bundles):
    pairs={}; values=[]
    for a,b in zip(SCALES[:-1],SCALES[1:]):
        qa,qb=bundles[a],bundles[b]
        align=np.cos(2*(qa["bundle_q_angle"]-qb["bundle_q_angle"]))
        value=float(np.mean(align[qa["valid"]&qb["valid"]])); values.append(value)
        pairs[f"{a}_to_{b}"]={"orientation_persistence":value,
          "magnitude_pearson":pearson(qa["bundle_q_abs"],qb["bundle_q_abs"])}
    score=float(np.mean(values)); cls="SCALE_COHERENT" if score>=.8 else ("SCALE_TRANSITION" if score>=.3 else "SCALE_INCOHERENT")
    return {"neighboring_scales":pairs,"mean_orientation_persistence":score,"classification":cls}


def effective_rank(x):
    x=np.asarray(x,float).reshape(-1,np.asarray(x).shape[-1]); x=x[np.all(np.isfinite(x),1)]
    if len(x)<2:return 0.
    x=x-x.mean(0); s=np.linalg.svd(x,compute_uv=False); p=s*s/(np.sum(s*s)+1e-30)
    return float(np.exp(-np.sum(p[p>0]*np.log(p[p>0]))))


def json_dump(path, value):
    def safe(x):
        if isinstance(x,np.generic): return x.item()
        if isinstance(x,np.ndarray): return x.tolist()
        if isinstance(x,float) and not np.isfinite(x): return None
        raise TypeError(type(x).__name__)
    Path(path).write_text(json.dumps(value,indent=2,sort_keys=True,default=safe)+"\n")
