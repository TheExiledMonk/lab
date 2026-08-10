"""Target-blind construction and portable storage of a mixed 3-D observer volume."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import numpy as np

AXES = ("u", "v", "w")


def observer_normal(rays):
    n = np.cross(np.asarray(rays["e1"], float), np.asarray(rays["e2"], float))
    n /= np.linalg.norm(n)
    velocity = np.column_stack([rays["dx"], rays["dy"], rays["dz"]])
    # dx/dy/dz are already observer components: forward is positive dz.
    if np.nanmedian(velocity[:, 2]) < 0:
        n = -n
    return n


def ray_features(rays):
    """Return receipt coordinates and uncompressed primitive per-ray features."""
    u0, v0, uf, vf = (np.asarray(rays[k], np.float64) for k in ("u0", "v0", "uf", "vf"))
    e1, e2, n = np.asarray(rays["e1"], float), np.asarray(rays["e2"], float), observer_normal(rays)
    receipt = np.column_stack([rays["rx"], rays["ry"], rays["rz"]]).astype(np.float64)
    launch = np.column_stack([rays["launch_x"], rays["launch_y"], np.zeros_like(u0)])
    rw = receipt @ n
    lw = launch @ n
    delta = np.column_stack([uf-u0, vf-v0, rw-lw])
    direction = np.column_stack([rays["dx"], rays["dy"], rays["dz"]]).astype(np.float64)
    values = {}
    for i, a in enumerate(AXES):
        values[f"delta_{a}"] = delta[:, i]
        values[f"delta_{a}2"] = delta[:, i] ** 2
        values[f"dir_{a}"] = direction[:, i]
        values[f"dir_{a}2"] = direction[:, i] ** 2
    for i, a in enumerate(AXES):
        for j, b in enumerate(AXES[i+1:], i+1):
            values[f"delta_{a}{b}"] = delta[:, i] * delta[:, j]
            values[f"dir_{a}{b}"] = direction[:, i] * direction[:, j]
    for i, a in enumerate(AXES):
        for j, b in enumerate(AXES):
            values[f"delta_{a}_dir_{b}"] = delta[:, i] * direction[:, j]
    values["launch_receipt_separation"] = np.linalg.norm(receipt-launch, axis=1)
    values["received_depth"] = rw
    values["depth_change"] = rw-lw
    return np.column_stack([uf, vf, rw]), values, {"e1": e1, "e2": e2, "normal": n}


def _cic_indices(coords, bounds, shape):
    lo = np.array([bounds[a][0] for a in AXES]); hi = np.array([bounds[a][1] for a in AXES])
    scale = (np.array(shape)-1) / np.maximum(hi-lo, np.finfo(float).eps)
    x = (coords-lo)*scale; base = np.floor(x).astype(np.int64); frac = x-base
    for bits in range(8):
        off = np.array([(bits >> i) & 1 for i in range(3)])
        ijk = base + off
        weight = np.prod(np.where(off, frac, 1-frac), axis=1)
        valid = np.all((ijk >= 0) & (ijk < np.array(shape)), axis=1) & np.isfinite(weight)
        yield ijk[valid], weight[valid], valid


def construct_volume(rays, shape=(64,64,64), uv_bounds=(-8.,8.), method="trilinear_cic_3d"):
    if method not in ("trilinear_cic_3d", "nearest_voxel_3d"):
        raise ValueError(method)
    coords, features, basis = ray_features(rays)
    finite_w = coords[np.isfinite(coords[:,2]), 2]
    bounds = {"u": [float(uv_bounds[0]), float(uv_bounds[1])],
              "v": [float(uv_bounds[0]), float(uv_bounds[1])],
              "w": [float(finite_w.min()), float(finite_w.max())]}
    size = int(np.prod(shape)); occupancy = np.zeros(size); count = np.zeros(size, np.int64)
    sums = {k: np.zeros(size) for k in features}; sumsq = {k: np.zeros(size) for k in features}
    if method == "nearest_voxel_3d":
        lo=np.array([bounds[a][0] for a in AXES]); hi=np.array([bounds[a][1] for a in AXES])
        ijk=np.rint((coords-lo)*(np.array(shape)-1)/np.maximum(hi-lo,np.finfo(float).eps)).astype(int)
        valid=np.all((ijk>=0)&(ijk<np.array(shape)),axis=1)&np.all(np.isfinite(coords),axis=1)
        deposits=[(ijk[valid],np.ones(valid.sum()),valid)]
    else: deposits=_cic_indices(coords,bounds,shape)
    for ijk, weight, valid in deposits:
        flat=np.ravel_multi_index(ijk.T,shape); occupancy += np.bincount(flat,weight,minlength=size)
        for name, x in features.items():
            good=np.isfinite(x[valid]); f=flat[good]; wx=weight[good]*x[valid][good]
            sums[name] += np.bincount(f,wx,minlength=size); sumsq[name] += np.bincount(f,weight[good]*x[valid][good]**2,minlength=size)
    # Raw count is deliberately one ray/one voxel, independent of CIC support.
    lo=np.array([bounds[a][0] for a in AXES]); hi=np.array([bounds[a][1] for a in AXES])
    nearest=np.rint((coords-lo)*(np.array(shape)-1)/np.maximum(hi-lo,np.finfo(float).eps)).astype(int)
    valid=np.all((nearest>=0)&(nearest<np.array(shape)),axis=1)&np.all(np.isfinite(coords),axis=1)
    count += np.bincount(np.ravel_multi_index(nearest[valid].T,shape),minlength=size)
    arrays={"occupancy":occupancy.reshape(shape),"ray_count":count.reshape(shape)}
    for name in features:
        arrays[f"sum_{name}"]=sums[name].reshape(shape); arrays[f"sum_squared_{name}"]=sumsq[name].reshape(shape)
        arrays[f"mean_{name}"]=np.divide(sums[name],occupancy,out=np.full(size,np.nan),where=occupancy>0).reshape(shape)
    meta={"shape":list(shape),"bounds":bounds,"robust_w_percentiles":np.percentile(finite_w,[0,1,5,50,95,99,100]).tolist(),
          "deposition":method,"basis":{k:v.tolist() for k,v in basis.items()},"ray_count":int(len(coords)),
          "target_used_for_volume_construction":False,"coordinate_order":"u,v,w"}
    return arrays,meta


def manifest(arrays):
    out=[]
    for name,a in arrays.items():
        if name in ("occupancy","ray_count"): family="occupancy"
        elif "delta_" in name and "dir_" in name: family="displacement_direction_cross"
        elif "delta" in name: family="displacement"
        elif "dir" in name: family="direction"
        else: family="launch_receipt"
        norm="raw_count" if name=="ray_count" else ("weighted_sum" if name.startswith("sum_") else "occupancy_weighted_mean")
        sources=[]
        if "delta" in name: sources += ["u0", "v0", "uf", "vf", "received_position", "observer_normal"]
        if "dir" in name: sources += ["dx", "dy", "dz"]
        if "depth" in name or "separation" in name: sources += ["launch_position", "received_position"]
        if name in ("occupancy","ray_count"): sources=["received_observer_coordinates"]
        out.append({"name":name,"family":family,"units":"native observer coordinate","spin_class":"tensor_or_component","normalization":norm,"source_fields":sorted(set(sources)),"dtype":str(a.dtype)})
    return out


def save_volume(directory, arrays, metadata):
    directory=Path(directory); directory.mkdir(parents=True,exist_ok=True)
    path=directory/"observer_volume.npz"; np.savez_compressed(path,**arrays)
    (directory/"metadata.json").write_text(json.dumps(metadata,indent=2,sort_keys=True)+"\n")
    (directory/"channel_manifest.json").write_text(json.dumps(manifest(arrays),indent=2,sort_keys=True)+"\n")
    summary={k:{"finite":int(np.isfinite(v).sum()),"min":float(np.nanmin(v)) if np.isfinite(v).any() else None,"max":float(np.nanmax(v)) if np.isfinite(v).any() else None} for k,v in arrays.items()}
    (directory/"volume_summary.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    return hashlib.sha256(path.read_bytes()).hexdigest()
