"""M08/M09/M10 verified second-review implementation.

Scientific contract
-------------------
PS1-A : diagnostic source-local vector, v = P_i n.
PS1   : antisymmetrised difference, v = (P_i n - P_j n)/2.
PS1-B : backward-compatible alias of PS1. It is NOT an independent
        physics candidate until a distinct frozen equation exists.
PS2   : midpoint-symmetrised sum, v = (P_i n + P_j n)/2
        = ((P_i + P_j)/2) n.

PM1 normalises the selected vector before multiplying by A_ij.
PM2 preserves its raw magnitude.

M09 endpoint assembly assigns +R_ij at i and -R_ij at j, therefore
sum_i R_i = 0 while local energy may remain non-zero.

M10 interface rasterisation assigns +R_ij/2 to both adjacent voxels.
The rasteriser records boolean traversal masks so consumed-pair audits
measure the actual execution path and do not infer counts from non-zero
field values.
"""
from __future__ import annotations
import numpy as np

from .conventions import PS_LANES, PM_LANES

__all__ = [
    "build_pair_responses", "build_pair_responses_reference",
    "assemble_endpoint_field", "assemble_endpoint_field_reference",
    "rasterize_interface_field", "rasterize_interface_field_reference",
    "expected_interface_pair_count", "consumed_interface_pair_count",
    "interface_pair_count_audit", "PS_LANES", "PM_LANES",
    "PM1", "PM2", "PS1_A", "PS1", "PS1_B", "PS2",
    "PS_EQUIVALENCE_CLASS", "PS_PHYSICS_CANDIDATE",
    "PairTransferError",
]

class PairTransferError(ValueError):
    pass

PM1 = "PM1"
PM2 = "PM2"
PS1_A = "PS1-A"
PS1 = "PS1"
PS1_B = "PS1-B"
PS2 = "PS2"

PS_EQUIVALENCE_CLASS = {
    PS1_A: "PS1-A_DIAGNOSTIC",
    PS1: "PS1_EQ_PS1-B",
    PS1_B: "PS1_EQ_PS1-B",
    PS2: "PS2",
}
PS_PHYSICS_CANDIDATE = {
    PS1_A: False,
    PS1: True,
    PS1_B: False,
    PS2: True,
}

_REQUIRED_RESPONSE_KEYS = (
    "R_ij_xp", "R_ij_y_xp", "R_ij_z_xp",
    "R_ij_yp", "R_ij_y_yp", "R_ij_z_yp",
    "R_ij_zp", "R_ij_y_zp", "R_ij_z_zp",
)


def _validate_projector(projector_field):
    if len(projector_field) != 6:
        raise PairTransferError("projector_field must contain 6 symmetric-tensor components")
    arrays = tuple(np.asarray(a, dtype=np.float64) for a in projector_field)
    shape = arrays[0].shape
    if any(a.shape != shape for a in arrays):
        raise PairTransferError("all projector components must share shape")
    return arrays, shape


def _validate_amplitudes(pair_amplitudes, shape):
    out = {}
    for key in ("A_xp", "A_yp", "A_zp"):
        if key not in pair_amplitudes:
            raise PairTransferError(f"missing pair amplitude {key}")
        a = np.asarray(pair_amplitudes[key], dtype=np.float64)
        if a.shape != shape:
            raise PairTransferError(f"{key} shape {a.shape} != projector shape {shape}")
        out[key] = a
    return out


def _partner(arr, array_axis):
    """Return arr(j=i+1) stored at source i; invalid boundary is zero."""
    out = np.zeros_like(arr)
    src = [slice(None)] * 3
    dst = [slice(None)] * 3
    src[array_axis] = slice(0, arr.shape[array_axis] - 1)
    dst[array_axis] = slice(1, arr.shape[array_axis])
    out[tuple(src)] = arr[tuple(dst)]
    return out


def _project_n(P, axis):
    Pxx, Pxy, Pxz, Pyy, Pyz, Pzz = P
    if axis == "xp":
        return Pxx, Pxy, Pxz
    if axis == "yp":
        return Pxy, Pyy, Pyz
    if axis == "zp":
        return Pxz, Pyz, Pzz
    raise PairTransferError(f"unknown pair axis {axis!r}")


def _normalise(vx, vy, vz):
    """Normalise only non-zero vectors; exact zero remains exact zero."""
    mag = np.sqrt(vx * vx + vy * vy + vz * vz)
    nz = mag > 0.0
    ox = np.zeros_like(vx); oy = np.zeros_like(vy); oz = np.zeros_like(vz)
    np.divide(vx, mag, out=ox, where=nz)
    np.divide(vy, mag, out=oy, where=nz)
    np.divide(vz, mag, out=oz, where=nz)
    return ox, oy, oz


def _pair_vector(P, axis, lane):
    array_axis = {"xp": 2, "yp": 1, "zp": 0}[axis]
    vi = _project_n(P, axis)
    Pj = tuple(_partner(c, array_axis) for c in P)
    vj = _project_n(Pj, axis)

    if lane == PS1_A:
        return vi
    if lane in (PS1, PS1_B):
        return tuple(0.5 * (a - b) for a, b in zip(vi, vj))
    if lane == PS2:
        return tuple(0.5 * (a + b) for a, b in zip(vi, vj))
    raise PairTransferError(f"unknown pair symmetrisation {lane!r}")


def _valid_source_mask(shape, axis):
    mask = np.zeros(shape, dtype=bool)
    if axis == "xp" and shape[2] >= 2:
        mask[:, :, :-1] = True
    elif axis == "yp" and shape[1] >= 2:
        mask[:, :-1, :] = True
    elif axis == "zp" and shape[0] >= 2:
        mask[:-1, :, :] = True
    return mask


def build_pair_responses(pair_registry, pair_amplitudes, projector_field,
                         magnitude_formulation=PM1, pair_symmetrization=PS2):
    """Construct the positive-N6 pair response arrays.

    ``pair_registry`` is accepted for API compatibility and checked
    against the geometric pair count when supplied. The array layout is
    (z,y,x); vector components are (x,y,z).
    """
    P, shape = _validate_projector(projector_field)
    A = _validate_amplitudes(pair_amplitudes, shape)
    if magnitude_formulation not in (PM1, PM2):
        raise PairTransferError(f"unknown magnitude formulation {magnitude_formulation!r}")
    if pair_symmetrization not in (PS1_A, PS1, PS1_B, PS2):
        raise PairTransferError(f"unknown pair symmetrisation {pair_symmetrization!r}")
    if pair_registry is not None and len(pair_registry) != expected_interface_pair_count(shape):
        raise PairTransferError(
            f"pair_registry has {len(pair_registry)} pairs; expected {expected_interface_pair_count(shape)}"
        )

    result = {}
    stats = {}
    for axis, akey, outkeys in (
        ("xp", "A_xp", ("R_ij_xp", "R_ij_y_xp", "R_ij_z_xp")),
        ("yp", "A_yp", ("R_ij_yp", "R_ij_y_yp", "R_ij_z_yp")),
        ("zp", "A_zp", ("R_ij_zp", "R_ij_y_zp", "R_ij_z_zp")),
    ):
        vx, vy, vz = _pair_vector(P, axis, pair_symmetrization)
        if magnitude_formulation == PM1:
            vx, vy, vz = _normalise(vx, vy, vz)
        valid = _valid_source_mask(shape, axis)
        amp = A[akey]
        rx = np.where(valid, amp * vx, 0.0)
        ry = np.where(valid, amp * vy, 0.0)
        rz = np.where(valid, amp * vz, 0.0)
        for key, arr in zip(outkeys, (rx, ry, rz)):
            result[key] = arr
        mag = np.sqrt(rx * rx + ry * ry + rz * rz)
        stats[axis] = {
            "R_min": float(mag.min()),
            "R_max": float(mag.max()),
            "R_rms": float(np.sqrt(np.mean(mag * mag))),
            "R_abs_sum": float(np.sum(mag)),
            "n_nonzero": int(np.count_nonzero(mag)),
            "valid_pair_slots": int(np.count_nonzero(valid)),
        }
    result["statistics"] = stats
    result["pair_symmetrization"] = pair_symmetrization
    result["candidate_equivalence_class"] = PS_EQUIVALENCE_CLASS[pair_symmetrization]
    result["physics_candidate"] = PS_PHYSICS_CANDIDATE[pair_symmetrization]
    return result


def build_pair_responses_reference(pair_registry, pair_amplitudes, projector_field,
                                   magnitude_formulation=PM1, pair_symmetrization=PS2):
    """Independent explicit pair-by-pair reference implementation."""
    P, shape = _validate_projector(projector_field)
    A = _validate_amplitudes(pair_amplitudes, shape)
    if pair_registry is None:
        raise PairTransferError("reference implementation requires pair_registry")
    arrays = {k: np.zeros(shape, dtype=np.float64) for k in _REQUIRED_RESPONSE_KEYS}
    Pxx, Pxy, Pxz, Pyy, Pyz, Pzz = P

    def mat_at(idx):
        z,y,x = idx
        return np.array([[Pxx[z,y,x],Pxy[z,y,x],Pxz[z,y,x]],
                         [Pxy[z,y,x],Pyy[z,y,x],Pyz[z,y,x]],
                         [Pxz[z,y,x],Pyz[z,y,x],Pzz[z,y,x]]], dtype=np.float64)

    amap = {"xp":"A_xp", "yp":"A_yp", "zp":"A_zp"}
    outmap = {
        "xp": ("R_ij_xp","R_ij_y_xp","R_ij_z_xp"),
        "yp": ("R_ij_yp","R_ij_y_yp","R_ij_z_yp"),
        "zp": ("R_ij_zp","R_ij_y_zp","R_ij_z_zp"),
    }
    for pair in pair_registry:
        i = pair.i_index; j = pair.j_index
        n = np.asarray(pair.direction_xyz, dtype=np.float64)
        vi = mat_at(i) @ n
        vj = mat_at(j) @ n
        if pair_symmetrization == PS1_A:
            v = vi
        elif pair_symmetrization in (PS1, PS1_B):
            v = 0.5 * (vi - vj)
        elif pair_symmetrization == PS2:
            v = 0.5 * (vi + vj)
        else:
            raise PairTransferError(pair_symmetrization)
        if magnitude_formulation == PM1:
            m = float(np.linalg.norm(v))
            v = v / m if m > 0.0 else np.zeros(3)
        elif magnitude_formulation != PM2:
            raise PairTransferError(magnitude_formulation)
        amp = float(A[amap[pair.axis]][i])
        r = amp * v
        for key, value in zip(outmap[pair.axis], r):
            arrays[key][i] = value
    arrays["statistics"] = {}
    arrays["pair_symmetrization"] = pair_symmetrization
    arrays["candidate_equivalence_class"] = PS_EQUIVALENCE_CLASS[pair_symmetrization]
    arrays["physics_candidate"] = PS_PHYSICS_CANDIDATE[pair_symmetrization]
    return arrays


def _validate_pair_responses(pair_responses, shape):
    out = {}
    for key in _REQUIRED_RESPONSE_KEYS:
        if key not in pair_responses:
            raise PairTransferError(f"missing response array {key}")
        a = np.asarray(pair_responses[key], dtype=np.float64)
        if a.shape != tuple(shape):
            raise PairTransferError(f"{key} shape {a.shape} != {tuple(shape)}")
        out[key] = a
    return out


def assemble_endpoint_field(pair_responses, shape):
    """Assign +R_ij at source i and -R_ij at partner j."""
    r = _validate_pair_responses(pair_responses, shape)
    Rx=np.zeros(shape); Ry=np.zeros(shape); Rz=np.zeros(shape)
    for keys, sax in (
        (("R_ij_xp","R_ij_y_xp","R_ij_z_xp"),2),
        (("R_ij_yp","R_ij_y_yp","R_ij_z_yp"),1),
        (("R_ij_zp","R_ij_y_zp","R_ij_z_zp"),0),
    ):
        src=[slice(None)]*3; dst=[slice(None)]*3
        src[sax]=slice(0,shape[sax]-1); dst[sax]=slice(1,shape[sax])
        src=tuple(src); dst=tuple(dst)
        for out,key in ((Rx,keys[0]),(Ry,keys[1]),(Rz,keys[2])):
            out[src] += r[key][src]
            out[dst] -= r[key][src]
    mag=np.sqrt(Rx*Rx+Ry*Ry+Rz*Rz)
    sum_vec=(float(Rx.sum()),float(Ry.sum()),float(Rz.sum()))
    return {"Rx_3d":Rx,"Ry_3d":Ry,"Rz_3d":Rz,"statistics":{
        "Rx_rms":float(np.sqrt(np.mean(Rx*Rx))),
        "Ry_rms":float(np.sqrt(np.mean(Ry*Ry))),
        "Rz_rms":float(np.sqrt(np.mean(Rz*Rz))),
        "total_rms":float(np.sqrt(np.mean(mag*mag))),
        "max_vector_norm":float(mag.max()),"n_nonzero":int(np.count_nonzero(mag)),
        "global_vector_sum":sum_vec,
        "global_vector_sum_norm":float(np.linalg.norm(sum_vec)),
        "endpoint_energy":float(np.sum(mag*mag)),
    }}


def assemble_endpoint_field_reference(pair_responses, shape):
    """Explicit source/partner loop reference for M09."""
    r=_validate_pair_responses(pair_responses,shape)
    Rx=np.zeros(shape); Ry=np.zeros(shape); Rz=np.zeros(shape)
    nz,ny,nx=shape
    for z in range(nz):
        for y in range(ny):
            for x in range(nx):
                for ok,j,keys in (
                    (x<nx-1,(z,y,x+1),("R_ij_xp","R_ij_y_xp","R_ij_z_xp")),
                    (y<ny-1,(z,y+1,x),("R_ij_yp","R_ij_y_yp","R_ij_z_yp")),
                    (z<nz-1,(z+1,y,x),("R_ij_zp","R_ij_y_zp","R_ij_z_zp")),
                ):
                    if not ok: continue
                    i=(z,y,x)
                    vals=(r[keys[0]][i],r[keys[1]][i],r[keys[2]][i])
                    for out,val in zip((Rx,Ry,Rz),vals):
                        out[i]+=val; out[j]-=val
    mag=np.sqrt(Rx*Rx+Ry*Ry+Rz*Rz)
    sv=(float(Rx.sum()),float(Ry.sum()),float(Rz.sum()))
    return {"Rx_3d":Rx,"Ry_3d":Ry,"Rz_3d":Rz,"statistics":{
        "global_vector_sum":sv,"global_vector_sum_norm":float(np.linalg.norm(sv)),
        "endpoint_energy":float(np.sum(mag*mag))}}


def rasterize_interface_field(pair_responses, shape):
    """Assign +R_ij/2 to both adjacent voxels and expose traversal masks."""
    r=_validate_pair_responses(pair_responses,shape)
    Rx=np.zeros(shape); Ry=np.zeros(shape); Rz=np.zeros(shape)
    masks={a:np.zeros(shape,dtype=bool) for a in ("xp","yp","zp")}
    for axis,keys,sax in (
        ("xp",("R_ij_xp","R_ij_y_xp","R_ij_z_xp"),2),
        ("yp",("R_ij_yp","R_ij_y_yp","R_ij_z_yp"),1),
        ("zp",("R_ij_zp","R_ij_y_zp","R_ij_z_zp"),0),
    ):
        src=[slice(None)]*3; dst=[slice(None)]*3
        src[sax]=slice(0,shape[sax]-1); dst[sax]=slice(1,shape[sax])
        src=tuple(src); dst=tuple(dst); masks[axis][src]=True
        for out,key in ((Rx,keys[0]),(Ry,keys[1]),(Rz,keys[2])):
            out[src]+=0.5*r[key][src]; out[dst]+=0.5*r[key][src]
    mag=np.sqrt(Rx*Rx+Ry*Ry+Rz*Rz)
    sv=(float(Rx.sum()),float(Ry.sum()),float(Rz.sum()))
    counts={a:int(np.count_nonzero(m)) for a,m in masks.items()}
    return {"Rx_3d_interface":Rx,"Ry_3d_interface":Ry,"Rz_3d_interface":Rz,
            "consumed_pair_masks":masks,"statistics":{
                "Rx_rms":float(np.sqrt(np.mean(Rx*Rx))),
                "Ry_rms":float(np.sqrt(np.mean(Ry*Ry))),
                "Rz_rms":float(np.sqrt(np.mean(Rz*Rz))),
                "total_rms":float(np.sqrt(np.mean(mag*mag))),
                "max_vector_norm":float(mag.max()),"n_nonzero":int(np.count_nonzero(mag)),
                "global_vector_sum":sv,"interface_energy":float(np.sum(mag*mag)),
                "consumed_pair_count_xp":counts["xp"],
                "consumed_pair_count_yp":counts["yp"],
                "consumed_pair_count_zp":counts["zp"],
                "consumed_pair_count_total":sum(counts.values()),
            }}


def rasterize_interface_field_reference(pair_responses, shape):
    """Explicit pair-slot loop reference for M10."""
    r=_validate_pair_responses(pair_responses,shape)
    Rx=np.zeros(shape); Ry=np.zeros(shape); Rz=np.zeros(shape)
    nz,ny,nx=shape
    for z in range(nz):
        for y in range(ny):
            for x in range(nx):
                for ok,j,keys in (
                    (x<nx-1,(z,y,x+1),("R_ij_xp","R_ij_y_xp","R_ij_z_xp")),
                    (y<ny-1,(z,y+1,x),("R_ij_yp","R_ij_y_yp","R_ij_z_yp")),
                    (z<nz-1,(z+1,y,x),("R_ij_zp","R_ij_y_zp","R_ij_z_zp")),
                ):
                    if not ok: continue
                    i=(z,y,x); vals=(r[keys[0]][i],r[keys[1]][i],r[keys[2]][i])
                    for out,val in zip((Rx,Ry,Rz),vals):
                        out[i]+=0.5*val; out[j]+=0.5*val
    mag=np.sqrt(Rx*Rx+Ry*Ry+Rz*Rz); sv=(float(Rx.sum()),float(Ry.sum()),float(Rz.sum()))
    return {"Rx_3d_interface":Rx,"Ry_3d_interface":Ry,"Rz_3d_interface":Rz,
            "statistics":{"global_vector_sum":sv,"interface_energy":float(np.sum(mag*mag))}}


def expected_interface_pair_count(shape):
    nz,ny,nx=shape
    return int(nz*ny*max(nx-1,0)+nz*max(ny-1,0)*nx+max(nz-1,0)*ny*nx)


def consumed_interface_pair_count(pair_responses, shape):
    iface=rasterize_interface_field(pair_responses,shape)
    return int(sum(np.count_nonzero(m) for m in iface["consumed_pair_masks"].values()))


def interface_pair_count_audit(pair_responses, shape):
    nz,ny,nx=shape
    iface=rasterize_interface_field(pair_responses,shape)
    expected={"xp":nz*ny*max(nx-1,0),"yp":nz*max(ny-1,0)*nx,"zp":max(nz-1,0)*ny*nx}
    rows=[]
    for axis in ("xp","yp","zp"):
        got=int(np.count_nonzero(iface["consumed_pair_masks"][axis])); exp=int(expected[axis])
        rows.append({"axis":axis,"expected_pair_count":exp,"consumed_pair_count":got,
                     "omitted_pair_count":max(exp-got,0),"duplicated_pair_count":max(got-exp,0),
                     "passes":got==exp})
    exp=sum(expected.values()); got=sum(r["consumed_pair_count"] for r in rows)
    rows.append({"axis":"TOTAL","expected_pair_count":int(exp),"consumed_pair_count":int(got),
                 "omitted_pair_count":max(int(exp-got),0),"duplicated_pair_count":max(int(got-exp),0),
                 "passes":got==exp and all(r["passes"] for r in rows)})
    return rows
