"""M14 — LOS Projection.

Correction pass FOUNDATION-001-CORRECTION-001
---------------------------------------------
* The previous docstring said "the depth component is discarded" but
  the function returned all three summed components.  This
  contradiction is resolved by splitting the API into two functions:
    - ``project_vector_los_full`` returns the FULL globally-ordered
      components (Rx_sum, Ry_sum, Rz_sum).
    - ``project_vector_to_image_plane`` returns ONLY the two
      transverse components in a documented order.
* Metadata (``los_axis``, ``depth_array_axis``, ``image_component_1``,
  ``image_component_2``, ``output_plane_axis_order``) is attached to
  every projection.
"""
from __future__ import annotations
import numpy as np

__all__ = [
    "project_vector_los_full",
    "project_vector_los_full_reference",
    "project_vector_to_image_plane",
    "project_vector_to_image_plane_reference",
    "LOSProjectionError",
]


class LOSProjectionError(ValueError):
    pass


# Mapping from los_axis to:
#   - the array axis along which we sum (depth_array_axis)
#   - the two kept components and their order
#   - the output plane axis order
_LOS_CONFIG = {
    "z": {
        "depth_array_axis": 0,
        "image_component_1": ("x", 0),  # (label, stacked-component-index)
        "image_component_2": ("y", 1),
        "output_plane_axis_order": ("y", "x"),
    },
    "y": {
        "depth_array_axis": 1,
        "image_component_1": ("x", 0),
        "image_component_2": ("z", 2),
        "output_plane_axis_order": ("z", "x"),
    },
    "x": {
        "depth_array_axis": 2,
        "image_component_1": ("y", 1),
        "image_component_2": ("z", 2),
        "output_plane_axis_order": ("z", "y"),
    },
}


def project_vector_los_full(Rx, Ry, Rz, los_axis="z"):
    """Full LOS diagnostic projection.

    Sums all three globally-ordered components over the depth axis.
    Returns ``(Rx_sum, Ry_sum, Rz_sum)`` together with metadata
    describing the projection.
    """
    if los_axis not in _LOS_CONFIG:
        raise LOSProjectionError(f"unsupported los_axis: {los_axis!r}")
    Rx = np.asarray(Rx, dtype=np.float64)
    Ry = np.asarray(Ry, dtype=np.float64)
    Rz = np.asarray(Rz, dtype=np.float64)
    if Rx.shape != Ry.shape or Rx.shape != Rz.shape:
        raise LOSProjectionError("Rx, Ry, Rz must share the same shape")
    depth = _LOS_CONFIG[los_axis]["depth_array_axis"]
    Rx_sum = np.sum(Rx, axis=depth)
    Ry_sum = np.sum(Ry, axis=depth)
    Rz_sum = np.sum(Rz, axis=depth)
    return {
        "Rx_sum": Rx_sum, "Ry_sum": Ry_sum, "Rz_sum": Rz_sum,
        "los_axis": los_axis,
        "depth_array_axis": depth,
        "image_component_1": _LOS_CONFIG[los_axis]["image_component_1"][0],
        "image_component_2": _LOS_CONFIG[los_axis]["image_component_2"][0],
        "output_plane_axis_order": _LOS_CONFIG[los_axis]["output_plane_axis_order"],
    }


def project_vector_los_full_reference(Rx, Ry, Rz, los_axis="z"):
    """Reference full LOS projection using an explicit depth-loop."""
    if los_axis not in _LOS_CONFIG:
        raise LOSProjectionError(f"unsupported los_axis: {los_axis!r}")
    Rx = np.asarray(Rx, dtype=np.float64)
    Ry = np.asarray(Ry, dtype=np.float64)
    Rz = np.asarray(Rz, dtype=np.float64)
    if Rx.shape != Ry.shape or Rx.shape != Rz.shape:
        raise LOSProjectionError("Rx, Ry, Rz must share the same shape")
    nz, ny, nx = Rx.shape
    if los_axis == "z":
        Rx_sum = np.zeros((ny, nx))
        Ry_sum = np.zeros((ny, nx))
        Rz_sum = np.zeros((ny, nx))
        for iz in range(nz):
            Rx_sum += Rx[iz]; Ry_sum += Ry[iz]; Rz_sum += Rz[iz]
    elif los_axis == "y":
        Rx_sum = np.zeros((nz, nx))
        Ry_sum = np.zeros((nz, nx))
        Rz_sum = np.zeros((nz, nx))
        for iy in range(ny):
            Rx_sum += Rx[:, iy]; Ry_sum += Ry[:, iy]; Rz_sum += Rz[:, iy]
    else:
        Rx_sum = np.zeros((nz, ny))
        Ry_sum = np.zeros((nz, ny))
        Rz_sum = np.zeros((nz, ny))
        for ix in range(nx):
            Rx_sum += Rx[:, :, ix]; Ry_sum += Ry[:, :, ix]; Rz_sum += Rz[:, :, ix]
    return {
        "Rx_sum": Rx_sum, "Ry_sum": Ry_sum, "Rz_sum": Rz_sum,
        "los_axis": los_axis,
        "depth_array_axis": _LOS_CONFIG[los_axis]["depth_array_axis"],
        "image_component_1": _LOS_CONFIG[los_axis]["image_component_1"][0],
        "image_component_2": _LOS_CONFIG[los_axis]["image_component_2"][0],
        "output_plane_axis_order": _LOS_CONFIG[los_axis]["output_plane_axis_order"],
    }


def project_vector_to_image_plane(Rx, Ry, Rz, los_axis="z"):
    """Project a 3D vector field to its image plane.

    Returns only the two transverse components in a documented order:

      LOS axis "z" → (Rx_sum, Ry_sum)    # discard Rz (depth)
      LOS axis "y" → (Rx_sum, Rz_sum)    # discard Ry
      LOS axis "x" → (Ry_sum, Rz_sum)    # discard Rx
    """
    full = project_vector_los_full(Rx, Ry, Rz, los_axis)
    c1 = full["image_component_1"]
    c2 = full["image_component_2"]
    comp_map = {"x": full["Rx_sum"], "y": full["Ry_sum"], "z": full["Rz_sum"]}
    return {
        "comp_1": comp_map[c1],  # first kept component
        "comp_1_label": c1,
        "comp_2": comp_map[c2],  # second kept component
        "comp_2_label": c2,
        "los_axis": los_axis,
        "depth_array_axis": full["depth_array_axis"],
        "output_plane_axis_order": full["output_plane_axis_order"],
    }


def project_vector_to_image_plane_reference(Rx, Ry, Rz, los_axis="z"):
    """Reference image-plane projection using an explicit depth-loop
    (only the two kept components are accumulated)."""
    if los_axis not in _LOS_CONFIG:
        raise LOSProjectionError(f"unsupported los_axis: {los_axis!r}")
    Rx = np.asarray(Rx, dtype=np.float64)
    Ry = np.asarray(Ry, dtype=np.float64)
    Rz = np.asarray(Rz, dtype=np.float64)
    nz, ny, nx = Rx.shape
    cfg = _LOS_CONFIG[los_axis]
    c1_lbl, c1_idx = cfg["image_component_1"]
    c2_lbl, c2_idx = cfg["image_component_2"]
    comps = (Rx, Ry, Rz)
    if los_axis == "z":
        out_shape = (ny, nx)
    elif los_axis == "y":
        out_shape = (nz, nx)
    else:
        out_shape = (nz, ny)
    c1_sum = np.zeros(out_shape)
    c2_sum = np.zeros(out_shape)
    if los_axis == "z":
        for iz in range(nz):
            c1_sum += comps[c1_idx][iz]
            c2_sum += comps[c2_idx][iz]
    elif los_axis == "y":
        for iy in range(ny):
            c1_sum += comps[c1_idx][:, iy]
            c2_sum += comps[c2_idx][:, iy]
    else:
        for ix in range(nx):
            c1_sum += comps[c1_idx][:, :, ix]
            c2_sum += comps[c2_idx][:, :, ix]
    return {
        "comp_1": c1_sum, "comp_1_label": c1_lbl,
        "comp_2": c2_sum, "comp_2_label": c2_lbl,
        "los_axis": los_axis,
        "depth_array_axis": cfg["depth_array_axis"],
        "output_plane_axis_order": cfg["output_plane_axis_order"],
    }


# ----------------------------------------------------------------------
# Self-check
# ----------------------------------------------------------------------
def _constant_field_full_test():
    nz, ny, nx = 4, 5, 6
    c = 0.3
    Rx = c * np.ones((nz, ny, nx))
    Ry = c * np.ones((nz, ny, nx))
    Rz = c * np.ones((nz, ny, nx))
    out = project_vector_los_full(Rx, Ry, Rz, "z")
    expected = c * nz
    err = max(np.max(np.abs(out["Rx_sum"] - expected)),
              np.max(np.abs(out["Ry_sum"] - expected)),
              np.max(np.abs(out["Rz_sum"] - expected)))
    return {"los_max_err": err, "passes": err == 0.0}


def _image_plane_constant_test():
    nz, ny, nx = 4, 5, 6
    c = 0.3
    Rx = c * np.ones((nz, ny, nx))
    Ry = c * np.ones((nz, ny, nx))
    Rz = c * np.ones((nz, ny, nx))
    out = project_vector_to_image_plane(Rx, Ry, Rz, "z")
    expected = c * nz
    err = max(np.max(np.abs(out["comp_1"] - expected)),
              np.max(np.abs(out["comp_2"] - expected)))
    return {"image_plane_err": err, "passes": err == 0.0}


def _image_plane_components_test():
    """For LOS axis 'z', the kept components are Rx and Ry.
    For LOS axis 'y', the kept components are Rx and Rz.
    For LOS axis 'x', the kept components are Ry and Rz."""
    rows = []
    nz, ny, nx = 3, 4, 5
    rng = np.random.RandomState(0)
    Rx = rng.randn(nz, ny, nx)
    Ry = rng.randn(nz, ny, nx)
    Rz = rng.randn(nz, ny, nx)
    full_z = project_vector_los_full(Rx, Ry, Rz, "z")
    full_y = project_vector_los_full(Rx, Ry, Rz, "y")
    full_x = project_vector_los_full(Rx, Ry, Rz, "x")
    ip_z = project_vector_to_image_plane(Rx, Ry, Rz, "z")
    ip_y = project_vector_to_image_plane(Rx, Ry, Rz, "y")
    ip_x = project_vector_to_image_plane(Rx, Ry, Rz, "x")
    rows.append({
        "los_axis": "z",
        "comp_1_label": ip_z["comp_1_label"],
        "comp_2_label": ip_z["comp_2_label"],
        "comp_1_matches_full": np.allclose(ip_z["comp_1"], full_z["Rx_sum"]),
        "comp_2_matches_full": np.allclose(ip_z["comp_2"], full_z["Ry_sum"]),
        "passes": (ip_z["comp_1_label"] == "x"
                    and ip_z["comp_2_label"] == "y"
                    and np.allclose(ip_z["comp_1"], full_z["Rx_sum"])
                    and np.allclose(ip_z["comp_2"], full_z["Ry_sum"])),
    })
    rows.append({
        "los_axis": "y",
        "comp_1_label": ip_y["comp_1_label"],
        "comp_2_label": ip_y["comp_2_label"],
        "comp_1_matches_full": np.allclose(ip_y["comp_1"], full_y["Rx_sum"]),
        "comp_2_matches_full": np.allclose(ip_y["comp_2"], full_y["Rz_sum"]),
        "passes": (ip_y["comp_1_label"] == "x"
                    and ip_y["comp_2_label"] == "z"
                    and np.allclose(ip_y["comp_1"], full_y["Rx_sum"])
                    and np.allclose(ip_y["comp_2"], full_y["Rz_sum"])),
    })
    rows.append({
        "los_axis": "x",
        "comp_1_label": ip_x["comp_1_label"],
        "comp_2_label": ip_x["comp_2_label"],
        "comp_1_matches_full": np.allclose(ip_x["comp_1"], full_x["Ry_sum"]),
        "comp_2_matches_full": np.allclose(ip_x["comp_2"], full_x["Rz_sum"]),
        "passes": (ip_x["comp_1_label"] == "y"
                    and ip_x["comp_2_label"] == "z"
                    and np.allclose(ip_x["comp_1"], full_x["Ry_sum"])
                    and np.allclose(ip_x["comp_2"], full_x["Rz_sum"])),
    })
    return {"results": rows,
            "passes": all(r["passes"] for r in rows)}


def _antisymmetric_depth_test():
    nz, ny, nx = 5, 4, 6
    rng = np.random.RandomState(0)
    Rx = rng.randn(nz, ny, nx)
    Rx = Rx - np.flip(Rx, axis=0)
    Ry = rng.randn(nz, ny, nx)
    Ry = Ry - np.flip(Ry, axis=0)
    Rz = rng.randn(nz, ny, nx)
    Rz = Rz - np.flip(Rz, axis=0)
    full = project_vector_los_full(Rx, Ry, Rz, "z")
    return {"max_los": float(max(np.max(np.abs(full["Rx_sum"])),
                                    np.max(np.abs(full["Ry_sum"])),
                                    np.max(np.abs(full["Rz_sum"])))),
            "passes": float(max(np.max(np.abs(full["Rx_sum"])),
                                  np.max(np.abs(full["Ry_sum"])),
                                  np.max(np.abs(full["Rz_sum"])))) < 1e-12}


def _single_slice_test():
    nz, ny, nx = 5, 4, 6
    rng = np.random.RandomState(0)
    Rx = rng.randn(nz, ny, nx); Rx[1:] = 0
    Ry = rng.randn(nz, ny, nx); Ry[1:] = 0
    Rz = rng.randn(nz, ny, nx); Rz[1:] = 0
    full = project_vector_los_full(Rx, Ry, Rz, "z")
    err = max(np.max(np.abs(full["Rx_sum"] - Rx[0])),
              np.max(np.abs(full["Ry_sum"] - Ry[0])),
              np.max(np.abs(full["Rz_sum"] - Rz[0])))
    return {"err": float(err), "passes": float(err) == 0.0}


def _zero_field_test():
    Rx = np.zeros((4, 5, 6)); Ry = np.zeros((4, 5, 6)); Rz = np.zeros((4, 5, 6))
    full = project_vector_los_full(Rx, Ry, Rz, "z")
    return {"passes": (full["Rx_sum"].sum() == 0.0
                        and full["Ry_sum"].sum() == 0.0
                        and full["Rz_sum"].sum() == 0.0)}


def _production_vs_reference_test():
    nz, ny, nx = 5, 6, 7
    rng = np.random.RandomState(0)
    Rx = rng.randn(nz, ny, nx); Ry = rng.randn(nz, ny, nx); Rz = rng.randn(nz, ny, nx)
    errs = []
    for axis in ("z", "y", "x"):
        p = project_vector_los_full(Rx, Ry, Rz, axis)
        r = project_vector_los_full_reference(Rx, Ry, Rz, axis)
        errs.append(max(float(np.max(np.abs(p[k] - r[k])))
                         for k in ("Rx_sum", "Ry_sum", "Rz_sum")))
        # image-plane agreement
        p_ip = project_vector_to_image_plane(Rx, Ry, Rz, axis)
        r_ip = project_vector_to_image_plane_reference(Rx, Ry, Rz, axis)
        errs.append(max(float(np.max(np.abs(p_ip["comp_1"] - r_ip["comp_1"]))),
                         float(np.max(np.abs(p_ip["comp_2"] - r_ip["comp_2"])))))
    err = max(errs)
    return {"max_diff": err, "passes": err < 1e-14}


def _known_cancellation_test():
    nz, ny, nx = 5, 4, 6
    rng = np.random.RandomState(0)
    R0 = rng.randn(ny, nx)
    Rx = np.zeros((nz, ny, nx))
    Ry = np.zeros((nz, ny, nx))
    Rz = np.zeros((nz, ny, nx))
    Rx[2] = R0; Rx[3] = -R0
    full = project_vector_los_full(Rx, Ry, Rz, "z")
    return {"max_los": float(max(np.max(np.abs(full["Rx_sum"])),
                                    np.max(np.abs(full["Ry_sum"])),
                                    np.max(np.abs(full["Rz_sum"])))),
            "passes": float(max(np.max(np.abs(full["Rx_sum"])),
                                  np.max(np.abs(full["Ry_sum"])),
                                  np.max(np.abs(full["Rz_sum"])))) < 1e-12}


if __name__ == "__main__":
    r = _constant_field_full_test(); assert r["passes"], r
    print("M14 full constant: pass")
    r = _image_plane_constant_test(); assert r["passes"], r
    print("M14 image-plane constant: pass")
    r = _image_plane_components_test(); assert r["passes"], r
    for row in r["results"]:
        print(f"M14 image-plane LOS={row['los_axis']}: "
              f"comp_1={row['comp_1_label']}, comp_2={row['comp_2_label']}")
    r = _antisymmetric_depth_test(); assert r["passes"], r
    print("M14 antisymmetric depth: pass")
    r = _single_slice_test(); assert r["passes"], r
    print("M14 single slice: pass")
    r = _zero_field_test(); assert r["passes"], r
    print("M14 zero field: pass")
    r = _production_vs_reference_test(); assert r["passes"], r
    print(f"M14 production vs reference: max_diff={r['max_diff']:.3e}")
    r = _known_cancellation_test(); assert r["passes"], r
    print("M14 known cancellation: pass")
    print("M14 LOS projection: all checks passed")
