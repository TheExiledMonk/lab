"""Target-blind coordinate and spin-2 basis utilities for Dev Doc 116."""
from __future__ import annotations

import numpy as np


def coordinate_inventory() -> dict:
    return {
        "launch horizontal axis": "launch.x0 (global +x)",
        "launch vertical axis": "launch.y0 (global +y)",
        "observer u axis": "e1: projection of global +x onto detector plane",
        "observer v axis": "e2=cross(screen_normal,e1)",
        "array row axis": "observer v; increasing row is increasing v",
        "array column axis": "observer u; increasing column is increasing u",
        "received x/y/z convention": "global Cartesian rx,ry,rz; dx,dy,dz are e1,e2,normal projections",
        "FITS axis ordering": "numpy [row=NAXIS2, column=NAXIS1]; WCS axis 1=RA, axis 2=DEC",
        "reference gamma1 array ordering": "FITS primary array [DEC pixel, RA pixel]",
        "reference gamma2 array ordering": "FITS primary array [DEC pixel, RA pixel]",
        "screen handedness": "right-handed because e2=normal cross e1 and e1 cross e2=normal",
        "image origin convention": "FITS is 1-based; numpy/deposited maps are 0-based lower-index origin",
        "increasing-row physical direction": "increasing v; FITS WCS increasing DEC",
        "increasing-column physical direction": "increasing u; FITS WCS decreasing RA (east is left)",
    }


def benchmark_basis_from_headers(headers: dict[str, object]) -> dict:
    """Inventory WCS orientation without examining target array values."""
    mats = {}
    for cluster, h in headers.items():
        cd = np.array([[float(h.get("CD1_1", h.get("CDELT1"))), float(h.get("CD1_2", 0.0))],
                       [float(h.get("CD2_1", 0.0)), float(h.get("CD2_2", h.get("CDELT2")))]])
        scale = np.linalg.norm(cd, axis=0)
        mats[cluster] = {"pixel_to_world_linear": cd.tolist(), "orthogonal_orientation": (cd / scale).tolist(),
                         "determinant": float(np.linalg.det(cd / scale)),
                         "ctype": [str(h.get("CTYPE1", "")), str(h.get("CTYPE2", ""))]}
    # WCS locates pixels, but these files do not declare whether gamma is in an
    # east/north, RA/DEC, or pixel-column/row component convention.
    return {"per_cluster": mats, "unique_transform": None,
            "reason": "FITS WCS specifies pixel orientation but no gamma-component basis convention",
            "target_data_used": False}


def spin2_matrix(R) -> np.ndarray:
    """Derive component map by S' = R S R^T (works for reflections)."""
    R = np.asarray(R, float)
    cols = []
    for g1, g2 in ((1.0, 0.0), (0.0, 1.0)):
        S = np.array([[g1, g2], [g2, -g1]])
        T = R @ S @ R.T
        cols.append([0.5 * (T[0, 0] - T[1, 1]), 0.5 * (T[0, 1] + T[1, 0])])
    return np.asarray(cols).T


def apply_spin2(gamma1, gamma2, R):
    T = spin2_matrix(R)
    return T[0, 0]*gamma1 + T[0, 1]*gamma2, T[1, 0]*gamma1 + T[1, 1]*gamma2

