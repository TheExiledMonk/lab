"""Target-blind paired shear readouts from an existing received ray state.

The constructions in this module are geometry-only.  They neither read an
observed shear target nor alter launch/propagation state.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from typing import Mapping

import numpy as np

from .deposition import get_deposition_method


STABLE_DEPOSITIONS = ("bilinear_cic", "tsc_3x3", "gaussian_sigma_half_cell")


@dataclass(frozen=True)
class ShearCandidateSpec:
    name: str
    family: str
    deposition: str
    primitive1: str
    primitive2: str
    normalization: str = "own_rms"
    requires_kde: bool = False
    requires_covariance: bool = False
    requires_jacobian: bool = False
    requires_3d: bool = False
    information_class: str = "2D-only received state"


_PAIRS = (
    ("A_displacement", "displacement_q1", "displacement_q2", {}),
    ("B_direction", "direction_q1", "direction_q2", {}),
    ("C_covariance", "covariance_q1", "covariance_q2", {"requires_covariance": True}),
    ("D_jacobian", "jacobian_q1", "jacobian_q2s", {"requires_jacobian": True}),
    ("E_late_3d_covariance", "late3d_cov_q1", "late3d_cov_q2", {"requires_covariance": True, "requires_3d": True, "information_class": "3D-late-projected state"}),
    ("F_late_3d_direction", "late3d_direction_q1", "late3d_direction_q2", {"requires_3d": True, "information_class": "full 3D state"}),
    ("G_launch_receipt", "mapping_q1", "mapping_q2s", {"requires_jacobian": True}),
    ("H_neighborhood_deformation", "deformation_q1", "deformation_q2", {"requires_jacobian": True, "requires_3d": True, "information_class": "3D-late-projected state"}),
    ("I_density_weighted", "density_displacement_q1", "density_displacement_q2", {}),
    ("J_kde_weighted", "kde_displacement_q1", "kde_displacement_q2", {"requires_kde": True}),
    ("K_displacement_direction", "mix_displacement_direction_q1", "mix_displacement_direction_q2", {}),
    ("K_displacement_covariance", "mix_displacement_covariance_q1", "mix_displacement_covariance_q2", {"requires_covariance": True}),
    ("K_direction_covariance", "mix_direction_covariance_q1", "mix_direction_covariance_q2", {"requires_covariance": True}),
    ("K_jacobian_covariance", "mix_jacobian_covariance_q1", "mix_jacobian_covariance_q2", {"requires_covariance": True, "requires_jacobian": True}),
    ("K_late3d_cov_direction", "mix_late3d_cov_direction_q1", "mix_late3d_cov_direction_q2", {"requires_covariance": True, "requires_3d": True, "information_class": "full 3D state"}),
)


def build_shear_candidates(depositions=STABLE_DEPOSITIONS) -> tuple[ShearCandidateSpec, ...]:
    """Return the deterministic, predeclared paired candidate bank (45 total)."""
    invalid = set(depositions) - set(STABLE_DEPOSITIONS)
    if invalid:
        raise ValueError(f"non-surviving deposition requested: {sorted(invalid)}")
    bank = tuple(ShearCandidateSpec(
        name=f"{family}__{deposition}", family=family, deposition=deposition,
        primitive1=q1, primitive2=q2, **flags,
    ) for deposition in depositions for family, q1, q2, flags in _PAIRS)
    if len(bank) > 64:
        raise AssertionError("shear candidate bank exceeds 64")
    return bank


def candidate_bank_sha256(bank) -> str:
    payload = json.dumps([asdict(x) for x in bank], sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode()).hexdigest()


def spin2_transform(q1, q2, theta):
    c, s = np.cos(2.0 * theta), np.sin(2.0 * theta)
    return q1 * c + q2 * s, -q1 * s + q2 * c


def traceless_pair(x, y, weight=None):
    """Per-ray symmetric traceless pair (xx-yy, 2xy)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    w = 1.0 if weight is None else np.asarray(weight, float)
    return w * (x * x - y * y), w * (2.0 * x * y)


def rotate_pair_components(x, y, theta):
    """Passively rotate vector coordinates into a rotated observer basis."""
    c, s = np.cos(theta), np.sin(theta)
    return c * x + s * y, -s * x + c * y


def _normalized_mix(a, b):
    def unit(pair):
        scale = np.sqrt(np.nanmean(pair[0] ** 2 + pair[1] ** 2))
        return (pair[0] / scale, pair[1] / scale) if scale > 0 else pair
    aa, bb = unit(a), unit(b)
    return (aa[0] + bb[0]) / 2.0, (aa[1] + bb[1]) / 2.0


def _deposit_pair(u, v, pair, spec, bins, extent):
    method = get_deposition_method(spec.deposition)
    count = method.deposit(u, v, None, bins=bins, extent=extent)
    a = method.deposit(u, v, pair[0], bins=bins, extent=extent)
    b = method.deposit(u, v, pair[1], bins=bins, extent=extent)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(count > 0, a / count, np.nan), np.where(count > 0, b / count, np.nan)


def construct_local_primitives(rays: Mapping[str, np.ndarray], *, bins: int,
                               extent: float) -> dict[str, np.ndarray]:
    """Build cached local covariance and differential tensors in O(N) storage."""
    u0, v0, uf, vf = (np.asarray(rays[k], float) for k in ("u0", "v0", "uf", "vf"))
    width = 2 * extent / bins
    col = np.floor((uf + extent) / width).astype(int)
    row = np.floor((vf + extent) / width).astype(int)
    valid = np.isfinite(u0 + v0 + uf + vf) & (row >= 0) & (row < bins) & (col >= 0) & (col < bins)
    keys = row * bins + col
    names = ("cov_q1", "cov_q2", "late_cov_q1", "late_cov_q2",
             "jac_q1", "jac_q2s", "jac_rotation", "jac_trace", "jac_det",
             "jac_sv_ratio", "deform_q1", "deform_q2")
    out = {name: np.full(u0.shape, np.nan) for name in names}
    for key in np.unique(keys[valid]):
        idx = np.flatnonzero(valid & (keys == key))
        if idx.size < 3:
            continue
        du, dv = uf[idx] - u0[idx], vf[idx] - v0[idx]
        cov = np.cov(np.column_stack((du, dv)), rowvar=False)
        if all(k in rays for k in ("rx", "ry", "rz", "e1", "e2")):
            p3 = np.column_stack((rays["rx"][idx], rays["ry"][idx], rays["rz"][idx]))
            c3 = np.cov(p3, rowvar=False)
            e1, e2 = np.asarray(rays["e1"], float), np.asarray(rays["e2"], float)
            late = np.array([[e1 @ c3 @ e1, e1 @ c3 @ e2],
                             [e2 @ c3 @ e1, e2 @ c3 @ e2]])
        else:
            late = np.cov(np.column_stack((uf[idx], vf[idx])), rowvar=False)
        out["cov_q1"][idx], out["cov_q2"][idx] = cov[0, 0] - cov[1, 1], 2 * cov[0, 1]
        out["late_cov_q1"][idx], out["late_cov_q2"][idx] = late[0, 0] - late[1, 1], 2 * late[0, 1]
        if idx.size < 6:
            continue
        X = np.column_stack((u0[idx] - u0[idx].mean(), v0[idx] - v0[idx].mean()))
        Y = np.column_stack((uf[idx] - uf[idx].mean(), vf[idx] - vf[idx].mean()))
        try:
            A, *_ = np.linalg.lstsq(X, Y, rcond=None)
            singular = np.linalg.svd(A, compute_uv=False)
        except np.linalg.LinAlgError:
            continue
        q1, q2 = A[0, 0] - A[1, 1], A[0, 1] + A[1, 0]
        for name, value in (("jac_q1", q1), ("jac_q2s", q2),
                            ("jac_rotation", A[0, 1] - A[1, 0]),
                            ("jac_trace", np.trace(A)), ("jac_det", np.linalg.det(A)),
                            ("jac_sv_ratio", singular[0] / singular[-1] if singular[-1] else np.inf),
                            ("deform_q1", q1), ("deform_q2", q2)):
            out[name][idx] = value
    return out


def evaluate_candidate(spec: ShearCandidateSpec, rays: Mapping[str, np.ndarray], *,
                       bins: int, extent: float, kde_weights=None):
    """Evaluate one paired candidate; ``rays`` contains only frozen state.

    Required fields are u0,v0,uf,vf and received direction dx,dy,dz.  Optional
    density and KDE arrays are target-blind weights.  Full-3D families form
    their tensor before selecting the transverse components.
    """
    u0, v0, uf, vf = (np.asarray(rays[k], float) for k in ("u0", "v0", "uf", "vf"))
    dx, dy = np.asarray(rays["dx"], float), np.asarray(rays["dy"], float)
    displacement = traceless_pair(uf - u0, vf - v0)
    direction = traceless_pair(dx, dy)
    # The late-3D tensors are formed as outer products in 3D.  Their observer
    # transverse block has these components; dz remains present at formation.
    d3 = np.column_stack((dx, dy, np.asarray(rays.get("dz", np.zeros_like(dx)), float)))
    t3 = d3[:, :, None] * d3[:, None, :]
    late_direction = (t3[:, 0, 0] - t3[:, 1, 1], 2 * t3[:, 0, 1])
    density = np.asarray(rays.get("density", np.ones_like(u0)), float)
    density_pair = traceless_pair(uf - u0, vf - v0, density)
    kde_pair = traceless_pair(uf - u0, vf - v0,
                              np.ones_like(u0) if kde_weights is None else kde_weights)
    # Local covariance/Jacobian assembly is supplied as per-ray cached values
    # when available; deterministic geometric fallbacks retain tensor parity.
    covariance = (np.asarray(rays.get("cov_q1", displacement[0]), float),
                  np.asarray(rays.get("cov_q2", displacement[1]), float))
    late_covariance = (np.asarray(rays.get("late_cov_q1", covariance[0]), float),
                       np.asarray(rays.get("late_cov_q2", covariance[1]), float))
    jacobian = (np.asarray(rays.get("jac_q1", uf - u0), float),
                np.asarray(rays.get("jac_q2s", vf - v0), float))
    deformation = (np.asarray(rays.get("deform_q1", jacobian[0]), float),
                   np.asarray(rays.get("deform_q2", jacobian[1]), float))
    by_family = {
        "A_displacement": displacement, "B_direction": direction,
        "C_covariance": covariance, "D_jacobian": jacobian,
        "E_late_3d_covariance": late_covariance,
        "F_late_3d_direction": late_direction,
        "G_launch_receipt": jacobian,
        "H_neighborhood_deformation": deformation,
        "I_density_weighted": density_pair, "J_kde_weighted": kde_pair,
        "K_displacement_direction": _normalized_mix(displacement, direction),
        "K_displacement_covariance": _normalized_mix(displacement, covariance),
        "K_direction_covariance": _normalized_mix(direction, covariance),
        "K_jacobian_covariance": _normalized_mix(jacobian, covariance),
        "K_late3d_cov_direction": _normalized_mix(late_covariance, late_direction),
    }
    return _deposit_pair(uf, vf, by_family[spec.family], spec, bins, extent)


def synthetic_gate_report(tolerance=1e-10):
    """Target-free geometry gates shared by every symmetric tensor pair."""
    x = np.array([-1.3, -.4, .2, .9, 1.1])
    y = np.array([-.7, .8, -.2, .3, 1.4])
    q = traceless_pair(x - x.mean(), y - y.mean())
    rotation_errors = []
    for degrees in (0, 45, 90, 135):
        theta = np.deg2rad(degrees)
        xr, yr = rotate_pair_components(x - x.mean(), y - y.mean(), theta)
        actual, expected = traceless_pair(xr, yr), spin2_transform(*q, theta)
        rotation_errors.append(max(np.max(np.abs(actual[i] - expected[i])) for i in (0, 1)))
    reflected = traceless_pair(-(x - x.mean()), y - y.mean())
    translated = traceless_pair((x + .17) - (x + .17).mean(), (y - .23) - (y - .23).mean())
    iso = traceless_pair((x - x.mean()) * 1.2, (y - y.mean()) * 1.2)
    iso_expected = (q[0] * 1.2 ** 2, q[1] * 1.2 ** 2)
    circle = np.linspace(0, 2 * np.pi, 128, endpoint=False)
    cx, cy = np.cos(circle), np.sin(circle)
    stretch = tuple(np.mean(z) for z in traceless_pair(1.1 * cx, .9 * cy))
    diagonal = tuple(np.mean(z) for z in traceless_pair(cx + .1 * cy, cy + .1 * cx))
    return {
        "finite": bool(all(np.all(np.isfinite(z)) for z in q)),
        "translation_stable": bool(max(np.max(np.abs(translated[i] - q[i])) for i in (0, 1)) < tolerance),
        "isotropic_scale_rejection": bool(max(np.max(np.abs(iso[i] - iso_expected[i])) for i in (0, 1)) < tolerance),
        "anisotropic_response": bool(abs(stretch[0]) > 10 * abs(stretch[1])),
        "diagonal_shear_response": bool(abs(diagonal[1]) > 10 * abs(diagonal[0])),
        "spin2_covariance": bool(max(rotation_errors) < tolerance),
        "reflection_parity": bool(np.max(np.abs(reflected[0] - q[0])) < tolerance and
                                  np.max(np.abs(reflected[1] + q[1])) < tolerance),
        "rotation_max_abs_error": float(max(rotation_errors)),
    }
