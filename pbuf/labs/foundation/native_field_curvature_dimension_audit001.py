#!/usr/bin/env python3
"""PBUF FOUNDATION — NATIVE FIELD CURVATURE DIMENSION AUDIT 001.

Fact-finding only.

Purpose
-------
Identify how the existing native PBUF fields scale with source density and
source size before assigning them a physical role in the mass -> spacetime
response bridge.

The previous mass-spacetime-response audit established the macroscopic anchor
but deliberately left the native PBUF mapping open.  This lab does NOT insert
G, h00, Newtonian potential, kappa, shear, HST data, or the historical 0.18
coefficient.  Instead it asks a narrower structural question:

    Under the frozen PBUF medium dynamics, which native field behaves like a
    local source/curvature quantity, a one-length integrated response, or a
    two-length integrated dimensionless deformation?

Synthetic experiment
--------------------
Use centered uniform 3D spheres on a fixed grid and run two independent ladders:

1. DENSITY ladder: fixed sphere radius, varying source density amplitude.
   A linear physical response should have log-log slope ~= 1 with density.

2. RADIUS ladder: fixed density, varying sphere radius.
   Compare the measured radius exponent against predeclared structural classes:

       R^0 : local source / curvature-like scaling
       R^1 : one-length integrated / gradient-connection-like scaling
       R^2 : two-length integrated / strain-metric-deformation-like scaling

These labels are dimensional/scaling diagnostics only.  They do NOT prove that
any numerical field already carries SI curvature, connection, strain, or metric
units.  The frozen transport equations contain their own grid scale, so the
result is a classification of the current implementation's response behavior.

Native quantities audited
--------------------------
- c_state scalar field
- |grad c_state|
- |laplacian c_state|
- M10 interface-vector magnitude
- |div M10|
- Frobenius norm of symmetrized spatial gradient of M10

For scalar/local fields we report center/interior/shell statistics where
appropriate.  For derivative/vector fields we use shell statistics around the
source boundary so spherical cancellation at the exact center does not erase
an otherwise real response.

Hard guardrails
---------------
- no G anywhere in the calculation;
- no lensing target, kappa, shear, HST, morphology, or observer comparison;
- no legacy strength=0.18 input;
- no fitting, tuning, optimization, or coefficient solving;
- no Quantum Engine or Planck-scale input;
- no random injection noise;
- no clipping allowed during the accepted runs;
- no production code modifications;
- stdout only; no run directory;
- classification thresholds and ladders are frozen before results are computed.
"""
from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from pbuf.core import los_projection as M14
from pbuf.models import a8_state as M06_state
import pbuf.labs.foundation.m10_coverage_25pct_science001 as BASE

LAB_ID = "PBUF-FOUNDATION-NATIVE-FIELD-CURVATURE-DIMENSION-AUDIT-001"

# Fixed synthetic geometry.  Coordinates are grid units; DX is deliberately
# explicit so derivative order is transparent.
SHAPE = (33, 33, 33)  # z, y, x
DX = 1.0
CENTER_INDEX = tuple(n // 2 for n in SHAPE)

# Ladders are predeclared and do not depend on any observed gravitational value.
FIXED_DENSITY = 0.05
RADIUS_LADDER = (2.5, 3.5, 4.5, 5.5, 6.5)
FIXED_RADIUS = 4.5
DENSITY_LADDER = (0.0125, 0.025, 0.05, 0.10, 0.20)

# Shell around the source boundary used for derivative/vector observables.
SHELL_INNER = 0.75
SHELL_OUTER = 1.25

# Scaling-class diagnostics.  These are not fitted physical coefficients.
RADIUS_CLASSES = {
    "R0_LOCAL_SOURCE_CURVATURE_LIKE": 0.0,
    "R1_ONE_LENGTH_INTEGRATED_GRADIENT_CONNECTION_LIKE": 1.0,
    "R2_TWO_LENGTH_INTEGRATED_STRAIN_METRIC_LIKE": 2.0,
}
CLASS_WARN_DISTANCE = 0.40
LINEAR_DENSITY_WARN_DISTANCE = 0.15
ALG_TOL = 1.0e-12


def _git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=str(ROOT), text=True).strip()


def _repo_state() -> dict:
    return {
        "repository": "TheExiledMonk/lab",
        "branch": _git("rev-parse", "--abbrev-ref", "HEAD"),
        "head_sha": _git("rev-parse", "HEAD"),
        "tracked_changes": _git("diff", "--name-only"),
        "staged_changes": _git("diff", "--name-only", "--cached"),
    }


def _rms(x) -> float:
    a = np.asarray(x, dtype=np.float64)
    a = a[np.isfinite(a)]
    return float(np.sqrt(np.mean(a * a))) if a.size else float("nan")


def _mean_abs(x) -> float:
    a = np.asarray(x, dtype=np.float64)
    a = a[np.isfinite(a)]
    return float(np.mean(np.abs(a))) if a.size else float("nan")


def _coords():
    z = (np.arange(SHAPE[0], dtype=np.float64) - CENTER_INDEX[0]) * DX
    y = (np.arange(SHAPE[1], dtype=np.float64) - CENTER_INDEX[1]) * DX
    x = (np.arange(SHAPE[2], dtype=np.float64) - CENTER_INDEX[2]) * DX
    Z, Y, X = np.meshgrid(z, y, x, indexing="ij")
    R = np.sqrt(X * X + Y * Y + Z * Z)
    return X, Y, Z, R


XGRID, YGRID, ZGRID, RGRID = _coords()


def _sphere(radius: float, density: float) -> dict:
    mask = RGRID <= float(radius)
    rho = np.zeros(SHAPE, dtype=np.float64)
    rho[mask] = float(density)
    voxel_volume = DX ** 3
    volume = float(np.count_nonzero(mask)) * voxel_volume
    r_eff = (3.0 * volume / (4.0 * math.pi)) ** (1.0 / 3.0)
    return {
        "rho": rho,
        "mask": mask,
        "nominal_radius": float(radius),
        "effective_radius": float(r_eff),
        "density": float(density),
        "volume": volume,
        "integrated_source": float(np.sum(rho) * voxel_volume),
    }


def _shell_mask(r_eff: float) -> np.ndarray:
    return (RGRID >= SHELL_INNER * r_eff) & (RGRID <= SHELL_OUTER * r_eff)


def _laplacian(a: np.ndarray) -> np.ndarray:
    a = np.asarray(a, dtype=np.float64)
    gz, gy, gx = np.gradient(a, DX, DX, DX, edge_order=2)
    gzz = np.gradient(gz, DX, axis=0, edge_order=2)
    gyy = np.gradient(gy, DX, axis=1, edge_order=2)
    gxx = np.gradient(gx, DX, axis=2, edge_order=2)
    return gxx + gyy + gzz


def _gradient_mag(a: np.ndarray) -> np.ndarray:
    gz, gy, gx = np.gradient(np.asarray(a, dtype=np.float64), DX, DX, DX, edge_order=2)
    return np.sqrt(gx * gx + gy * gy + gz * gz)


def _divergence(vx: np.ndarray, vy: np.ndarray, vz: np.ndarray) -> np.ndarray:
    dvx_dx = np.gradient(vx, DX, axis=2, edge_order=2)
    dvy_dy = np.gradient(vy, DX, axis=1, edge_order=2)
    dvz_dz = np.gradient(vz, DX, axis=0, edge_order=2)
    return dvx_dx + dvy_dy + dvz_dz


def _strain_frobenius(vx: np.ndarray, vy: np.ndarray, vz: np.ndarray) -> np.ndarray:
    # Symmetric spatial gradient E_ij = 1/2(d_i v_j + d_j v_i).
    dvx_dz, dvx_dy, dvx_dx = np.gradient(vx, DX, DX, DX, edge_order=2)
    dvy_dz, dvy_dy, dvy_dx = np.gradient(vy, DX, DX, DX, edge_order=2)
    dvz_dz, dvz_dy, dvz_dx = np.gradient(vz, DX, DX, DX, edge_order=2)

    exx = dvx_dx
    eyy = dvy_dy
    ezz = dvz_dz
    exy = 0.5 * (dvx_dy + dvy_dx)
    exz = 0.5 * (dvx_dz + dvz_dx)
    eyz = 0.5 * (dvy_dz + dvz_dy)
    return np.sqrt(exx * exx + eyy * eyy + ezz * ezz + 2.0 * (exy * exy + exz * exz + eyz * eyz))


def _noise_free_state(rho: np.ndarray) -> dict:
    # Unit source-loading diagnostic.  The historical 0.18 coefficient is not
    # called, imported, multiplied, divided, or inferred here.
    u_slow0 = np.asarray(rho, dtype=np.float64).copy()
    u_fast0 = np.asarray(rho, dtype=np.float64).copy()
    us, uf, history = M06_state.evolve_a8_transport_3d(
        u_slow0, u_fast0, stencil="N6", boundary="reflective"
    )
    c = np.asarray(history[-1], dtype=np.float64)
    max_abs = max(float(np.max(np.abs(us))), float(np.max(np.abs(uf))))
    if max_abs >= M06_state.A8_INIT_CLIP - 1.0e-12:
        raise RuntimeError(f"clipping gate failed: max_abs={max_abs}, clip={M06_state.A8_INIT_CLIP}")
    return {
        "rho_3d": np.asarray(rho, dtype=np.float64).copy(),
        "u_slow": np.asarray(us, dtype=np.float64),
        "u_fast": np.asarray(uf, dtype=np.float64),
        "c_state": c,
        "max_abs_state": max_abs,
    }


def _native_metrics(source: dict) -> dict:
    state = _noise_free_state(source["rho"])
    candidate = BASE._candidate(state)
    vx, vy, vz = BASE._interface_vector(candidate)
    vx = np.asarray(vx, dtype=np.float64)
    vy = np.asarray(vy, dtype=np.float64)
    vz = np.asarray(vz, dtype=np.float64)

    c = state["c_state"]
    grad_c = _gradient_mag(c)
    lap_c = _laplacian(c)
    vmag = np.sqrt(vx * vx + vy * vy + vz * vz)
    div_v = _divergence(vx, vy, vz)
    strain = _strain_frobenius(vx, vy, vz)

    mask = source["mask"]
    shell = _shell_mask(source["effective_radius"])
    ci = CENTER_INDEX

    # M14 projection is included only as a representation check, not an
    # observer/lensing calculation.  It confirms the 3D M10 vector is finite.
    los = M14.project_vector_to_image_plane(vx, vy, vz, los_axis="z")
    los_mag = np.hypot(np.asarray(los["comp_1"]), np.asarray(los["comp_2"]))

    return {
        "nominal_radius": source["nominal_radius"],
        "effective_radius": source["effective_radius"],
        "density": source["density"],
        "voxel_volume": DX ** 3,
        "source_volume": source["volume"],
        "integrated_source": source["integrated_source"],
        "source_voxel_count": int(np.count_nonzero(mask)),
        "shell_voxel_count": int(np.count_nonzero(shell)),
        "max_abs_state": state["max_abs_state"],
        "c_state_center_abs": float(abs(c[ci])),
        "c_state_interior_mean_abs": _mean_abs(c[mask]),
        "c_state_shell_mean_abs": _mean_abs(c[shell]),
        "grad_c_shell_mean": _mean_abs(grad_c[shell]),
        "laplacian_c_interior_mean_abs": _mean_abs(lap_c[mask]),
        "laplacian_c_shell_mean_abs": _mean_abs(lap_c[shell]),
        "m10_vector_shell_mean": _mean_abs(vmag[shell]),
        "m10_vector_rms": _rms(vmag),
        "m10_divergence_interior_mean_abs": _mean_abs(div_v[mask]),
        "m10_divergence_shell_mean_abs": _mean_abs(div_v[shell]),
        "m10_symgrad_fro_shell_mean": _mean_abs(strain[shell]),
        "m10_projected_mag_rms_representation_only": _rms(los_mag),
        "candidate_gradient_rms_internal": float(candidate["gradient_rms"]),
        "all_native_arrays_finite": bool(
            np.all(np.isfinite(c)) and np.all(np.isfinite(grad_c)) and
            np.all(np.isfinite(lap_c)) and np.all(np.isfinite(vmag)) and
            np.all(np.isfinite(div_v)) and np.all(np.isfinite(strain))
        ),
    }


AUDIT_METRICS = (
    "c_state_center_abs",
    "c_state_interior_mean_abs",
    "c_state_shell_mean_abs",
    "grad_c_shell_mean",
    "laplacian_c_interior_mean_abs",
    "laplacian_c_shell_mean_abs",
    "m10_vector_shell_mean",
    "m10_vector_rms",
    "m10_divergence_interior_mean_abs",
    "m10_divergence_shell_mean_abs",
    "m10_symgrad_fro_shell_mean",
    "m10_projected_mag_rms_representation_only",
)


def _loglog_fit(xs, ys) -> dict:
    x = np.asarray(xs, dtype=np.float64)
    y = np.asarray(ys, dtype=np.float64)
    m = np.isfinite(x) & np.isfinite(y) & (x > 0.0) & (y > 0.0)
    if int(np.count_nonzero(m)) < 3:
        return {"slope": float("nan"), "r2": float("nan"), "count": int(np.count_nonzero(m))}
    lx = np.log(x[m])
    ly = np.log(y[m])
    A = np.column_stack((lx, np.ones_like(lx)))
    beta, *_ = np.linalg.lstsq(A, ly, rcond=None)
    pred = A @ beta
    ss_res = float(np.sum((ly - pred) ** 2))
    ss_tot = float(np.sum((ly - np.mean(ly)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else float("nan")
    return {"slope": float(beta[0]), "r2": r2, "count": int(np.count_nonzero(m))}


def _classify_radius_slope(slope: float) -> dict:
    if not math.isfinite(slope):
        return {"closest_class": "UNRESOLVED", "distance": float("nan"), "within_predeclared_window": False}
    distances = {name: abs(slope - exponent) for name, exponent in RADIUS_CLASSES.items()}
    name = min(distances, key=distances.get)
    d = float(distances[name])
    return {
        "closest_class": name,
        "class_exponent": RADIUS_CLASSES[name],
        "distance": d,
        "within_predeclared_window": bool(d <= CLASS_WARN_DISTANCE),
    }


def _scaling_summary(radius_rows: list[dict], density_rows: list[dict]) -> dict:
    out = {}
    radius_x = [r["effective_radius"] for r in radius_rows]
    density_x = [r["density"] for r in density_rows]
    for metric in AUDIT_METRICS:
        rf = _loglog_fit(radius_x, [r[metric] for r in radius_rows])
        df = _loglog_fit(density_x, [r[metric] for r in density_rows])
        out[metric] = {
            "radius_scaling": rf,
            "radius_classification_diagnostic_only": _classify_radius_slope(rf["slope"]),
            "density_scaling": df,
            "density_linear_distance_from_1": abs(df["slope"] - 1.0) if math.isfinite(df["slope"]) else float("nan"),
            "density_linear_within_predeclared_window": bool(
                math.isfinite(df["slope"]) and abs(df["slope"] - 1.0) <= LINEAR_DENSITY_WARN_DISTANCE
            ),
        }
    return out


def _source_geometry_checks(radius_rows: list[dict], density_rows: list[dict]) -> dict:
    # The synthetic source itself provides controls: at fixed density, integrated
    # source should scale approximately R^3 after effective-radius correction; at
    # fixed radius it should scale exactly linearly with density.
    rfit = _loglog_fit(
        [r["effective_radius"] for r in radius_rows],
        [r["integrated_source"] for r in radius_rows],
    )
    dfit = _loglog_fit(
        [r["density"] for r in density_rows],
        [r["integrated_source"] for r in density_rows],
    )
    return {
        "integrated_source_radius_slope": rfit,
        "integrated_source_density_slope": dfit,
        "radius_volume_scaling_R3_pass": bool(math.isfinite(rfit["slope"]) and abs(rfit["slope"] - 3.0) <= 2.0e-12),
        "density_linearity_pass": bool(math.isfinite(dfit["slope"]) and abs(dfit["slope"] - 1.0) <= ALG_TOL),
    }


def _json_default(x):
    if isinstance(x, np.floating): return float(x)
    if isinstance(x, np.integer): return int(x)
    if isinstance(x, np.bool_): return bool(x)
    if isinstance(x, Path): return str(x)
    return str(x)


def main() -> None:
    repo = _repo_state()

    radius_rows = [_native_metrics(_sphere(r, FIXED_DENSITY)) for r in RADIUS_LADDER]
    density_rows = [_native_metrics(_sphere(FIXED_RADIUS, d)) for d in DENSITY_LADDER]

    geometry_checks = _source_geometry_checks(radius_rows, density_rows)
    scaling = _scaling_summary(radius_rows, density_rows)

    all_finite = all(r["all_native_arrays_finite"] for r in radius_rows + density_rows)
    no_clip = all(r["max_abs_state"] < M06_state.A8_INIT_CLIP - 1.0e-12 for r in radius_rows + density_rows)
    density_linear_metrics = [
        k for k, v in scaling.items() if v["density_linear_within_predeclared_window"]
    ]

    result = {
        "lab_id": LAB_ID,
        "status": "FACT_FINDING_ONLY",
        "repo_state": repo,
        "policy": {
            "gravity_fundamental_in_PBUF": False,
            "G_used": False,
            "legacy_0p18_used": False,
            "lensing_target_used": False,
            "observer_comparison_used": False,
            "random_noise_used": False,
            "fit_or_tuning_used": False,
            "classification_role": "scaling_behavior_diagnostic_only_not_literal_SI_dimension_proof",
        },
        "synthetic_design": {
            "shape_zyx": SHAPE,
            "dx_grid_units": DX,
            "fixed_density_radius_ladder": FIXED_DENSITY,
            "radius_ladder_nominal_grid_units": RADIUS_LADDER,
            "fixed_radius_density_ladder": FIXED_RADIUS,
            "density_ladder_dimensionless_loading_units": DENSITY_LADDER,
            "noise_free_initialization": "u_slow0=u_fast0=rho3",
            "shell_fraction": [SHELL_INNER, SHELL_OUTER],
            "predeclared_radius_classes": RADIUS_CLASSES,
        },
        "radius_rows": radius_rows,
        "density_rows": density_rows,
        "source_geometry_controls": geometry_checks,
        "scaling": scaling,
        "summary": {
            "density_linear_metric_count": len(density_linear_metrics),
            "density_linear_metrics": density_linear_metrics,
            "interpretation_guardrail": (
                "A radius class identifies the scaling rank of the current fixed-grid implementation. "
                "It does not by itself establish physical SI units or close the mass-to-metric map."
            ),
            "next_if_clean": (
                "Use the identified native scaling rank to construct the dimensional mapping from the empirical "
                "mass-induced curvature source onto the native PBUF variable without importing 0.18 or fitting lensing."
            ),
        },
        "checks": {
            "all_native_arrays_finite": all_finite,
            "no_state_clipping": no_clip,
            "synthetic_source_radius_volume_control_R3_pass": geometry_checks["radius_volume_scaling_R3_pass"],
            "synthetic_source_density_linearity_control_pass": geometry_checks["density_linearity_pass"],
            "G_not_used": True,
            "legacy_0p18_not_used": True,
            "no_lensing_target_input": True,
            "no_observer_comparison": True,
            "no_random_injection_noise": True,
            "no_fit_or_tuning": True,
            "no_quantum_engine_input": True,
            "no_planck_scale_input": True,
            "no_tracked_or_staged_changes_created_by_lab": repo["tracked_changes"] == "" and repo["staged_changes"] == "",
            "stdout_only_no_run_directory_created": True,
        },
    }

    print(LAB_ID)
    print("status=FACT_FINDING_ONLY")
    print(f"head_sha={repo['head_sha']}")
    print("G_used=false")
    print("gravity_fundamental_in_PBUF=false")
    print("legacy_0p18_used=false")
    print("lensing_target_used=false")
    print("fit_or_tuning_used=false")
    print("random_noise_used=false")
    print()

    print("RADIUS_LADDER_NATIVE_METRICS")
    print("R_eff | density | c_center | c_inside | grad_c_shell | lap_c_inside | M10_shell | div_M10_inside | symgrad_M10_shell")
    for r in radius_rows:
        print(
            f"{r['effective_radius']:.17e} | {r['density']:.17e} | "
            f"{r['c_state_center_abs']:.17e} | {r['c_state_interior_mean_abs']:.17e} | "
            f"{r['grad_c_shell_mean']:.17e} | {r['laplacian_c_interior_mean_abs']:.17e} | "
            f"{r['m10_vector_shell_mean']:.17e} | {r['m10_divergence_interior_mean_abs']:.17e} | "
            f"{r['m10_symgrad_fro_shell_mean']:.17e}"
        )
    print()

    print("DENSITY_LADDER_NATIVE_METRICS")
    print("R_eff | density | c_center | c_inside | grad_c_shell | lap_c_inside | M10_shell | div_M10_inside | symgrad_M10_shell")
    for r in density_rows:
        print(
            f"{r['effective_radius']:.17e} | {r['density']:.17e} | "
            f"{r['c_state_center_abs']:.17e} | {r['c_state_interior_mean_abs']:.17e} | "
            f"{r['grad_c_shell_mean']:.17e} | {r['laplacian_c_interior_mean_abs']:.17e} | "
            f"{r['m10_vector_shell_mean']:.17e} | {r['m10_divergence_interior_mean_abs']:.17e} | "
            f"{r['m10_symgrad_fro_shell_mean']:.17e}"
        )
    print()

    print("SCALING_CLASSIFICATION")
    print("metric | density_slope | density_R2 | radius_slope | radius_R2 | closest_radius_class | class_distance | class_window")
    for metric in AUDIT_METRICS:
        s = scaling[metric]
        d = s["density_scaling"]
        r = s["radius_scaling"]
        c = s["radius_classification_diagnostic_only"]
        print(
            f"{metric} | {d['slope']:.17e} | {d['r2']:.17e} | "
            f"{r['slope']:.17e} | {r['r2']:.17e} | {c['closest_class']} | "
            f"{c['distance']:.17e} | {c['within_predeclared_window']}"
        )
    print()

    print("SOURCE_GEOMETRY_CONTROLS")
    print(f"integrated_source_radius_slope={geometry_checks['integrated_source_radius_slope']['slope']:.17e}")
    print(f"integrated_source_density_slope={geometry_checks['integrated_source_density_slope']['slope']:.17e}")
    print(f"radius_volume_scaling_R3_pass={str(geometry_checks['radius_volume_scaling_R3_pass']).lower()}")
    print(f"density_linearity_pass={str(geometry_checks['density_linearity_pass']).lower()}")
    print()

    print("CHECKS")
    for k, v in result["checks"].items():
        print(f"{k}={str(v).lower() if isinstance(v, bool) else v}")
    print("JSON=" + json.dumps(result, sort_keys=True, separators=(",", ":"), default=_json_default))

    if not all(result["checks"].values()):
        raise SystemExit(2)


if __name__ == "__main__":
    main()
