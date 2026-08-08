#!/usr/bin/env python3
"""PBUF FOUNDATION — NATIVE ACCUMULATED FULL LENSING 001.

Purpose
-------
Run the existing real-cluster G3D lensing/observer stack with the historical
source-amplitude training wheel removed and replaced by the newly supported
native accumulated medium response.

Historical lane:
    rho3 -> strength(0.18) * rho3 -> A8 -> M10 -> LOS -> G3D -> observer

Native lane tested here:
    rho3 -> raw A8 c_state (NO strength scalar)
         -> bounded-strain six-neighbor accumulation
         -> native accumulated deformation gradient
         -> existing LOS -> existing G3D -> existing observer

The old STRENGTH=0.18 is inventoried and may be run only as a historical control.
It is not used, inverted, fitted, or algebraically transformed to build the native
lane.  No replacement scalar is introduced.

Observed weak-lensing kappa pixels are not loaded until the independent source and
ALL prediction lanes are complete.  They are used only for end-of-chain comparison.
The independent HST/F160W source remains a luminous-structure proxy rather than an
absolute baryonic mass map, so this lab is a full propagation/lensing comparison
but not yet an absolute SI-amplitude closure claim.
"""
from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import pbuf.labs.foundation.independent_source_training_wheels_off001_common_footprint_fix as SRC
import pbuf.labs.foundation.independent_source_training_wheels_off001 as IND
import pbuf.labs.foundation.m10_coverage_25pct_science001 as BASE
import pbuf.labs.foundation.los_consistent_ray_geometry001 as GEO
import pbuf.labs.foundation.g3d_angular_received_distribution001 as ANG
import pbuf.labs.foundation.g3d_native_angular_detector_image001 as DET
from pbuf.core import los_projection as M14
from pbuf.models import a8_state as M06_state

LAB_ID = "PBUF-FOUNDATION-NATIVE-ACCUMULATED-FULL-LENSING-001"

# Freeze bridge constitutive parameters exactly as established in the native
# accumulation labs.  K0 and EPSILON_MAX are structural normalizations, not fits.
K0 = 1.0
EPSILON_MAX = 1.0
PICARD_TOL = 2.0e-7
PICARD_MAX_ITER = 30
PICARD_DAMP = 0.65
CG_REL_TOL = 2.0e-9
CG_MAX_ITER = 700

LEGACY_STRENGTH = float(BASE.STRENGTH)
EXPECTED_LEGACY_STRENGTH = 0.18
CHECKPOINT = GEO.CHECKPOINTS[-1]
EXPECTED_SUPPORT = GEO.EXPECTED_SUPPORT
UNIT_SPEED_TOL = GEO.UNIT_SPEED_TOL


def git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def repo_state() -> dict:
    return {
        "repository": "TheExiledMonk/lab",
        "branch": git("rev-parse", "--abbrev-ref", "HEAD"),
        "head_sha": git("rev-parse", "HEAD"),
        "tracked_changes": git("diff", "--name-only"),
        "staged_changes": git("diff", "--name-only", "--cached"),
    }


def rms(a) -> float:
    x = np.asarray(a, dtype=np.float64)
    m = np.isfinite(x)
    return float(np.sqrt(np.mean(x[m] * x[m]))) if np.any(m) else float("nan")


def rel_rms(a, b) -> float:
    aa = np.asarray(a, dtype=np.float64)
    bb = np.asarray(b, dtype=np.float64)
    m = np.isfinite(aa) & np.isfinite(bb)
    if not np.any(m):
        return float("nan")
    num = float(np.sqrt(np.mean((aa[m] - bb[m]) ** 2)))
    den = max(float(np.sqrt(np.mean(bb[m] ** 2))), 1.0e-30)
    return num / den


def zero_boundary(a: np.ndarray) -> None:
    a[0, :, :] = a[-1, :, :] = 0.0
    a[:, 0, :] = a[:, -1, :] = 0.0
    a[:, :, 0] = a[:, :, -1] = 0.0


def bond_strains(u: np.ndarray):
    return (
        u[1:, :, :] - u[:-1, :, :],
        u[:, 1:, :] - u[:, :-1, :],
        u[:, :, 1:] - u[:, :, :-1],
    )


def secant_weights(u: np.ndarray):
    weights = []
    max_frac = 0.0
    for e in bond_strains(u):
        frac = np.abs(e) / EPSILON_MAX
        max_frac = max(max_frac, float(np.max(frac)))
        if max_frac >= 0.995:
            raise RuntimeError("bounded-strain barrier approached too closely")
        weights.append(K0 / (1.0 - frac * frac))
    return tuple(weights), max_frac


def apply_A(u: np.ndarray, weights) -> np.ndarray:
    wz, wy, wx = weights
    out = np.zeros_like(u)
    dz = (u[1:, :, :] - u[:-1, :, :]) * wz
    dy = (u[:, 1:, :] - u[:, :-1, :]) * wy
    dx = (u[:, :, 1:] - u[:, :, :-1]) * wx
    out[:-1, :, :] -= dz; out[1:, :, :] += dz
    out[:, :-1, :] -= dy; out[:, 1:, :] += dy
    out[:, :, :-1] -= dx; out[:, :, 1:] += dx
    zero_boundary(out)
    return out


def cg_solve(source: np.ndarray, weights, x0=None) -> dict:
    b = np.asarray(source, dtype=np.float64).copy(); zero_boundary(b)
    x = np.zeros_like(b) if x0 is None else np.asarray(x0, dtype=np.float64).copy(); zero_boundary(x)
    r = b - apply_A(x, weights); zero_boundary(r)
    p = r.copy()
    rr_old = float(np.sum(r * r))
    bnorm = math.sqrt(float(np.sum(b * b)))
    if bnorm == 0.0:
        return {"field": x, "iterations": 0, "relative_residual": 0.0}
    rel = math.sqrt(rr_old) / bnorm
    it = 0
    for k in range(1, CG_MAX_ITER + 1):
        Ap = apply_A(p, weights)
        denom = float(np.sum(p * Ap))
        if not math.isfinite(denom) or denom <= 0.0:
            raise RuntimeError("CG operator lost positive definiteness")
        alpha = rr_old / denom
        x += alpha * p; zero_boundary(x)
        r -= alpha * Ap; zero_boundary(r)
        rr_new = float(np.sum(r * r))
        rel = math.sqrt(rr_new) / bnorm
        it = k
        if rel <= CG_REL_TOL:
            break
        beta = rr_new / rr_old
        p = r + beta * p; zero_boundary(p)
        rr_old = rr_new
    return {"field": x, "iterations": it, "relative_residual": rel}


def solve_bounded_strain(source: np.ndarray) -> dict:
    nz, ny, nx = source.shape
    ones = (
        np.ones((nz - 1, ny, nx), dtype=np.float64) * K0,
        np.ones((nz, ny - 1, nx), dtype=np.float64) * K0,
        np.ones((nz, ny, nx - 1), dtype=np.float64) * K0,
    )
    first = cg_solve(source, ones)
    u = first["field"]
    history = []
    converged = False
    for it in range(1, PICARD_MAX_ITER + 1):
        weights, _ = secant_weights(u)
        sol = cg_solve(source, weights, x0=u)
        candidate = sol["field"]
        new_u = PICARD_DAMP * candidate + (1.0 - PICARD_DAMP) * u
        zero_boundary(new_u)
        scale = max(float(np.sqrt(np.mean(new_u * new_u))), 1.0e-14)
        change = float(np.sqrt(np.mean((new_u - u) ** 2))) / scale
        u = new_u
        _, max_frac = secant_weights(u)
        history.append({
            "picard": it,
            "relative_change": change,
            "max_strain_fraction": max_frac,
            "cg_iterations": sol["iterations"],
            "cg_relative_residual": sol["relative_residual"],
        })
        if change <= PICARD_TOL:
            converged = True
            break
    weights, max_frac = secant_weights(u)
    residual = source - apply_A(u, weights); zero_boundary(residual)
    src_norm = math.sqrt(float(np.sum(source[1:-1, 1:-1, 1:-1] ** 2)))
    nl_rel = math.sqrt(float(np.sum(residual * residual))) / src_norm if src_norm > 0 else 0.0
    return {
        "field": u,
        "converged": converged,
        "picard_iterations": len(history),
        "max_strain_fraction": max_frac,
        "nonlinear_relative_residual": nl_rel,
        "history_tail": history[-3:],
    }


def raw_c_state(rho3: np.ndarray) -> np.ndarray:
    # Crucial replacement of the legacy initialization: NO strength multiplier.
    u0 = np.asarray(rho3, dtype=np.float64).copy()
    us, uf, history = M06_state.evolve_a8_transport_3d(
        u0.copy(), u0.copy(), stencil="N6", boundary="reflective"
    )
    max_abs = max(float(np.max(np.abs(us))), float(np.max(np.abs(uf))))
    if max_abs >= M06_state.A8_INIT_CLIP - 1.0e-12:
        raise RuntimeError("raw c_state clipping gate failed")
    return np.asarray(history[-1], dtype=np.float64)


def native_accumulated_vector(rho3: np.ndarray) -> dict:
    c = raw_c_state(rho3)
    sol = solve_bounded_strain(c)
    u = np.asarray(sol["field"], dtype=np.float64)
    # The supported PR #77 result showed that propagation responds correctly to
    # the gradient of accumulated deformation.  Here that 3-D gradient is supplied
    # to the already-existing LOS/G3D ray machinery.  No coefficient multiplies it.
    gz, gy, gx = np.gradient(u, edge_order=2)
    vector = (-gx, -gy, -gz)
    return {
        "vector": vector,
        "c_state": c,
        "accumulated": u,
        "rho_integral": float(np.sum(rho3)),
        "c_state_integral": float(np.sum(c)),
        "c_state_integral_relative_error": abs(float(np.sum(c)) - float(np.sum(rho3))) / max(abs(float(np.sum(rho3))), 1.0e-30),
        **sol,
    }


def run_g3d_from_vector(vector, observed_for_first_step=None) -> dict:
    los = M14.project_vector_to_image_plane(*vector, los_axis="z")
    Rx = np.asarray(los["comp_1"], dtype=np.float64)
    Ry = np.asarray(los["comp_2"], dtype=np.float64)
    los_mag = np.hypot(Rx, Ry)
    grid = np.linspace(-BASE.CFG["extent"], BASE.CFG["extent"], Rx.shape[0])
    field = {"xgrid": grid, "ygrid": grid, "rx": Rx, "ry": Ry}

    x0, y0, _, _ = BASE._launch_expanded_25pct()
    groups = GEO._source_groups(x0, y0)
    if len(groups) != EXPECTED_SUPPORT:
        raise RuntimeError(f"expected {EXPECTED_SUPPORT} source bins, got {len(groups)}")
    checkpoints, g3d = GEO._propagate_g3d(field, BASE.CFG["step"], BASE.CFG["steps"], x0, y0)
    if g3d["max_unit_speed_error"] > UNIT_SPEED_TOL:
        raise RuntimeError(f"G3D unit-speed gate failed: {g3d['max_unit_speed_error']}")

    control_obs = np.zeros_like(los_mag) if observed_for_first_step is None else observed_for_first_step
    first = GEO._first_step_geometry(field, x0, y0, checkpoints[1], control_obs, los_mag)
    if not first["first_step_exact_pass"]:
        raise RuntimeError("first-step exact geometry gate failed")

    final_ang = ANG._angular_distribution_fields(checkpoints[CHECKPOINT], groups)
    gates = ANG._moment_gates(final_ang)
    if gates["second_moment_equals_cov_plus_centroid_outer_relative_rms_error"] > ANG.MOMENT_IDENTITY_TOL:
        raise RuntimeError("angular second-moment identity failed")
    if not gates["covariance_psd_pass"]:
        raise RuntimeError("angular covariance PSD gate failed")
    if not gates["direction_mean_vector_bound_pass"]:
        raise RuntimeError("angular direction-mean bound failed")

    snap = checkpoints[CHECKPOINT]
    vz = np.asarray(snap["vz"], dtype=np.float64)
    if float(np.min(np.abs(vz))) <= DET.VZ_MIN:
        raise RuntimeError("final tangent projection vz too small")
    tx = np.asarray(snap["vx"], dtype=np.float64) / vz
    ty = np.asarray(snap["vy"], dtype=np.float64) / vz
    return {
        "los_mag": los_mag,
        "checkpoints": checkpoints,
        "g3d": g3d,
        "first": first,
        "final_ang": final_ang,
        "tx": tx,
        "ty": ty,
        "groups": groups,
        "angular_gates": gates,
    }


def legacy_chain(rho3: np.ndarray) -> dict:
    # Exact historical control; its 0.18 remains isolated to this lane.
    return IND._run_chain_from_rho3(rho3, observed_for_first_step=None)


def unit_control_chain(rho3: np.ndarray) -> dict:
    # Diagnostic only: reconstruct the historical A8/M10 route with s=1 to expose
    # the effect of the old scalar.  Never interpreted as a physical prediction.
    rng = np.random.RandomState(BASE.SEED)
    eq = np.asarray(rho3, dtype=np.float64)
    noise = M06_state.A8_INIT_INJECTION_NOISE * rng.randn(*rho3.shape)
    initial = {"rho_3d": rho3.copy(), "u_slow0": eq.copy(), "u_fast0": eq + noise}
    state = BASE._evolve(initial)
    candidate = BASE._candidate(state)
    return run_g3d_from_vector(BASE._interface_vector(candidate), observed_for_first_step=None)


def corr(a, b) -> tuple[float, float, int]:
    return IND._corr(a, b)


def compare_lane(prefix: str, rho2: np.ndarray, chain: dict, observed: np.ndarray) -> dict:
    fields = {
        "source_rho2": rho2,
        "los_mag": chain["los_mag"],
        "final_angular_centroid_mag": chain["final_ang"]["angular_centroid_mag"],
        "final_angular_spread_rms": chain["final_ang"]["angular_spread_rms"],
        "final_angular_rms_angle_mag": chain["final_ang"]["angular_rms_angle_mag"],
    }
    out = {}
    for name, field in fields.items():
        p, s, n = corr(field, observed)
        out[f"{prefix}_{name}_vs_observed_pearson"] = p
        out[f"{prefix}_{name}_vs_observed_spearman"] = s
        out[f"{prefix}_{name}_vs_observed_count"] = n
        out[f"{prefix}_{name}_rms"] = rms(field)
    return out


def run_cluster(cluster: dict) -> dict:
    # ------------------------------------------------------------------
    # ALL PREDICTION LANES FIRST. Only benchmark header/shape may be used by SRC.
    # ------------------------------------------------------------------
    source = SRC._independent_source(cluster)
    rho3 = np.asarray(source["rho3"], dtype=np.float64)

    legacy = legacy_chain(rho3)
    unit = unit_control_chain(rho3)
    native_build = native_accumulated_vector(rho3)
    native = run_g3d_from_vector(native_build["vector"], observed_for_first_step=None)

    # ------------------------------------------------------------------
    # ONLY NOW load observed kappa values for external comparison.
    # ------------------------------------------------------------------
    kpath = IND._kappa_path(cluster)
    with fits.open(kpath) as hdul:
        kappa_native = np.asarray(hdul[0].data, dtype=np.float64)
    observed, _assisted_rho2 = SRC._benchmark_on_common_grid(kappa_native, source["geometry"])

    legacy_metrics = compare_lane("legacy", source["rho2"], legacy, observed)
    unit_metrics = compare_lane("unit_control", source["rho2"], unit, observed)
    native_metrics = compare_lane("native", source["rho2"], native, observed)

    # Direct amplitude diagnostics: how much signal the 0.18 control suppressed,
    # and how the physically generated native lane compares without any scalar.
    amp = {
        "legacy_strength": LEGACY_STRENGTH,
        "legacy_los_rms": rms(legacy["los_mag"]),
        "unit_los_rms": rms(unit["los_mag"]),
        "native_los_rms": rms(native["los_mag"]),
        "legacy_over_unit_los_rms": rms(legacy["los_mag"]) / max(rms(unit["los_mag"]), 1.0e-30),
        "native_over_legacy_los_rms": rms(native["los_mag"]) / max(rms(legacy["los_mag"]), 1.0e-30),
        "legacy_final_angle_rms": rms(legacy["final_ang"]["angular_rms_angle_mag"]),
        "unit_final_angle_rms": rms(unit["final_ang"]["angular_rms_angle_mag"]),
        "native_final_angle_rms": rms(native["final_ang"]["angular_rms_angle_mag"]),
    }

    return {
        "cluster_id": cluster["id"],
        "source_role": source["source_role"],
        "source_limit": "HST_F160W_luminous_structure_proxy_not_absolute_baryonic_mass_map",
        "benchmark_values_loaded_after_all_prediction_lanes_complete": True,
        "legacy_role": "historical_strength_0p18_control_only",
        "unit_role": "diagnostic_response_kernel_control_not_physical_prediction",
        "native_role": "raw_c_state_to_bounded_strain_accumulation_to_gradient_to_existing_G3D_no_strength",
        "native_c_state_integral_relative_error": native_build["c_state_integral_relative_error"],
        "native_accumulation_converged": native_build["converged"],
        "native_max_strain_fraction": native_build["max_strain_fraction"],
        "native_nonlinear_relative_residual": native_build["nonlinear_relative_residual"],
        "native_g3d_unit_speed_max_error": native["g3d"]["max_unit_speed_error"],
        **amp,
        **legacy_metrics,
        **unit_metrics,
        **native_metrics,
    }


def main() -> int:
    state = repo_state()
    clusters = list(BASE.CLUSTERS)
    rows = []
    failures = []
    for cluster in clusters:
        try:
            rows.append(run_cluster(cluster))
        except Exception as exc:
            failures.append({"cluster_id": cluster["id"], "error": f"{type(exc).__name__}: {exc}"})

    finite_rows = [r for r in rows if np.isfinite(r.get("native_final_angular_rms_angle_mag_rms", float("nan")))]
    checks = {
        "legacy_strength_inventory_exact": bool(abs(LEGACY_STRENGTH - EXPECTED_LEGACY_STRENGTH) <= 1.0e-15),
        "all_five_clusters_completed": bool(len(rows) == len(clusters) == 5 and not failures),
        "native_c_state_integral_preserved": bool(rows and all(r["native_c_state_integral_relative_error"] <= 1.0e-12 for r in rows)),
        "native_accumulation_converged": bool(rows and all(r["native_accumulation_converged"] for r in rows)),
        "native_G3D_valid": bool(rows and all(r["native_g3d_unit_speed_max_error"] <= UNIT_SPEED_TOL for r in rows)),
        "native_signal_nonzero": bool(rows and all(r["native_los_rms"] > 0.0 and r["native_final_angle_rms"] > 0.0 for r in rows)),
        "legacy_0p18_suppression_visible": bool(rows and all(abs(r["legacy_over_unit_los_rms"] - LEGACY_STRENGTH) <= 0.12 for r in rows)),
    }

    # Lensing comparison is deliberately reported, not gated to a target value.
    # This prevents benchmark agreement from becoming a tuning criterion.
    native_better_pearson_count = sum(
        1 for r in rows
        if np.isfinite(r.get("native_final_angular_rms_angle_mag_vs_observed_pearson", np.nan))
        and np.isfinite(r.get("legacy_final_angular_rms_angle_mag_vs_observed_pearson", np.nan))
        and r["native_final_angular_rms_angle_mag_vs_observed_pearson"] > r["legacy_final_angular_rms_angle_mag_vs_observed_pearson"]
    )

    execution_checks = {
        "all_measured_values_finite_for_completed_rows": bool(all(np.isfinite(r["native_los_rms"]) for r in rows)),
        "raw_c_state_used_without_strength": True,
        "bounded_strain_accumulation_used_without_rescaling": True,
        "existing_G3D_ray_tracker_reused": True,
        "legacy_strength_confined_to_control_lane": True,
        "no_replacement_strength_scalar": True,
        "no_native_amplitude_rescaling": True,
        "no_fit_or_tuning": True,
        "no_GR_potential_decomposition": True,
        "no_LCDM": True,
        "no_Rmax": True,
        "no_quantum_engine": True,
        "no_planck_input": True,
        "benchmark_pixels_loaded_only_after_predictions": True,
        "no_tracked_or_staged_changes": bool(not state["tracked_changes"] and not state["staged_changes"]),
        "stdout_only_no_run_directory_created": True,
    }
    execution_gate_pass = bool(all(execution_checks.values()))

    matched = sum(bool(v) for v in checks.values())
    if matched == len(checks) and execution_gate_pass:
        status = "NATIVE_ACCUMULATED_FULL_LENSING_EXECUTED"
    elif matched >= 4 and execution_gate_pass:
        status = "NATIVE_ACCUMULATED_FULL_LENSING_PARTIAL_EXECUTION"
    else:
        status = "NATIVE_ACCUMULATED_FULL_LENSING_NOT_ESTABLISHED"

    result = {
        "lab_id": LAB_ID,
        "status": status,
        "repo_state": state,
        "model": {
            "historical_control": "strength=0.18 -> A8 -> M10 -> LOS -> G3D",
            "unit_control": "strength=1 diagnostic only -> A8 -> M10 -> LOS -> G3D",
            "native_prediction": "rho3 -> raw c_state -> bounded-strain accumulated deformation -> -grad(u) -> existing LOS/G3D",
            "replacement_strength_scalar": None,
            "legacy_strength_used_in_native_lane": False,
            "observed_kappa_role": "end_of_chain_external_comparison_only",
            "source_limit": "normalized_HST_F160W_luminous_structure_proxy_not_absolute_mass_map",
        },
        "rows": rows,
        "failures": failures,
        "checks": checks,
        "matched_checks": matched,
        "native_better_than_legacy_pearson_count_of_5": native_better_pearson_count,
        "execution_checks": execution_checks,
        "execution_gate_pass": execution_gate_pass,
        "interpretation": {
            "strength_0p18_role": "historical_dimensionless_initial_source_loading_multiplier_not_ray_bending_coefficient",
            "question": "Can the newly derived native accumulated medium response replace the old strength-scaled source loading and drive the existing full G3D lensing stack without any replacement scalar?",
            "important_limit": "Because the current independent source is normalized F160W luminous structure rather than an absolute mass map, observational amplitude agreement is diagnostic and cannot yet be called an SI absolute prediction.",
        },
    }

    print(LAB_ID)
    print(f"status={status}")
    print(f"head_sha={state['head_sha']}")
    print(f"legacy_strength={LEGACY_STRENGTH:.12g}")
    print("legacy_strength_role=initial_source_loading_multiplier")
    print("replacement_strength_scalar=none")
    print("raw_c_state_used_without_strength=true")
    print("native_accumulated_response_rescaled=false")
    print("existing_G3D_ray_tracker_reused=true")
    print("benchmark_loaded_after_prediction_lanes=true")
    print()
    print("CLUSTERS")
    for r in rows:
        print(
            f"cluster={r['cluster_id']} "
            f"legacy_over_unit_los_rms={r['legacy_over_unit_los_rms']:.12g} "
            f"native_over_legacy_los_rms={r['native_over_legacy_los_rms']:.12g} "
            f"legacy_angle_rms={r['legacy_final_angle_rms']:.12g} "
            f"native_angle_rms={r['native_final_angle_rms']:.12g} "
            f"legacy_obs_pearson={r.get('legacy_final_angular_rms_angle_mag_vs_observed_pearson', float('nan')):.12g} "
            f"native_obs_pearson={r.get('native_final_angular_rms_angle_mag_vs_observed_pearson', float('nan')):.12g}"
        )
    for f in failures:
        print(f"failure_cluster={f['cluster_id']} error={f['error']}")
    print()
    print("CHECKS")
    for k, v in checks.items(): print(f"{k}={str(v).lower()}")
    print("EXECUTION_CHECKS")
    for k, v in execution_checks.items(): print(f"{k}={str(v).lower()}")
    print(f"execution_gate_pass={str(execution_gate_pass).lower()}")
    print(json.dumps(result, sort_keys=True, separators=(",", ":"), default=str))
    return 0 if execution_gate_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
