#!/usr/bin/env python3
"""PBUF FOUNDATION — NATIVE PROPAGATION INTERFACE AMPLITUDE AUDIT 001.

Purpose
-------
Locate where the ~10^2 native/unit amplitude mismatch first appears without
changing any physics.  On the same five canonical local benchmark sources, this
lab reconstructs the old unit-loading A8/M10 route and the new native
zero-flux-c_state/bounded-strain route, then records amplitudes at each interface:

    source rho3
      -> A8 c_state
      -> propagation-input 3-vector
      -> LOS projected 2-vector
      -> G3D final angular response

The native route additionally reports accumulated deformation u itself between
c_state and its gradient.  No normalization, fitting, rescaling, replacement
strength, or observed-amplitude target is used.  Observed kappa values are not
needed at all: this is an internal field-interface audit.
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

from a8_three_dimensional_projection_lab001 import construct_rho_3d
from pbuf.core import benchmark_data as BENCH
from pbuf.core import los_projection as M14
from pbuf.models import a8_state as M06_state
import pbuf.labs.foundation.m10_coverage_25pct_science001 as BASE
import pbuf.labs.foundation.native_accumulated_full_lensing001 as FULL
import pbuf.labs.foundation.native_accumulated_full_lensing_local_benchmark001 as LOCAL

LAB_ID = "PBUF-FOUNDATION-NATIVE-PROPAGATION-INTERFACE-AMPLITUDE-AUDIT-001"
EXPECTED_CLUSTER_IDS = (
    "Abell2744", "MACS0416", "MACS1149", "AbellS1063", "Abell370"
)


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
    return float(np.sqrt(np.mean(x * x))) if x.size else float("nan")


def vector_rms(vector) -> float:
    comps = [np.asarray(v, dtype=np.float64) for v in vector]
    mag2 = np.zeros_like(comps[0])
    for comp in comps:
        mag2 += comp * comp
    return float(np.sqrt(np.mean(mag2)))


def max_vector_mag(vector) -> float:
    comps = [np.asarray(v, dtype=np.float64) for v in vector]
    mag2 = np.zeros_like(comps[0])
    for comp in comps:
        mag2 += comp * comp
    return float(np.sqrt(np.max(mag2)))


def ratio(a: float, b: float) -> float:
    return float(a / max(abs(b), 1.0e-30))


def local_rho3(cluster: dict) -> np.ndarray:
    kappa = BENCH.load_kappa(cluster)
    rho2 = BASE.construct_common_proxy(
        kappa, bins=BASE.OBS_BINS, extent=BASE.CFG["extent"]
    )
    return np.asarray(
        construct_rho_3d(rho2, BASE.NZ, profile=BASE.PROFILE),
        dtype=np.float64,
    )


def unit_intermediates(rho3: np.ndarray) -> dict:
    """Reconstruct exactly the historical route with source loading set to one."""
    rng = np.random.RandomState(BASE.SEED)
    eq = np.asarray(rho3, dtype=np.float64)
    noise = M06_state.A8_INIT_INJECTION_NOISE * rng.randn(*rho3.shape)
    initial = {
        "rho_3d": rho3.copy(),
        "u_slow0": eq.copy(),
        "u_fast0": eq + noise,
    }
    state = BASE._evolve(initial)
    candidate = BASE._candidate(state)
    vector = tuple(np.asarray(v, dtype=np.float64) for v in BASE._interface_vector(candidate))
    los = M14.project_vector_to_image_plane(*vector, los_axis="z")
    los_vector = (
        np.asarray(los["comp_1"], dtype=np.float64),
        np.asarray(los["comp_2"], dtype=np.float64),
    )
    chain = FULL.run_g3d_from_vector(vector, observed_for_first_step=None)
    return {
        "c_state": np.asarray(state["c_state"], dtype=np.float64),
        "vector": vector,
        "los_vector": los_vector,
        "chain": chain,
    }


def native_intermediates(rho3: np.ndarray) -> dict:
    """Use the frozen zero-flux native bridge from the successful five-cluster run."""
    build = LOCAL.native_accumulated_vector_zero_flux(rho3)
    vector = tuple(np.asarray(v, dtype=np.float64) for v in build["vector"])
    los = M14.project_vector_to_image_plane(*vector, los_axis="z")
    los_vector = (
        np.asarray(los["comp_1"], dtype=np.float64),
        np.asarray(los["comp_2"], dtype=np.float64),
    )
    chain = FULL.run_g3d_from_vector(vector, observed_for_first_step=None)
    return {
        "c_state": np.asarray(build["c_state"], dtype=np.float64),
        "accumulated": np.asarray(build["accumulated"], dtype=np.float64),
        "vector": vector,
        "los_vector": los_vector,
        "chain": chain,
        "c_state_integral_relative_error": float(build["c_state_integral_relative_error"]),
        "max_strain_fraction": float(build["max_strain_fraction"]),
        "accumulation_converged": bool(build["converged"]),
    }


def final_angle_rms(chain: dict) -> float:
    return FULL.rms(chain["final_ang"]["angular_rms_angle_mag"])


def run_cluster(cluster: dict) -> dict:
    rho3 = local_rho3(cluster)
    unit = unit_intermediates(rho3)
    native = native_intermediates(rho3)

    source_rms = rms(rho3)
    unit_c_rms = rms(unit["c_state"])
    native_c_rms = rms(native["c_state"])
    native_u_rms = rms(native["accumulated"])
    unit_vec_rms = vector_rms(unit["vector"])
    native_vec_rms = vector_rms(native["vector"])
    unit_los_rms = vector_rms(unit["los_vector"])
    native_los_rms = vector_rms(native["los_vector"])
    unit_angle = final_angle_rms(unit["chain"])
    native_angle = final_angle_rms(native["chain"])

    c_ratio = ratio(native_c_rms, unit_c_rms)
    vec_ratio = ratio(native_vec_rms, unit_vec_rms)
    los_ratio = ratio(native_los_rms, unit_los_rms)
    angle_ratio = ratio(native_angle, unit_angle)

    stage_ratios = {
        "c_state": c_ratio,
        "propagation_input_vector": vec_ratio,
        "los_vector": los_ratio,
        "final_angle": angle_ratio,
    }
    first_gt10 = next((name for name, value in stage_ratios.items() if value >= 10.0), "none")

    return {
        "cluster_id": cluster["id"],
        "source_rho3_rms": source_rms,
        "unit_c_state_rms": unit_c_rms,
        "native_c_state_rms": native_c_rms,
        "native_over_unit_c_state_rms": c_ratio,
        "native_accumulated_u_rms": native_u_rms,
        "native_accumulated_u_max_abs": float(np.max(np.abs(native["accumulated"]))),
        "unit_m10_vector_rms": unit_vec_rms,
        "native_gradient_vector_rms": native_vec_rms,
        "native_over_unit_interface_vector_rms": vec_ratio,
        "unit_m10_vector_max_mag": max_vector_mag(unit["vector"]),
        "native_gradient_vector_max_mag": max_vector_mag(native["vector"]),
        "unit_los_vector_rms": unit_los_rms,
        "native_los_vector_rms": native_los_rms,
        "native_over_unit_los_vector_rms": los_ratio,
        "unit_final_angle_rms": unit_angle,
        "native_final_angle_rms": native_angle,
        "native_over_unit_final_angle_rms": angle_ratio,
        "native_gradient_rms_over_native_c_state_rms": ratio(native_vec_rms, native_c_rms),
        "unit_m10_vector_rms_over_unit_c_state_rms": ratio(unit_vec_rms, unit_c_rms),
        "native_gradient_rms_over_accumulated_u_rms": ratio(native_vec_rms, native_u_rms),
        "first_native_over_unit_stage_ge_10x": first_gt10,
        "native_c_state_integral_relative_error": native["c_state_integral_relative_error"],
        "native_max_strain_fraction": native["max_strain_fraction"],
        "native_accumulation_converged": native["accumulation_converged"],
        "native_g3d_unit_speed_max_error": float(native["chain"]["g3d"]["max_unit_speed_error"]),
        "unit_g3d_unit_speed_max_error": float(unit["chain"]["g3d"]["max_unit_speed_error"]),
    }


def main() -> int:
    state = repo_state()
    clusters = list(BENCH.clusters())
    inventory = BENCH.inventory()
    ids = tuple(c["id"] for c in clusters)

    rows = []
    failures = []
    if ids == EXPECTED_CLUSTER_IDS and len(inventory) == 5 and all(x["exists"] for x in inventory):
        for cluster in clusters:
            try:
                rows.append(run_cluster(cluster))
            except Exception as exc:
                failures.append({
                    "cluster_id": cluster["id"],
                    "error": f"{type(exc).__name__}: {exc}",
                })

    checks = {
        "canonical_five_cluster_inventory": ids == EXPECTED_CLUSTER_IDS,
        "all_five_local_fits_present": bool(len(inventory) == 5 and all(x["exists"] for x in inventory)),
        "all_five_clusters_completed": bool(len(rows) == 5 and not failures),
        "all_measured_stage_values_finite": bool(rows and all(
            all(math.isfinite(float(r[k])) for k in (
                "source_rho3_rms", "unit_c_state_rms", "native_c_state_rms",
                "native_accumulated_u_rms", "unit_m10_vector_rms",
                "native_gradient_vector_rms", "unit_los_vector_rms",
                "native_los_vector_rms", "unit_final_angle_rms",
                "native_final_angle_rms",
            )) for r in rows
        )),
        "native_c_state_integral_preserved": bool(rows and all(r["native_c_state_integral_relative_error"] <= 1.0e-12 for r in rows)),
        "native_accumulation_converged": bool(rows and all(r["native_accumulation_converged"] for r in rows)),
        "unit_and_native_G3D_valid": bool(rows and all(
            r["native_g3d_unit_speed_max_error"] <= FULL.UNIT_SPEED_TOL
            and r["unit_g3d_unit_speed_max_error"] <= FULL.UNIT_SPEED_TOL
            for r in rows
        )),
        "no_tracked_or_staged_changes": bool(not state["tracked_changes"] and not state["staged_changes"]),
    }
    execution_gate_pass = bool(all(checks.values()))

    if execution_gate_pass:
        status = "NATIVE_PROPAGATION_INTERFACE_AMPLITUDE_AUDIT_EXECUTED"
    elif rows:
        status = "NATIVE_PROPAGATION_INTERFACE_AMPLITUDE_AUDIT_PARTIAL_EXECUTION"
    else:
        status = "NATIVE_PROPAGATION_INTERFACE_AMPLITUDE_AUDIT_NOT_ESTABLISHED"

    first_stage_counts = {}
    for row in rows:
        key = row["first_native_over_unit_stage_ge_10x"]
        first_stage_counts[key] = first_stage_counts.get(key, 0) + 1

    result = {
        "lab_id": LAB_ID,
        "status": status,
        "repo_state": state,
        "model": {
            "unit_route": "rho3 -> unit A8 c_state -> M10 interface vector -> LOS -> existing G3D",
            "native_route": "rho3 -> zero-flux raw c_state -> bounded-strain accumulated u -> -grad(u) -> LOS -> existing G3D",
            "observed_lensing_values_used": False,
            "replacement_strength_scalar": None,
            "normalization_or_rescaling": False,
            "fit_or_tuning": False,
        },
        "rows": rows,
        "failures": failures,
        "first_stage_ge_10x_counts": first_stage_counts,
        "checks": checks,
        "execution_gate_pass": execution_gate_pass,
        "interpretation_rule": "Report where the native/unit ratio first exceeds 10x; this threshold is diagnostic only and does not tune or gate the model.",
    }

    print(LAB_ID)
    print(f"status={status}")
    print(f"head_sha={state['head_sha']}")
    print("benchmark_loader=pbuf.core.benchmark_data")
    print("network_access_used=false")
    print("observed_lensing_values_used=false")
    print("replacement_strength_scalar=none")
    print("native_response_rescaled=false")
    print("fit_or_tuning=false")
    print()
    print("CLUSTERS")
    for r in rows:
        print(
            f"cluster={r['cluster_id']} "
            f"c_ratio={r['native_over_unit_c_state_rms']:.12g} "
            f"interface_ratio={r['native_over_unit_interface_vector_rms']:.12g} "
            f"los_ratio={r['native_over_unit_los_vector_rms']:.12g} "
            f"angle_ratio={r['native_over_unit_final_angle_rms']:.12g} "
            f"native_u_rms={r['native_accumulated_u_rms']:.12g} "
            f"native_grad_over_c={r['native_gradient_rms_over_native_c_state_rms']:.12g} "
            f"unit_m10_over_c={r['unit_m10_vector_rms_over_unit_c_state_rms']:.12g} "
            f"first_ge_10x={r['first_native_over_unit_stage_ge_10x']}"
        )
    for failure in failures:
        print(f"failure_cluster={failure['cluster_id']} error={failure['error']}")
    print()
    print("FIRST_STAGE_GE_10X_COUNTS")
    for key in sorted(first_stage_counts):
        print(f"stage={key} count={first_stage_counts[key]}")
    print("CHECKS")
    for key, value in checks.items():
        print(f"{key}={str(value).lower()}")
    print(f"execution_gate_pass={str(execution_gate_pass).lower()}")
    print(json.dumps(result, sort_keys=True, separators=(",", ":"), default=str))
    return 0 if execution_gate_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
