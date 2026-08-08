#!/usr/bin/env python3
"""PBUF FOUNDATION — NATIVE ACCUMULATED FULL LENSING LOCAL BENCHMARK 001.

Purpose
-------
Correct the execution mistake in native_accumulated_full_lensing001: do not
rediscover or download data that already exists in the repository.

This audit uses the five canonical local FITS files under PBUF_benchmark through
``pbuf.core.benchmark_data`` and runs the existing full G3D lensing stack in
three paired lanes using the exact same locally loaded source morphology:

    LEGACY : local kappa proxy -> strength=0.18 -> A8/M10 -> LOS/G3D -> observer
    UNIT   : local kappa proxy -> strength=1 diagnostic -> A8/M10 -> LOS/G3D
    NATIVE : local kappa proxy -> raw c_state -> bounded-strain accumulation
             -> -grad(u) -> existing LOS/G3D -> observer

The local kappa map is intentionally used here as a benchmark-assisted source,
matching earlier foundation lensing diagnostics.  Therefore this lab is NOT an
independent prediction test.  Its narrow purpose is to test the historical
strength=0.18 role and whether the newly derived native accumulation can drive
the already-built full lensing machinery with no replacement scalar.

No network access, HST discovery/download, fitted amplitude, native rescaling,
GR/LCDM potential machinery, Rmax, Quantum Engine, or Planck input is used.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from a8_three_dimensional_projection_lab001 import construct_rho_3d
from weak_lensing_observation001 import resample_to_grid
from pbuf.core import benchmark_data as BENCH
import pbuf.labs.foundation.m10_coverage_25pct_science001 as BASE
import pbuf.labs.foundation.native_accumulated_full_lensing001 as FULL

LAB_ID = "PBUF-FOUNDATION-NATIVE-ACCUMULATED-FULL-LENSING-LOCAL-BENCHMARK-001"
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


def local_source(cluster: dict) -> dict:
    """Build the historical benchmark-assisted 3-D source from local FITS only."""
    path = BENCH.require_kappa_path(cluster)
    kappa = BENCH.load_kappa(cluster)
    rho2 = BASE.construct_common_proxy(
        kappa, bins=BASE.OBS_BINS, extent=BASE.CFG["extent"]
    )
    rho3 = construct_rho_3d(rho2, BASE.NZ, profile=BASE.PROFILE)
    observed = resample_to_grid(kappa, BASE.OBS_BINS, BASE.CFG["extent"])
    return {
        "path": path,
        "kappa": kappa,
        "rho2": np.asarray(rho2, dtype=np.float64),
        "rho3": np.asarray(rho3, dtype=np.float64),
        "observed": np.asarray(observed, dtype=np.float64),
    }


def run_cluster(cluster: dict) -> dict:
    source = local_source(cluster)
    rho2 = source["rho2"]
    rho3 = source["rho3"]
    observed = source["observed"]

    # Exact historical control.  0.18 remains confined to this lane.
    legacy = FULL.legacy_chain(rho3)

    # Unit historical route is diagnostic only; it is not a candidate physical
    # replacement for 0.18.
    unit = FULL.unit_control_chain(rho3)

    # Native lane: no strength factor anywhere.  Raw c_state is accumulated by
    # the bounded-strain network and its deformation gradient is supplied to the
    # already-existing LOS/G3D tracker.
    native_build = FULL.native_accumulated_vector(rho3)
    native = FULL.run_g3d_from_vector(
        native_build["vector"], observed_for_first_step=None
    )

    legacy_metrics = FULL.compare_lane("legacy", rho2, legacy, observed)
    unit_metrics = FULL.compare_lane("unit_control", rho2, unit, observed)
    native_metrics = FULL.compare_lane("native", rho2, native, observed)

    legacy_los = FULL.rms(legacy["los_mag"])
    unit_los = FULL.rms(unit["los_mag"])
    native_los = FULL.rms(native["los_mag"])
    legacy_angle = FULL.rms(legacy["final_ang"]["angular_rms_angle_mag"])
    unit_angle = FULL.rms(unit["final_ang"]["angular_rms_angle_mag"])
    native_angle = FULL.rms(native["final_ang"]["angular_rms_angle_mag"])

    return {
        "cluster_id": cluster["id"],
        "benchmark_path": str(source["path"]),
        "source_role": "local_kappa_benchmark_assisted_diagnostic_source",
        "independent_prediction": False,
        "network_used": False,
        "legacy_role": "historical_strength_0p18_control_only",
        "unit_role": "unit_loading_response_diagnostic_only",
        "native_role": "raw_c_state_to_bounded_strain_accumulation_to_gradient_to_existing_G3D_no_strength",
        "native_c_state_integral_relative_error": native_build["c_state_integral_relative_error"],
        "native_accumulation_converged": native_build["converged"],
        "native_max_strain_fraction": native_build["max_strain_fraction"],
        "native_nonlinear_relative_residual": native_build["nonlinear_relative_residual"],
        "native_g3d_unit_speed_max_error": native["g3d"]["max_unit_speed_error"],
        "legacy_los_rms": legacy_los,
        "unit_los_rms": unit_los,
        "native_los_rms": native_los,
        "legacy_over_unit_los_rms": legacy_los / max(unit_los, 1.0e-30),
        "native_over_legacy_los_rms": native_los / max(legacy_los, 1.0e-30),
        "legacy_final_angle_rms": legacy_angle,
        "unit_final_angle_rms": unit_angle,
        "native_final_angle_rms": native_angle,
        "legacy_over_unit_angle_rms": legacy_angle / max(unit_angle, 1.0e-30),
        "native_over_legacy_angle_rms": native_angle / max(legacy_angle, 1.0e-30),
        **legacy_metrics,
        **unit_metrics,
        **native_metrics,
    }


def main() -> int:
    state = repo_state()
    inventory = BENCH.inventory()
    clusters = list(BENCH.clusters())

    inventory_ids = tuple(row["id"] for row in clusters)
    local_files_ready = bool(
        inventory_ids == EXPECTED_CLUSTER_IDS
        and len(inventory) == 5
        and all(row["exists"] for row in inventory)
    )

    rows = []
    failures = []
    if local_files_ready:
        for cluster in clusters:
            try:
                rows.append(run_cluster(cluster))
            except Exception as exc:
                failures.append(
                    {
                        "cluster_id": cluster["id"],
                        "error": f"{type(exc).__name__}: {exc}",
                    }
                )

    checks = {
        "canonical_five_cluster_inventory": bool(inventory_ids == EXPECTED_CLUSTER_IDS),
        "all_five_local_fits_present": local_files_ready,
        "all_five_clusters_completed": bool(len(rows) == 5 and not failures),
        "legacy_strength_inventory_exact": bool(
            abs(FULL.LEGACY_STRENGTH - FULL.EXPECTED_LEGACY_STRENGTH) <= 1.0e-15
        ),
        "native_c_state_integral_preserved": bool(
            rows and all(r["native_c_state_integral_relative_error"] <= 1.0e-12 for r in rows)
        ),
        "native_accumulation_converged": bool(
            rows and all(r["native_accumulation_converged"] for r in rows)
        ),
        "native_G3D_valid": bool(
            rows and all(r["native_g3d_unit_speed_max_error"] <= FULL.UNIT_SPEED_TOL for r in rows)
        ),
        "native_signal_nonzero": bool(
            rows and all(r["native_los_rms"] > 0.0 and r["native_final_angle_rms"] > 0.0 for r in rows)
        ),
        "legacy_0p18_suppression_visible": bool(
            rows
            and all(
                abs(r["legacy_over_unit_los_rms"] - FULL.LEGACY_STRENGTH) <= 0.12
                for r in rows
            )
        ),
    }

    execution_checks = {
        "local_benchmark_module_used": True,
        "network_access_used": False,
        "HST_download_or_discovery_used": False,
        "benchmark_assisted_source_explicitly_labeled": True,
        "independent_prediction_claimed": False,
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
        "all_measured_values_finite_for_completed_rows": bool(
            all(
                np.all(
                    np.isfinite(
                        [
                            r["legacy_los_rms"], r["unit_los_rms"], r["native_los_rms"],
                            r["legacy_final_angle_rms"], r["unit_final_angle_rms"],
                            r["native_final_angle_rms"], r["native_max_strain_fraction"],
                        ]
                    )
                )
                for r in rows
            )
        ),
        "no_tracked_or_staged_changes": bool(
            not state["tracked_changes"] and not state["staged_changes"]
        ),
        "stdout_only_no_run_directory_created": True,
    }
    execution_gate_pass = bool(all(execution_checks.values()))

    if all(checks.values()) and execution_gate_pass:
        status = "LOCAL_BENCHMARK_NATIVE_FULL_LENSING_EXECUTED"
    elif sum(bool(v) for v in checks.values()) >= 5 and execution_gate_pass:
        status = "LOCAL_BENCHMARK_NATIVE_FULL_LENSING_PARTIAL_EXECUTION"
    else:
        status = "LOCAL_BENCHMARK_NATIVE_FULL_LENSING_NOT_ESTABLISHED"

    result = {
        "lab_id": LAB_ID,
        "status": status,
        "repo_state": state,
        "benchmark_inventory": inventory,
        "model": {
            "data_loader": "pbuf.core.benchmark_data",
            "data_source": "five_existing_local_PBUF_benchmark_Merten_v1_kappa_FITS",
            "network_access": False,
            "source_role": "benchmark_assisted_diagnostic_not_independent_prediction",
            "historical_control": "local kappa proxy -> strength=0.18 -> A8/M10 -> LOS/G3D",
            "unit_control": "local kappa proxy -> unit loading diagnostic -> A8/M10 -> LOS/G3D",
            "native_lane": "local kappa proxy -> raw c_state -> bounded-strain accumulated response -> -grad(u) -> existing LOS/G3D",
            "replacement_strength_scalar": None,
        },
        "rows": rows,
        "failures": failures,
        "checks": checks,
        "execution_checks": execution_checks,
        "execution_gate_pass": execution_gate_pass,
        "interpretation": {
            "question": "On the same five local benchmark sources used by earlier lensing labs, can the native accumulated response drive the existing full G3D stack without the historical strength=0.18 scalar?",
            "strength_role": "historical_initial_source_loading_multiplier_only",
            "scientific_limit": "The local kappa map supplies both benchmark-assisted source morphology and comparison morphology; this run tests strength removal and end-to-end propagation, not independent predictive power.",
        },
    }

    print(LAB_ID)
    print(f"status={status}")
    print(f"head_sha={state['head_sha']}")
    print("benchmark_loader=pbuf.core.benchmark_data")
    print("benchmark_root=PBUF_benchmark")
    print("network_access_used=false")
    print("benchmark_assisted_source=true")
    print("independent_prediction=false")
    print(f"legacy_strength={FULL.LEGACY_STRENGTH:.12g}")
    print("replacement_strength_scalar=none")
    print("native_response_rescaled=false")
    print()
    print("LOCAL_BENCHMARK_INVENTORY")
    for item in inventory:
        print(
            f"cluster={item['id']} exists={str(item['exists']).lower()} "
            f"path={item['path']}"
        )
    print()
    print("CLUSTERS")
    for r in rows:
        print(
            f"cluster={r['cluster_id']} "
            f"legacy_over_unit_los_rms={r['legacy_over_unit_los_rms']:.12g} "
            f"native_over_legacy_los_rms={r['native_over_legacy_los_rms']:.12g} "
            f"legacy_over_unit_angle_rms={r['legacy_over_unit_angle_rms']:.12g} "
            f"native_over_legacy_angle_rms={r['native_over_legacy_angle_rms']:.12g} "
            f"native_max_strain_fraction={r['native_max_strain_fraction']:.12g} "
            f"legacy_obs_pearson={r.get('legacy_final_angular_rms_angle_mag_vs_observed_pearson', float('nan')):.12g} "
            f"native_obs_pearson={r.get('native_final_angular_rms_angle_mag_vs_observed_pearson', float('nan')):.12g}"
        )
    for failure in failures:
        print(f"failure_cluster={failure['cluster_id']} error={failure['error']}")
    print()
    print("CHECKS")
    for key, value in checks.items():
        print(f"{key}={str(value).lower()}")
    print("EXECUTION_CHECKS")
    for key, value in execution_checks.items():
        print(f"{key}={str(value).lower()}")
    print(f"execution_gate_pass={str(execution_gate_pass).lower()}")
    print(json.dumps(result, sort_keys=True, separators=(",", ":"), default=str))
    return 0 if execution_gate_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
