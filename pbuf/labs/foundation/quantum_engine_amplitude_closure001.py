#!/usr/bin/env python3
"""PBUF FOUNDATION — QUANTUM ENGINE AMPLITUDE CLOSURE 001.

Read the frozen Quantum Engine V11 outputs directly from pbuf/data/quantum/
and determine whether they close the absolute local PBUF response amplitude
without fitting. No archive transport layer is used.
"""
from __future__ import annotations

import hashlib
import json
import math
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "pbuf" / "data" / "quantum"
OUT = ROOT / "runs" / "quantum_engine_amplitude_closure001_raw_fix"
LAB_ID = "PBUF-FOUNDATION-QUANTUM-ENGINE-AMPLITUDE-CLOSURE-001-RAW"
REL_TOL = 1e-12
ABS_TOL = 1e-15

# Original raw thermal-table provenance from the supplied QE snapshot.
THERMAL_EXPECTED_SIZE = 444514
THERMAL_EXPECTED_SHA256 = "1de56b28474ba453eb9d1b3b12cf9e2b588cc0acb0d8a4193199fabdcdd6e2bd"

REQUIRED_FILES = {
    "config.json": DATA / "config.json",
    "defaults.yaml": DATA / "config" / "defaults.yaml",
    "micro_cache.json": DATA / "micro_cache.json",
    "optimisation_result.json": DATA / "optimisation_result.json",
    "thermal_table_cache.json": DATA / "thermal_table_cache.json",
}


def git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def write_json(name: str, obj) -> None:
    (OUT / name).write_text(json.dumps(obj, indent=2, sort_keys=False) + "\n")


def sha256(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def close(a, b) -> bool:
    return math.isclose(float(a), float(b), rel_tol=REL_TOL, abs_tol=ABS_TOL)


def direct_snapshot_integrity() -> tuple[dict, dict[str, bytes]]:
    checks = {}
    payloads: dict[str, bytes] = {}
    for name, path in REQUIRED_FILES.items():
        if not path.is_file():
            raise RuntimeError(f"required QE file missing: {path.relative_to(ROOT)}")
        raw = path.read_bytes()
        payloads[name] = raw
        checks[name] = {
            "path": str(path.relative_to(ROOT)),
            "size_bytes": len(raw),
            "sha256": sha256(raw),
            "present": True,
        }

    thermal = checks["thermal_table_cache.json"]
    thermal["expected_size_bytes"] = THERMAL_EXPECTED_SIZE
    thermal["expected_sha256"] = THERMAL_EXPECTED_SHA256
    thermal["size_provenance_pass"] = thermal["size_bytes"] == THERMAL_EXPECTED_SIZE
    thermal["sha256_provenance_pass"] = thermal["sha256"] == THERMAL_EXPECTED_SHA256
    thermal["provenance_pass"] = bool(
        thermal["size_provenance_pass"] and thermal["sha256_provenance_pass"]
    )

    return {
        "storage_mode": "direct_raw_repository_files",
        "archive_or_base64_layer_used": False,
        "checks": checks,
        "all_required_files_present": True,
        "thermal_original_provenance_verified": thermal["provenance_pass"],
        "all_pass": thermal["provenance_pass"],
    }, payloads


def thermal_audit(micro: dict, thermal: dict) -> dict:
    rows = thermal["rows"]
    alpha_qm = float(micro["alpha_qm"])
    t0 = float(micro["t_min"])
    alpha_rel = []
    aT_rel = []
    flags = {}
    for row in rows:
        eps = float(row["epsilon0_T"])
        alpha = float(row["alpha_T"])
        expected = alpha_qm * eps
        alpha_rel.append(abs(alpha - expected) / max(abs(expected), 1e-300))
        at = float(row["a"]) * float(row["T_K"])
        aT_rel.append(abs(at - t0) / max(abs(t0), 1e-300))
        flag = str(row.get("validity_flag", "missing"))
        flags[flag] = flags.get(flag, 0) + 1
    max_alpha = max(alpha_rel)
    max_at = max(aT_rel)
    return {
        "row_count": len(rows),
        "alpha_T_equals_alpha_qm_times_epsilon0_T": {
            "max_relative_error": max_alpha,
            "pass": max_alpha <= REL_TOL,
        },
        "a_times_T_equals_t_min": {
            "max_relative_error": max_at,
            "pass": max_at <= REL_TOL,
            "role": "bookkeeping_identity_not_local_stiffness_law",
        },
        "validity_flag_counts": flags,
    }


def main() -> int:
    started = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)

    repo = {
        "branch": git("rev-parse", "--abbrev-ref", "HEAD"),
        "head_sha": git("rev-parse", "HEAD"),
        "tracked_changes": git("diff", "--name-only"),
        "staged_changes": git("diff", "--name-only", "--cached"),
        "stash_list": git("stash", "list"),
    }
    if repo["branch"] != "main":
        raise RuntimeError(f"runner must execute on main, got {repo['branch']}")
    if repo["tracked_changes"] or repo["staged_changes"]:
        raise RuntimeError("tracked or staged repository changes present")

    integ, payloads = direct_snapshot_integrity()
    write_json("snapshot_integrity.json", integ)
    if not integ["all_pass"]:
        raise RuntimeError("direct Quantum Engine snapshot integrity/provenance gate failed")

    try:
        config = json.loads(payloads["config.json"].decode("utf-8"))
        micro = json.loads(payloads["micro_cache.json"].decode("utf-8"))
        optimization = json.loads(payloads["optimisation_result.json"].decode("utf-8"))
        thermal = json.loads(payloads["thermal_table_cache.json"].decode("utf-8"))
    except Exception as exc:
        raise RuntimeError(f"direct QE JSON parse failed: {exc}") from exc

    meta = thermal["metadata"]
    rows = thermal["rows"]

    consistency = {
        "engine_source_match": micro["engine_source"] == meta["micro_source"] == "quantum_engine_v11",
        "method_version_match": int(config["method_version"]) == int(micro["method_version"]) == int(meta["method_version"]),
        "table_version_match": int(config["table_version"]) == int(micro["table_version"]) == int(meta["table_version"]),
        "temperature_point_count_match": int(config["temperature_points"]) == int(micro["temperature_points"]) == int(meta["num_points"]) == len(rows),
        "alpha_qm_match": close(micro["alpha_qm"], meta["alpha_qm"]),
        "eps0_base_match": close(micro["eps0_base"], meta["eps0_base"]),
        "f_cut_match": close(micro["f_cut"], meta["f_cut_T"]),
        "f_coup_match": close(micro["f_coup"], meta["f_coup_T"]),
    }
    consistency["all_pass"] = all(consistency.values())

    thermal_result = thermal_audit(micro, thermal)
    write_json("thermal_identity_audit.json", thermal_result)

    numerator = float(micro["f_coup"]) * float(micro["f_cut"])**2 * float(micro["eps0_base"])
    c_eff = numerator / float(micro["alpha_qm"])

    inventory = {
        "source_directory": "pbuf/data/quantum",
        "storage_mode": "direct_raw_repository_files",
        "engine_source": micro["engine_source"],
        "regulator": micro["regulator"],
        "field_content": micro["field_content"],
        "alpha_qm": float(micro["alpha_qm"]),
        "epsilon0_base": float(micro["eps0_base"]),
        "f_cut": float(micro["f_cut"]),
        "f_coup": float(micro["f_coup"]),
        "mixing_strength": float(micro["mixing_strength"]),
        "beta": float(micro["beta"]),
        "power_index": float(micro["power_index"]),
        "T_star_K": float(micro["T_star"]),
        "C_eff_output_inferred": c_eff,
        "C_eff_formula": "f_coup*f_cut^2*eps0_base/alpha_qm",
        "C_definition_present_in_snapshot": False,
        "C_eff_authorized_as_physical_constant": False,
        "C_eff_authorized_as_local_modulus": False,
        "alpha_qm_authorized_as_local_modulus": False,
        "epsilon0_authorized_as_local_modulus": False,
        "absolute_local_energy_density_scale_present": False,
        "microscopic_action_normalization_present": False,
        "microscopic_strain_coupling_derivative_present": False,
        "local_tensor_Hessian_present": False,
        "retarded_response_kernel_present": False,
        "coarse_graining_to_metric_strain_present": False,
    }
    write_json("quantum_amplitude_inventory.json", inventory)

    opt = optimization["result"]
    cmb = opt["dataset_results"]["cmb"]
    quarantine = {
        "artifact_present": True,
        "phase6a_passed": bool(opt["phase6a_passed"]),
        "cmb_status": cmb["status"],
        "historical_optimization_used_as_physics_input": False,
        "optimization_best_parameters_used": False,
        "optimization_alpha_used": False,
        "optimization_Rmax_used": False,
        "quarantine_pass": bool((not opt["phase6a_passed"]) and cmb["status"] == "sanity_failed"),
    }

    closure = {
        "normalized_QE_amplitude_chain_present": True,
        "normalized_QE_thermal_state_present": True,
        "physical_amplitude_bridge_complete": False,
        "replace_legacy_strength_now_authorized": False,
        "physical_cluster_lensing_run_authorized": False,
        "missing_required_components": [
            "effective_action_normalization_C_definition_and_units",
            "absolute_microscopic_energy_or_action_density_scale",
            "microscopic_mode_to_metric_strain_coupling_derivative",
            "local_tensor_Hessian_K_or_retarded_kernel_G_R",
            "coarse_graining_normalization_to_finite_metric_strain",
        ],
    }
    write_json("local_kernel_closure_status.json", closure)

    next_interface = {
        "target": "expose_QE_effective_action_normalization_and_local_metric_strain_response_interface",
        "required": [
            "C definition/provenance/dimensions/regulator dependence",
            "absolute microscopic energy or action density scale",
            "delta S_micro/delta chi_mn or dE_mode/dchi_mn",
            "K_mn_ab or causal G_R",
            "dimension-preserving coarse-graining to finite chi_mn",
        ],
        "rule": "if absent in external QE, mark missing; never infer by fitting",
    }
    write_json("next_required_quantum_engine_interface.json", next_interface)

    gates = {
        "direct_snapshot_integrity_pass": integ["all_pass"],
        "snapshot_internal_consistency_pass": consistency["all_pass"],
        "thermal_alpha_identity_pass": thermal_result["alpha_T_equals_alpha_qm_times_epsilon0_T"]["pass"],
        "thermal_scale_factor_bookkeeping_pass": thermal_result["a_times_T_equals_t_min"]["pass"],
        "historical_failed_optimization_quarantined": quarantine["quarantine_pass"],
        "no_fitted_numbers_introduced": True,
        "legacy_strength_not_used": True,
        "replacement_strength_not_selected": True,
        "benchmark_pixel_values_not_loaded": True,
        "lensing_pipeline_not_executed": True,
        "hst_mass_conversion_not_executed": True,
        "C_eff_not_promoted_to_physics": True,
        "local_modulus_not_selected": True,
    }
    if not all(gates.values()):
        write_json("validation_failure.json", {
            "gates": gates,
            "consistency": consistency,
            "optimization_quarantine": quarantine,
        })
        raise RuntimeError("Quantum Engine amplitude closure gate failed")

    outcome = "Outcome B — QE NORMALIZED AMPLITUDE CHAIN VERIFIED; ABSOLUTE LOCAL TENSOR RESPONSE NORMALIZATION STILL OPEN"
    validation = {
        "lab_id": LAB_ID,
        "outcome": outcome,
        "head_sha": repo["head_sha"],
        "QE_source_directory": "pbuf/data/quantum",
        "QE_storage_mode": "direct_raw_repository_files",
        "archive_or_base64_layer_used": False,
        "benchmark_pixel_values_loaded": False,
        "lensing_pipeline_executed": False,
        "hst_mass_conversion_executed": False,
        "fitted_numbers_introduced": False,
        "legacy_strength_used": False,
        "replacement_strength_selected": False,
        "QE_snapshot_integrity_verified": True,
        "QE_snapshot_internal_consistency_verified": True,
        "QE_normalized_amplitude_chain_verified": True,
        "QE_thermal_alpha_identity_verified": True,
        "QE_alpha_qm": float(micro["alpha_qm"]),
        "QE_eps0_base": float(micro["eps0_base"]),
        "QE_f_cut": float(micro["f_cut"]),
        "QE_f_coup": float(micro["f_coup"]),
        "QE_C_eff_output_inferred": c_eff,
        "QE_C_eff_is_derived_physical_constant": False,
        "QE_C_eff_is_local_modulus": False,
        "historical_optimization_used_as_physics_input": False,
        "absolute_local_energy_density_scale_present": False,
        "microscopic_action_normalization_present": False,
        "microscopic_strain_coupling_derivative_present": False,
        "local_tensor_Hessian_present": False,
        "retarded_response_kernel_present": False,
        "coarse_graining_to_metric_strain_present": False,
        "physical_amplitude_bridge_complete": False,
        "physical_cluster_lensing_run_authorized": False,
        "replace_legacy_strength_now_authorized": False,
        "next_derivation_authorized": True,
        "next_derivation_target": next_interface["target"],
        "consistency": consistency,
        "thermal_identity_summary": thermal_result,
        "optimization_quarantine": quarantine,
        "gates": gates,
        "duration_seconds": time.perf_counter() - started,
    }
    write_json("validation.json", validation)
    write_json("repository_state.json", repo)

    print(json.dumps({
        "lab_id": LAB_ID,
        "outcome": outcome,
        "head_sha": repo["head_sha"],
        "QE_source_directory": "pbuf/data/quantum",
        "QE_alpha_qm": validation["QE_alpha_qm"],
        "QE_C_eff_output_inferred": c_eff,
        "physical_amplitude_bridge_complete": False,
        "next_derivation_target": validation["next_derivation_target"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
