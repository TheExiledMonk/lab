#!/usr/bin/env python3
"""PBUF FOUNDATION — QE AMPLITUDE CLOSURE 001 ROW-COUNT SEMANTICS FIX.

Narrow rerun of the raw Quantum Engine amplitude closure audit. The previous
run showed that all science/integrity gates passed except the row-count
consistency check. QE metadata defines 512 base temperature points while the
raw table contains 640 rows; config.fit_samples=128, so the generated table
contains base_grid_points + fit_sample_rows = 512 + 128 = 640.

No physics equation, tolerance, QE value, amplitude, or guardrail is changed.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

from pbuf.labs.foundation import quantum_engine_amplitude_closure001 as BASE

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "runs" / "quantum_engine_amplitude_closure001_row_count_fix"
LAB_ID = "PBUF-FOUNDATION-QUANTUM-ENGINE-AMPLITUDE-CLOSURE-001-ROW-COUNT-FIX"


def write_json(name: str, obj) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / name).write_text(json.dumps(obj, indent=2, sort_keys=False) + "\n")


def main() -> int:
    started = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)

    repo = {
        "branch": BASE.git("rev-parse", "--abbrev-ref", "HEAD"),
        "head_sha": BASE.git("rev-parse", "HEAD"),
        "tracked_changes": BASE.git("diff", "--name-only"),
        "staged_changes": BASE.git("diff", "--name-only", "--cached"),
        "stash_list": BASE.git("stash", "list"),
    }
    if repo["branch"] != "main":
        raise RuntimeError(f"runner must execute on main, got {repo['branch']}")
    if repo["tracked_changes"] or repo["staged_changes"]:
        raise RuntimeError("tracked or staged repository changes present")

    integ, payloads = BASE.direct_snapshot_integrity()
    write_json("snapshot_integrity.json", integ)
    if not integ["all_pass"]:
        raise RuntimeError("direct Quantum Engine snapshot integrity/provenance gate failed")

    config = json.loads(payloads["config.json"].decode("utf-8"))
    micro = json.loads(payloads["micro_cache.json"].decode("utf-8"))
    optimization = json.loads(payloads["optimisation_result.json"].decode("utf-8"))
    thermal = json.loads(payloads["thermal_table_cache.json"].decode("utf-8"))

    meta = thermal["metadata"]
    rows = thermal["rows"]

    base_temperature_points = int(config["temperature_points"])
    fit_samples = int(config["fit_samples"])
    metadata_base_points = int(meta["num_points"])
    micro_base_points = int(micro["temperature_points"])
    expected_total_rows = base_temperature_points + fit_samples
    observed_total_rows = len(rows)

    row_semantics = {
        "config_temperature_points": base_temperature_points,
        "micro_temperature_points": micro_base_points,
        "metadata_num_points": metadata_base_points,
        "config_fit_samples": fit_samples,
        "expected_total_rows_base_plus_fit_samples": expected_total_rows,
        "observed_total_rows": observed_total_rows,
        "base_temperature_point_metadata_match": bool(
            base_temperature_points == micro_base_points == metadata_base_points
        ),
        "total_row_count_match": bool(observed_total_rows == expected_total_rows),
        "interpretation": (
            "QE metadata num_points/temperature_points describe the 512-point base thermal grid; "
            "the serialized rows additionally contain the 128 configured fit-sample rows."
        ),
        "physics_change": False,
        "tolerance_change": False,
    }
    row_semantics["pass"] = bool(
        row_semantics["base_temperature_point_metadata_match"]
        and row_semantics["total_row_count_match"]
    )
    write_json("row_count_semantics.json", row_semantics)

    consistency = {
        "engine_source_match": micro["engine_source"] == meta["micro_source"] == "quantum_engine_v11",
        "method_version_match": int(config["method_version"]) == int(micro["method_version"]) == int(meta["method_version"]),
        "table_version_match": int(config["table_version"]) == int(micro["table_version"]) == int(meta["table_version"]),
        "base_temperature_point_metadata_match": row_semantics["base_temperature_point_metadata_match"],
        "serialized_total_row_count_match": row_semantics["total_row_count_match"],
        "alpha_qm_match": BASE.close(micro["alpha_qm"], meta["alpha_qm"]),
        "eps0_base_match": BASE.close(micro["eps0_base"], meta["eps0_base"]),
        "f_cut_match": BASE.close(micro["f_cut"], meta["f_cut_T"]),
        "f_coup_match": BASE.close(micro["f_coup"], meta["f_coup_T"]),
    }
    consistency["all_pass"] = all(consistency.values())

    thermal_result = BASE.thermal_audit(micro, thermal)
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
        "row_count_semantics_pass": row_semantics["pass"],
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
            "row_count_semantics": row_semantics,
            "optimization_quarantine": quarantine,
        })
        raise RuntimeError("Quantum Engine amplitude closure row-count-fix gate failed")

    outcome = "Outcome B — QE NORMALIZED AMPLITUDE CHAIN VERIFIED; ABSOLUTE LOCAL TENSOR RESPONSE NORMALIZATION STILL OPEN"
    validation = {
        "lab_id": LAB_ID,
        "outcome": outcome,
        "head_sha": repo["head_sha"],
        "QE_source_directory": "pbuf/data/quantum",
        "QE_storage_mode": "direct_raw_repository_files",
        "archive_or_base64_layer_used": False,
        "row_count_semantics_corrected": True,
        "row_count_semantics": row_semantics,
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
        "QE_alpha_qm": validation["QE_alpha_qm"],
        "QE_C_eff_output_inferred": c_eff,
        "row_count_semantics_pass": True,
        "physical_amplitude_bridge_complete": False,
        "next_derivation_target": validation["next_derivation_target"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
