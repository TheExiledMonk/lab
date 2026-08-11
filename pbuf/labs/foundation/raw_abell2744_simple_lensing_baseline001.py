#!/usr/bin/env python3
"""Dev160: fail-closed audit of the RAW Abell 2744 lensing baseline.

This lab deliberately does not turn detector counts into matter or a lens.  The
repository currently has no such transform.  It records that missing edge and
the separate, benchmark-assisted production chain without executing the latter.
"""
from __future__ import annotations

import ast
import json
import sys
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[3]
RAW = ROOT / "PBUF_raw_benchmark" / "WLRAW-001_Abell2744"
OUT = ROOT / "runs" / "raw_abell2744_simple_lensing_baseline001"


def dump(name: str, value: object) -> None:
    (OUT / name).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def _calls(path: Path) -> list[str]:
    tree = ast.parse(path.read_text())
    result = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            fn = node.func
            if isinstance(fn, ast.Name): result.append(fn.id)
            elif isinstance(fn, ast.Attribute): result.append(fn.attr)
    return sorted(set(result))


def raw_inventory() -> dict:
    selected_path = RAW / "provenance" / "selected_products.json"
    selected = json.loads(selected_path.read_text())
    products = selected["products"]
    classes = Counter(row["classification"] for row in products)
    examples = {}
    for row in products:
        cls = row["classification"]
        if cls in examples: continue
        folder = {"RAW_DETECTOR": "raw", "FLT_CALIBRATED_CONTROL": "flt",
                  "FLC_CTE_CALIBRATED_CONTROL": "flc"}[cls]
        path = RAW / folder / row["filename"]
        with fits.open(path, mode="readonly", memmap=False, lazy_load_hdus=True) as hdul:
            shapes = []
            for hdu in hdul:
                n = int(hdu.header.get("NAXIS", 0))
                shapes.append([int(hdu.header[f"NAXIS{i}"]) for i in range(n, 0, -1)] if n else None)
            examples[cls] = {
                "path": str(path.relative_to(ROOT)),
                "hdu_shapes": shapes,
                "bunit": [h.header.get("BUNIT") for h in hdul],
                "coordinate_metadata": {k: hdul[0].header.get(k) for k in
                    ("RA_TARG", "DEC_TARG", "PA_V3", "ORIENTAT")},
            }
    return {
        "target_cluster": "Abell2744", "data_mode": "RAW",
        "archive_root": str(RAW.relative_to(ROOT)), "archive_exists": RAW.is_dir(),
        "selection_manifest": str(selected_path.relative_to(ROOT)),
        "selection_sha256": selected["selection_sha256"],
        "selected_file_count": len(products), "classification_counts": dict(classes),
        "raw_exposure_count": classes.get("RAW_DETECTOR", 0),
        "fields": ["PRIMARY/SCI detector counts", "ERR", "DQ", "header metadata"],
        "representative_shapes_units_coordinates": examples,
        "raw_masks": "DQ extensions (archive input; not consumed by a lensing runner)",
        "raw_metadata": "FITS headers plus acquisition provenance manifests",
        "pbuf_internal_preprocessing": "NONE in acquisition; Dev126 performs read-only audit statistics",
        "raw_input_current_lensing_consumers": [],
        "classification_note": {
            "RAW_DETECTOR": "raw input", "FLT_CALIBRATED_CONTROL": "external calibrated control",
            "FLC_CTE_CALIBRATED_CONTROL": "external calibrated/CTE-corrected control",
            "DEV126_METRICS": "PBUF-derived audit intermediates, not lensing input"},
    }


def pipeline_inventory() -> dict:
    raw_consumers = [
        {"module": "pbuf.labs.data.frontier_abell2744_raw_acquisition001", "function": "main",
         "input": "MAST acquisition metadata/products", "output": "local RAW/FLT/FLC archive", "downstream": "Dev126/Dev134 audits"},
        {"module": "pbuf.labs.data.hst_acs_raw_flt_flc_audit001", "function": "main",
         "input": "RAW/FLT/FLC exposure families", "output": "read-only calibration audit", "downstream": "none in pbuf.wl"},
        {"module": "pbuf.labs.foundation.wl_hst_acs_detector_geometry001", "function": "main",
         "input": "FITS headers only", "output": "partial detector geometry inventory", "downstream": "no production mapping"},
    ]
    production = [
        {"stage": "D0", "module": "pbuf.core.benchmark_data", "function": "load_product",
         "input": "external/preprocessed Frontier Fields kappa/gamma FITS", "output": "2D observable arrays"},
        {"stage": "S0", "module": "pbuf.wl.source", "function": "load_cluster_source",
         "input": "five-cluster benchmark entry", "output": "benchmark-assisted rho3"},
        {"stage": "L0", "module": "pbuf.wl.native_response", "function": "build_native_response",
         "input": "rho3", "output": "M10 native vector response"},
        {"stage": "L0/P0", "module": "pbuf.wl.los", "function": "project_interface_to_los",
         "input": "M10 vector", "output": "Rx/Ry grid field"},
        {"stage": "P0", "module": "pbuf.wl.backends.cpu", "function": "CpuReferenceBackend.propagate",
         "input": "Rx/Ry field plus RayLaunch", "output": "ray checkpoints/final snapshot"},
        {"stage": "R0/O0", "module": "pbuf.labs.foundation._vulkan_g3d_common", "function": "downstream",
         "input": "launch plus propagation", "output": "screen, received state, decoded bank, candidates"},
    ]
    return {
        "raw_archive_consumers": raw_consumers, "separate_production_lensing_chain": production,
        "raw_to_native_source_edge": "ABSENT", "raw_to_native_lens_edge": "ABSENT",
        "active_raw_lensing_runner": None,
        "evidence": {
            "only_foundation_RAW_consumer": "pbuf/labs/foundation/wl_hst_acs_detector_geometry001.py",
            "geometry_runner_status": "production_mapping_performed=false; pixel access=false",
            "production_prepare_calls": _calls(ROOT / "pbuf/labs/foundation/_vulkan_g3d_common.py"),
            "production_source_import": "current_native_five_cluster_observable_benchmark001"},
        "five_cluster_benchmark_used_as_development_baseline": False,
        "five_cluster_chain_executed_by_dev160": False,
    }


def unavailable(stage: str, reason: str) -> dict:
    return {"stage": stage, "status": "NOT_DEFINED", "reason": reason}


def figure(name: str, title: str, lines: list[str]) -> None:
    fig, ax = plt.subplots(figsize=(10, 4.8)); ax.axis("off")
    ax.text(.5, .9, title, ha="center", va="top", fontsize=14, weight="bold")
    ax.text(.05, .72, "\n".join(lines), ha="left", va="top", family="monospace", fontsize=10)
    fig.tight_layout(); fig.savefig(OUT / name, dpi=120); plt.close(fig)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    raw = raw_inventory(); pipe = pipeline_inventory()
    missing = "No repository transform connects RAW detector pixels to rho3 or an M10 native lens field."
    dump("raw_input_inventory.json", raw)
    dump("active_pipeline_inventory.json", pipe)
    dump("native_source_support.json", {**unavailable("S0_NATIVE_SOURCE", missing),
        "source_native_support_cells": "NOT_DEFINED", "source_native_support_fraction": "NOT_DEFINED",
        "source_native_rms_radius": "NOT_DEFINED", "source_native_peak_location": "NOT_DEFINED"})
    dump("native_lens_support.json", {**unavailable("L0_NATIVE_LENS", missing),
        "lens_field_support_type": "NOT_DEFINED", "lens_field_domain_support": "NOT_DEFINED",
        "lens_field_rms_radius": "NOT_DEFINED", "lens_field_r50": "NOT_DEFINED", "lens_field_r68": "NOT_DEFINED",
        "lens_field_r90": "NOT_DEFINED", "lens_field_r95": "NOT_DEFINED", "lens_field_sign_structure": "NOT_DEFINED",
        "lens_field_anisotropy": "NOT_DEFINED", "positive_weight_manufactured": False})
    dump("propagation_object_contract.json", {
        "raw_lane_status": "NOT_REACHABLE", "current_production_propagation_object": "ZERO_WIDTH_RAY",
        "production_launch_container": "RAY_SHEET (independent zero-width rays on regular launch grid)",
        "per_element_state": ["x", "y", "z", "vx", "vy", "vz", "launch x0", "launch y0", "step index", "rx_sample", "ry_sample"],
        "implementation": "pbuf.labs.foundation.los_consistent_ray_geometry001._propagate_g3d",
        "dev159_state_used": False})
    dump("lens_sampling_contract.json", {"raw_lane_status": "NOT_REACHABLE",
        "current_production_method": "POINT_SAMPLE", "implementation": "_sample",
        "index_rule": "searchsorted then clipped preceding grid index", "interpolation": False,
        "finite_neighborhood": False, "integrated_footprint": False})
    dump("null_vs_loaded_response.json", {"status": "NOT_EXECUTED", "null_run": False, "loaded_run": False,
        "response_established": False, "reason": missing,
        "why_fail_closed": "Using benchmark rho3 would violate DATA_MODE=RAW; deriving rho3 from counts would add new lens/source physics."})
    dump("full_received_state_inventory.json", {**unavailable("R0_FULL_3D_RECEIPT", missing),
        "fields": [], "dimension": "NOT_DEFINED", "channel_count": 0, "support": "NOT_DEFINED",
        "production_checkpoint_fields_if_benchmark_lane_were_used": ["x", "y", "z", "vx", "vy", "vz", "rx_sample", "ry_sample"]})
    dump("projection_information_loss.json", {"status": "NOT_MEASURABLE_FOR_RAW_LANE", "J3": "NOT_DEFINED", "G3": "NOT_DEFINED",
        "J2": "NOT_DEFINED", "G2": "NOT_DEFINED", "delta_G": "NOT_DEFINED", "reason": missing})
    stages = {"D0_RAW": "EXPLICIT", "S0_NATIVE_SOURCE": "UNRESOLVED", "L0_NATIVE_LENS": "UNRESOLVED",
              "P0_PROPAGATION_HISTORY": "UNRESOLVED", "R0_FULL_3D_RECEIPT": "UNRESOLVED",
              "R1_PROJECTED_2D_RECEIPT": "UNRESOLVED", "O0_OBSERVER_PRODUCTS": "UNRESOLVED"}
    dump("reversal_stage_audit.json", {"lens_extent_information": stages, "R1_complete_history": "NOT_REACHABLE",
        "R2_full_3d": "NOT_REACHABLE", "R3_projected_2d": "NOT_REACHABLE", "R4_observer": "NOT_REACHABLE",
        "controlled_lens_extent_sweep": "NOT_JUSTIFIED", "truth_radius_used": False,
        "conclusion": "Extent recovery cannot be tested because no raw-derived lens extent exists upstream."})
    blocker = {"primary_reversal_blocker": "BLOCKER_UNRESOLVED", "exact_information_loss_point": "D0_TO_S0_EDGE_ABSENT",
        "old_propagation_lens_extent_identifiable": "UNRESOLVED", "reason": missing,
        "secondary_production_contract_observation": "When reached through the forbidden benchmark lane, rays use pointwise samples; this is not promoted to the RAW-lane blocker without a RAW run."}
    dump("blocker_classification.json", blocker)
    dump("dev159_relevance_contract.json", {"dev156_available": True, "dev157_available": True, "dev158_available": True,
        "dev159_available": True, "dev159_finite_state_used_in_lensing": False, "dev159_state_substituted_into_lensing": False,
        "dev159_expected_to_address_blocker": "UNRESOLVED", "reason": "Dev159 cannot supply the absent observational RAW-to-native source/lens bridge by itself."})
    matrix = {"RAW_ABELL2744_INPUT": "PARTIAL", "STATIC_NATIVE_LENS": "UNRESOLVED",
        "OLD_GEOMETRIC_PROPAGATION": "FROZEN_BASELINE", "OLD_RECEIVED_STATE": "UNRESOLVED",
        "OLD_OBSERVER": "UNRESOLVED", "DEV159_FINITE_PROPAGATING_STATE": "NOT_YET_USED"}
    dump("downstream_validity_matrix.json", matrix)
    contract = {"DEV160_AUDIT_COMPLETE": True, "TARGET_CLUSTER": "Abell2744", "DATA_MODE": "RAW",
        "RAW_ABELL2744_PIPELINE_LOCATED": "PARTIAL", "ACTIVE_RAW_LENSING_RUNNER_IDENTIFIED": False,
        "NATIVE_SOURCE_STATE_IDENTIFIED": False, "NATIVE_LENS_FIELD_IDENTIFIED": False,
        "LENS_FIELD_SUPPORT_MEASURABLE": False, "CURRENT_PROPAGATION_OBJECT": "ZERO_WIDTH_RAY",
        "LENS_SAMPLING_METHOD": "POINT_SAMPLE", "NULL_VS_LOADED_RESPONSE_ESTABLISHED": False,
        "FULL_3D_RECEIVED_STATE_AVAILABLE": False, "PROJECTED_2D_RECEIVED_STATE_AVAILABLE": False,
        "LENS_EXTENT_INFORMATION_AT_FULL_3D_RECEIPT": "UNRESOLVED",
        "LENS_EXTENT_INFORMATION_AFTER_2D_PROJECTION": "UNRESOLVED",
        "LENS_EXTENT_INFORMATION_AFTER_OBSERVER_REDUCTION": "UNRESOLVED",
        "OLD_PROPAGATION_LENS_EXTENT_IDENTIFIABLE": "UNRESOLVED", "PRIMARY_REVERSAL_BLOCKER": "BLOCKER_UNRESOLVED",
        "DEV159_EXPECTED_TO_ADDRESS_BLOCKER": "UNRESOLVED", "FIVE_CLUSTER_BENCHMARK_USED_AS_DEVELOPMENT_BASELINE": False,
        "DEV159_FINITE_STATE_USED_IN_LENSING": False, "DEV159_STATE_SUBSTITUTED_INTO_LENSING": False,
        "OBSERVER_MODIFIED": False, "LENS_PHYSICS_MODIFIED": False, "PROPAGATION_PHYSICS_MODIFIED": False,
        "OBSERVATIONAL_FIT_PERFORMED": False, "TARGET_DERIVED_RADIUS_USED": False,
        "ARBITRARY_LENS_RADIUS_INTRODUCED": False, "ARBITRARY_OBSERVER_KERNEL_INTRODUCED": False,
        "EXACT_FAILURE_POINT": "D0_RAW -> S0_NATIVE_SOURCE is not implemented"}
    dump("final_raw_a2744_lensing_contract.json", contract)
    figures = {
        "raw_source_support.png": ["D0: 116 RAW detector exposures present", "S0: NOT DEFINED (no RAW->rho3 transform)"],
        "native_lens_field.png": ["L0: NOT DEFINED for RAW lane", "No positive density or radius manufactured"],
        "lens_radial_support.png": ["RMS/R50/R68/R90/R95: NOT DEFINED"],
        "launch_sheet.png": ["P0: NOT REACHED", "Production contract: sheet of zero-width rays"],
        "received_3d_displacement.png": ["R0: NOT AVAILABLE for RAW lane"],
        "received_2d_displacement.png": ["R1: NOT AVAILABLE for RAW lane"],
        "propagation_delta_loaded_minus_null.png": ["C00/C01: NOT EXECUTED", "Frozen scope provides no valid RAW-derived loaded lens"],
        "information_loss_stage_summary.png": ["D0 RAW archive: EXPLICIT", "  |", "  X  missing RAW->native construction", "S0/L0/P0/R0/O0: UNRESOLVED"]}
    for name, lines in figures.items(): figure(name, name.removesuffix(".png").replace("_", " ").title(), lines)
    report = """DEV160 RAW ABELL 2744 SIMPLE-LENSING BASELINE AUDIT\n\nOutcome: fail-closed forensic baseline.\nExact failure point: D0_RAW -> S0_NATIVE_SOURCE is absent from the current repository.\nThe RAW archive is acquired and audited, but no active runner converts detector pixels into rho3 or a native M10 lens.\nThe actual production lensing chain instead imports current_native_five_cluster_observable_benchmark001. Dev160 did not execute it.\nConsequently C00/C01, lens support, receipt, projection loss, and reversal are UNRESOLVED for DATA_MODE=RAW.\nCreating those objects here would add a source/lens law and violate the frozen scope.\n\nPRIMARY_REVERSAL_BLOCKER=BLOCKER_UNRESOLVED\nEXACT_INFORMATION_LOSS_POINT=D0_TO_S0_EDGE_ABSENT\nFIVE_CLUSTER_BENCHMARK_USED_AS_DEVELOPMENT_BASELINE=false\nOBSERVER_MODIFIED=false\nLENS_PHYSICS_MODIFIED=false\nPROPAGATION_PHYSICS_MODIFIED=false\nDEV159_STATE_SUBSTITUTED_INTO_LENSING=false\n"""
    (OUT / "report.txt").write_text(report)
    print(report, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
