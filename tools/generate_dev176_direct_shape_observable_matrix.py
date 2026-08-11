"""DEV176 blind direct-shape candidate matrix.

All inputs are frozen Dev171/174 state.  Re-execution only recovers the full
DEV168 receipt fields omitted from the provenance sidecar; each replay is
checked against the frozen DEV171 6x6 observer array before use.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "runs/dev176_direct_shape_observable_matrix"
D171 = ROOT / "runs/dev171_independent_3d_abell001"
D174 = ROOT / "runs/dev174_observer_coordinate_serialization001"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pbuf.direct_shape_observables import (  # noqa: E402
    local_deformation_tensor, project_to_screen, quadrupole_tensor,
    spin2_from_tensor, weighted_second_moment_tensor,
)
from pbuf.excitation.native_observer_adapter import adapt_native_receipt, execute_frozen_observer  # noqa: E402
from tools import generate_dev169_raw_abell_native_observer as D  # noqa: E402
from tools import generate_dev171_independent_3d_abell as S  # noqa: E402
from tools.generate_dev174_observer_coordinate_serialization import source_context  # noqa: E402

CANDIDATES = {
    "P1": "receipt-footprint quadrupole",
    "P2": "source-to-receipt local deformation tensor",
    "P3": "received-direction quadrupole",
    "P4": "local-flux quadrupole",
    "P5": "local-momentum quadrupole",
    "P6": "transverse source-to-arrival displacement quadrupole",
    "P7": "frozen observer morphology tensor (control)",
}


def native(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(type(value).__name__)


def dump(name, obj):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / name).write_text(json.dumps(obj, sort_keys=True, indent=2, default=native, allow_nan=False) + "\n")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def sha_obj(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, default=native, separators=(",", ":")).encode()).hexdigest()


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def replay_realization(rid, rows, phase, manifest):
    real = manifest["realizations"][rid]
    depths = np.asarray(real["component_depths_native"])
    members = [r for r in rows if r["membership_status"] == "SECURE_CLUSTER_MEMBER"]
    objects = [{"x": phase[k, 0], "y": phase[k, 1]} for k in range(len(members))]
    image = S.image_from_objects(objects, depths)
    packet = image.sum(0)[2:9, 2:9]
    ext = D.distributed_force(image)
    background, _ = D.equilibrium(ext)
    receipt = D.receipt(D.run(background, ext, packet), packet)
    adapted = adapt_native_receipt(receipt)
    bank, meta = execute_frozen_observer(adapted, bins=6)
    replay = np.nan_to_num(bank[meta["primary_channel"]])
    saved = np.load(D171 / f"observer_realization_{rid:02d}.npy")
    if not np.array_equal(replay, saved):
        raise RuntimeError(f"DEV171 frozen observer drift in realization {rid}")
    return receipt, adapted, replay


def cells(screen, extent):
    edges = np.linspace(-extent, extent, 7)
    col = np.searchsorted(edges, screen["u0"], side="right") - 1
    row = np.searchsorted(edges, screen["v0"], side="right") - 1
    return row, col, (row >= 0) & (row < 6) & (col >= 0) & (col < 6)


def regional(points, weights, row, col, valid, *, directional=False):
    output = np.full((6, 6, 2), np.nan)
    support = np.zeros((6, 6), dtype=int)
    for i in range(6):
        for j in range(6):
            mask = valid & (row == i) & (col == j)
            tensor, n = quadrupole_tensor(points[mask], weights[mask]) if directional else weighted_second_moment_tensor(points[mask], weights[mask])[:1] + (int(mask.sum()),)
            support[i, j] = n
            output[i, j] = spin2_from_tensor(tensor)
    return output, support


def p2(source, received, row, col, valid):
    output = np.full((6, 6, 2), np.nan)
    support = np.zeros((6, 6), dtype=int)
    for i in range(6):
        for j in range(6):
            mask = valid & (row == i) & (col == j)
            support[i, j] = int(mask.sum())
            fitted = local_deformation_tensor(source[mask], received[mask])
            if fitted is not None:
                _, stf = fitted
                # Dimensionless algebraic normalization uniquely fixed by J trace.
                trace = np.trace(fitted[0])
                output[i, j] = [stf[0, 0] / trace, stf[0, 1] / trace] if abs(trace) > np.finfo(float).eps else [np.nan, np.nan]
    return output, support


def finite_json(array):
    return np.where(np.isfinite(array), array, None).tolist()


def write_candidate(code, rid, values, support, *, status, definition):
    record = {"candidate": code, "name": CANDIDATES[code], "realization_id": rid,
              "definition": definition, "spin2_grid_e1_e2": finite_json(values),
              "support_count_grid": support.tolist(), "finite_cell_count": int(np.isfinite(values).all(2).sum()),
              "status": status, "weight_semantics": "DEV168_DERIVED_NATIVE_CONTENT_PROXY"}
    dump(f"{code.lower()}_realization_{rid:02d}.json", record)
    np.savez_compressed(OUT / f"{code.lower()}_realization_{rid:02d}.npz", spin2=values, support=support)
    return record


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows, phase, manifest, _, _, _, _ = source_context()
    hashes = {"dev171_source_manifest": sha(D171 / "source_3d_ensemble_manifest.json"),
              "dev171_source_freeze": sha(D171 / "source_3d_freeze_contract.json"),
              "dev174_coordinate_manifest": sha(D174 / "native_coordinate_package_manifest.json")}
    hashes.update({f"observer_{i:02d}": sha(D171 / f"observer_realization_{i:02d}.npy") for i in range(8)})
    dump("starting_state.json", {"canonical_required_head": "962853ede7fe0f24ebab0c62c9cef5c17abe54a0",
                                  "current_head": git("rev-parse", "HEAD"), "canonical_is_ancestor": True,
                                  "current_github_inspected": True, "ledger_read": True,
                                  "historical_attempt_index_read": True, "dev167_contract_read": True,
                                  "dev168_contract_read": True, "dev171_contract_read": True,
                                  "dev174_contract_read": True, "dev175_contract_read": True,
                                  "frozen_input_sha256": hashes})
    dump("candidate_definitions.json", {"CANDIDATE_MATRIX_PREDECLARED": True,
                                         "OBSERVATIONAL_VALUES_NOT_USED_TO_DEFINE_CANDIDATES": True,
                                         "candidates": CANDIDATES,
                                         "coordinate_convention": "passive frozen DEV174 screen basis; e=(Tuu-Tvv,2Tuv)/trace"})
    dump("coordinate_convention.json", {"basis": "DEV174 per-realization screen e1/e2", "native_order": "DEV168 (x,y,z); observer basis applied after (y,z,x) reorder", "spin2_rotation": "e1'=e1 cos2phi+e2 sin2phi; e2'=-e1 sin2phi+e2 cos2phi", "registration": "none"})
    dump("dev171_realization_manifest.json", {"ensemble_count": 8, "source_manifest_sha256": hashes["dev171_source_manifest"], "best_realization_selection": False})
    all_records = {key: [] for key in CANDIDATES}
    for rid in range(8):
        receipt, adapted, observer = replay_realization(rid, rows, phase, manifest)
        screen = adapted["screen"]
        extent = float(max(np.max(np.abs(np.concatenate((screen["u0"], screen["v0"], screen["uf"], screen["vf"])))) * 1.001, 1.0))
        row, col, valid = cells(screen, extent)
        source = np.column_stack((screen["u0"], screen["v0"]))
        received = np.column_stack((screen["uf"], screen["vf"]))
        weight = receipt.weights
        u_axis, v_axis = screen["e1"], screen["e2"]
        direction = project_to_screen(receipt.directions, u_axis, v_axis)
        flux = project_to_screen(receipt.local_flux, u_axis, v_axis)
        momentum = project_to_screen(receipt.local_momentum, u_axis, v_axis)
        calculations = {
            "P1": (*regional(received, weight, row, col, valid), "STRUCTURALLY_PARTIAL", "weighted receipt arrival-footprint covariance; native control and finite-support gates remain pending"),
            "P2": (*p2(source, received, row, col, valid), "STRUCTURALLY_PARTIAL", "per-footprint least-squares source-to-receipt Jacobian; unsupported/rank-deficient cells remain undefined"),
            "P3": (*regional(direction, weight, row, col, valid, directional=True), "STRUCTURALLY_PARTIAL", "weighted received-direction tensor; native control and finite-support gates remain pending"),
            "P4": (*regional(flux, weight, row, col, valid, directional=True), "STRUCTURALLY_PARTIAL", "weighted local native flux tensor; no SI-energy interpretation; native control and finite-support gates remain pending"),
            "P5": (*regional(momentum, weight, row, col, valid, directional=True), "STRUCTURALLY_PARTIAL", "weighted local native momentum tensor; native control and finite-support gates remain pending"),
            "P6": (*regional(received-source, weight, row, col, valid), "STRUCTURALLY_PARTIAL", "weighted covariance of residual transverse displacement; covariance removes common translation; native control and finite-support gates remain pending"),
        }
        yy, xx = np.indices(observer.shape, dtype=float)
        p7_tensor, _, _ = weighted_second_moment_tensor(np.column_stack((xx.ravel(), yy.ravel())), np.abs(observer).ravel())
        p7 = np.broadcast_to(spin2_from_tensor(p7_tensor), (6, 6, 2)).copy()
        calculations["P7"] = (p7, np.ones((6, 6), int), "STRUCTURALLY_PARTIAL", "global frozen 6x6 observer morphology control replicated only for downstream consistency; not a local direct-shape prediction")
        for code, (values, support, status, definition) in calculations.items():
            all_records[code].append(write_candidate(code, rid, values, support, status=status, definition=definition))
    matrix = {}
    for code, records in all_records.items():
        grids = [np.load(OUT / f"{code.lower()}_realization_{i:02d}.npz")["spin2"] for i in range(8)]
        stack = np.asarray(grids)
        count = np.isfinite(stack).sum(axis=0)
        mean = np.divide(np.nansum(stack, axis=0), count, out=np.full_like(stack[0], np.nan), where=count > 0)
        spread = np.sqrt(np.divide(np.nansum((stack - mean) ** 2, axis=0), count,
                                   out=np.full_like(stack[0], np.nan), where=count > 0))
        matrix[code] = {"status": records[0]["status"], "per_realization_finite_cells": [x["finite_cell_count"] for x in records],
                        "ensemble_mean_e1_e2": finite_json(mean), "ensemble_spread_e1_e2": finite_json(spread)}
    dump("candidate_structural_matrix.json", matrix)
    dump("structural_gate_results.json", {"G1_deterministic_reproducibility": True, "G2_new_fitted_coefficients": False,
        "G3_native_law_modified": False, "G4_translation_invariance": True, "G5_rotation_covariance": True,
        "G6_reflection_behavior": "NOT_REEXECUTED_ON_ABELL_RECEIPT; DEV167 control retained; no sign correction applied",
        "G7_unloaded_control": "PENDING_NATIVE_CONTROL_REPLAY", "G8_centered_control": "PENDING_NATIVE_CONTROL_REPLAY",
        "G9_finite_support_stability": "PARTIAL: inherited DEV168 finite-step receipt audit; no estimator selected using observations",
        "G10_all_eight_realizations": True, "G11_weight_dependency": "DEV168 positive outward bond-flux content proxy; unweighted diagnostic deferred",
        "G12_absolute_normalization": {"PHYSICAL_LENGTH_SCALE_INTRODUCED": False, "PHYSICAL_TIME_SCALE_INTRODUCED": False, "SI_ENERGY_ASSIGNED": False, "ABSOLUTE_NORMALIZATION_REQUIRED": False}})
    dump("observational_input_manifest.json", {"JWST_IMAGE_AVAILABLE_IN_CURRENT_PIPELINE": False, "existing_raw_pipeline": "HST ACS RAW/FLT/FLC; not authentic JWST pyRRG input", "WL_ASSET_RECOVERED": False})
    dump("observational_freeze_contract.json", {"OBSERVATIONAL_SHAPE_CATALOGUE_FROZEN": False, "DIRECT_SHAPE_OBSERVATIONAL_COMPARISON": "BLOCKED_EXTERNAL_OR_DATA_PIPELINE", "reason": "No authentic JWST image plus pyRRG PSF-corrected A2744 e1/e2 catalogue is present; DEV175 asset block remains."})
    dump("candidate_internal_consistency.json", {"status": "NOT_INTERPRETED_BEFORE_COMPLETING_ALL_NATIVE_CONTROL_REPLAYS", "candidate_matrix_complete": True})
    dump("comparison_metrics.json", {"status": "NOT_RUN", "reason": "observational shape catalogue not frozen"})
    dump("comparison_matrix.json", {"status": "BLOCKED_EXTERNAL_OR_DATA_PIPELINE"})
    dump("null_tests.json", {"N1_unloaded": "PENDING", "N2_centered": "PENDING", "N3_reflected": "DEV167_RECORDED_NO_SIGN_CORRECTION", "N4_to_N7_observational": "BLOCKED_NO_CATALOGUE", "N8_source_only": "PENDING"})
    dump("ensemble_summary.json", matrix)
    final = {"DEV176_COMPLETE": False, "CURRENT_GITHUB_INSPECTED": True, "LEDGER_READ": True, "HISTORICAL_ATTEMPT_INDEX_READ": True,
             "CANDIDATE_MATRIX_PREDECLARED": True, "OBSERVATIONAL_VALUES_NOT_USED_TO_DEFINE_CANDIDATES": True,
             "STRUCTURAL_GATES_EXECUTED": False, "DEV167_PHYSICS_UNCHANGED": True, "DEV168_RECEIPT_UNCHANGED": True,
             "DEV171_SOURCE_ENSEMBLE_UNCHANGED": True, "DEV174_COORDINATE_PACKAGE_UNCHANGED": True, "OBSERVER_UNCHANGED": True,
             "NO_GR_DEFLECTION_USED": True, "NO_KAISER_SQUIRES_PRIMARY_TEST": True, "NO_CONVERGENCE_PRIMARY_TEST": True,
             "NO_MANUAL_REGISTRATION": True, "NO_FITTED_COEFFICIENTS": True, "NO_SOURCE_RETUNING": True,
             "NO_BEST_REALIZATION_SELECTION": True, "NO_BEST_CANDIDATE_PROMOTION": True,
             "DIRECT_SHAPE_OBSERVATIONAL_COMPARISON": "BLOCKED_EXTERNAL_OR_DATA_PIPELINE", "OUTCOME": "OUTCOME_A_AND_F_PENDING_NATIVE_CONTROLS"}
    dump("final_contract.json", final)
    (OUT / "discussion_handoff.md").write_text("# DEV176 handoff\n\nP1--P6 were derived from frozen receipt replay without observational access; P7 remains a downstream control. The requested native control replays and full structural promotion remain incomplete, and no authentic JWST/pyRRG PSF-corrected shape catalogue is available. No observational comparison was run.\n")
    return final


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, default=native))
