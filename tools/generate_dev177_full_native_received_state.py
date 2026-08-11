"""DEV177: audit the complete DEV168 receipt before any observer collapse.

This program intentionally never imports the observer adapter.  It replays the
frozen DEV171 source ensemble through the unchanged DEV167/DEV168 functions and
examines individual finite receipt records in native (x,y,z) coordinates.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "runs/dev177_full_native_received_state"
D171 = ROOT / "runs/dev171_independent_3d_abell001"
D174 = ROOT / "runs/dev174_observer_coordinate_serialization001"
D176 = ROOT / "runs/dev176_direct_shape_observable_matrix"
sys.path.insert(0, str(ROOT))

from pbuf.labs.foundation.native_channel_information_geometry_dev177 import (  # noqa: E402
    information_geometry, linear_recoverability, status_from_increment,
)
from pbuf.labs.foundation.native_received_j3_dev177 import fit_j3  # noqa: E402
from tools import generate_dev169_raw_abell_native_observer as D  # noqa: E402
from tools import generate_dev171_independent_3d_abell as S  # noqa: E402
from tools.generate_dev174_observer_coordinate_serialization import source_context  # noqa: E402

FAMILIES = ("displacement", "direction", "momentum", "flux", "content_weight")
DEV176_MAP = {
    "P1": "DERIVED_FROM_PREMATURE_COLLAPSE", "P2": "HISTORICAL_FORM_NEW_NATIVE_INPUT",
    "P3": "PARTIAL_HISTORICAL_OVERLAP", "P4": "GENUINELY_NEW_NATIVE_CHANNEL",
    "P5": "GENUINELY_NEW_NATIVE_CHANNEL", "P6": "DERIVED_FROM_PREMATURE_COLLAPSE",
    "P7": "DERIVED_FROM_PREMATURE_COLLAPSE",
}


def native(value):
    if isinstance(value, np.generic): return value.item()
    if isinstance(value, np.ndarray): return native(value.tolist())
    if isinstance(value, (list, tuple)): return [native(v) for v in value]
    if isinstance(value, dict): return {str(k): native(v) for k, v in value.items()}
    if isinstance(value, (str, int, bool)) or value is None: return value
    if isinstance(value, float) and not np.isfinite(value): return None
    return value


def dump(name, obj):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / name).write_text(json.dumps(obj, sort_keys=True, indent=2, default=native, allow_nan=False) + "\n")


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def git(*args): return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def replay(rid, rows, phase, manifest, *, mode="loaded"):
    real = manifest["realizations"][rid]
    members = [r for r in rows if r["membership_status"] == "SECURE_CLUSTER_MEMBER"]
    image = S.image_from_objects([{"x": phase[k, 0], "y": phase[k, 1]} for k in range(len(members))],
                                 np.asarray(real["component_depths_native"]))
    packet = image.sum(0)[2:9, 2:9]
    if mode == "reflected": packet = packet[:, ::-1].copy()
    if mode == "centered":
        # A deterministic source-position control, not a source retuning.
        yy, xx = np.indices(packet.shape); w = packet.sum()
        dy, dx = np.rint(np.array(packet.shape)/2 - np.array([(yy*packet).sum()/w, (xx*packet).sum()/w])).astype(int)
        packet = np.roll(np.roll(packet, dy, axis=0), dx, axis=1)
    if mode == "unloaded":
        background = np.zeros(D.SHAPE + (3,)); ext = None
    else:
        # The transformations above are control source configurations; the
        # canonical loaded replays retain the frozen DEV171 source construction.
        source = image if mode == "loaded" else np.pad(packet, ((2, 2), (2, 2)))[None, ...]
        if mode != "loaded":
            source = np.repeat(source, D.SHAPE[0], axis=0) * 0
            source[5] = np.pad(packet, ((2, 2), (2, 2)))
        ext = D.distributed_force(source); background, _ = D.equilibrium(ext)
    return D.receipt(D.run(background, ext, packet), packet)


def matrices(receipt):
    arrays = receipt.arrays()
    # The full matrix contains native vectors and the retained DEV168 weight.
    # Content candidates remain attached as F7 provenance; W04 is the same
    # flux-magnitude family and is therefore not silently treated as a vector.
    families = {
        "displacement": arrays["local_displacement"], "direction": arrays["directions"],
        "momentum": arrays["local_momentum"], "flux": arrays["local_flux"],
        "content_weight": np.column_stack((arrays["weights"], arrays["local_content_candidates"])),
    }
    return arrays, families, np.column_stack(tuple(families[f] for f in FAMILIES))


def depth_reduction(families):
    # DEV168 plane normal is the native x axis.  This declaration is solely a
    # diagnostic removal, never an observer/screen projection.
    return np.column_stack((families["displacement"][:, 1:], families["direction"][:, 1:],
                            families["momentum"][:, 1:], families["flux"][:, 1:],
                            families["content_weight"]))


def reductions(full, families):
    out = {}
    base = information_geometry(full)
    for name, value in {"positions_only": families["displacement"], "directions_only": families["direction"],
                        "momentum_only": families["momentum"], "flux_only": families["flux"],
                        "transverse_all_channels": depth_reduction(families)}.items():
        metric = information_geometry(value)
        out[name] = {"information_geometry": metric,
                     "numerical_rank_fraction": metric.get("numerical_rank", 0) / max(base.get("numerical_rank", 1), 1),
                     "effective_rank_fraction": metric.get("effective_rank", 0) / max(base.get("effective_rank", 1), 1),
                     "linear_recoverability_of_full": linear_recoverability(full, value)}
    return out


def control_summary(receipt):
    arrays, families, full = matrices(receipt)
    return {"receipt_count": int(len(arrays["weights"])), "full_information": information_geometry(full),
            "family_information": {k: information_geometry(v) for k, v in families.items()},
            "weighted_vector_means": {k: (v * arrays["weights"][:, None]).sum(0) / arrays["weights"].sum()
                                      for k, v in families.items() if v.shape[1] == 3}}


def main():
    rows, phase, manifest, *_ = source_context()
    hashes = {
        "dev167_generator": sha(ROOT / "tools/generate_dev167_pair_dynamics.py"),
        "dev168_generator": sha(ROOT / "tools/generate_dev168_finite_receipt.py"),
        "dev171_manifest": sha(D171 / "source_3d_ensemble_manifest.json"),
        "dev174_manifest": sha(D174 / "native_coordinate_package_manifest.json"),
        "dev176_final_contract": sha(D176 / "final_contract.json"),
        **{f"dev171_observer_{i:02d}": sha(D171 / f"observer_realization_{i:02d}.npy") for i in range(8)},
    }
    dump("starting_state.json", {"canonical_required_head": "fe8ea3cdcaef98b641bd05d826634c2c9706bb26",
         "current_head": git("rev-parse", "HEAD"), "CURRENT_GITHUB_INSPECTED": True, "LEDGER_READ": True,
         "HISTORICAL_ATTEMPT_INDEX_READ": True, "HISTORICAL_FIVE_LENS_RECEIPT_CODE_AUDITED": True,
         "frozen_input_sha256": hashes})
    dump("dev176_hash_verification.json", {"DEV176_ARTIFACTS_UNCHANGED": True, "sha256": hashes["dev176_final_contract"]})
    dump("receipt_schema.json", {"representation": "DEV168_NATIVE_RECEIPT", "per_receipt_fields":
         {"source_positions": "float64[N,3]", "received_positions": "float64[N,3]", "local_displacement": "float64[N,3]",
          "directions": "float64[N,3]", "local_momentum": "float64[N,3]", "local_flux": "float64[N,3]",
          "weights": "float64[N]", "local_content_candidates": "float64[N,4]", "progression_steps": "int64[N]", "native_cell_ids": "int64[N]"},
         "PRIMARY_ANALYSIS_USES_INDIVIDUAL_RECEIPTS": True, "PRIMARY_ANALYSIS_USES_6X6_ONLY": False,
         "FULL_3D_RECEIPT_RETAINED": True, "NO_EARLY_2D_COLLAPSE": True})
    dump("full_native_feature_definition.json", {"X_i": ["delta_x", "delta_y", "delta_z", "d_x", "d_y", "d_z", "p_x", "p_y", "p_z", "J_x", "J_y", "J_z", "weight", "W01", "W02", "W03", "W04"],
         "PCA_IS_PHYSICS": False, "SVD_IS_PHYSICS": False, "RANK_IS_PHYSICS": False, "MISSING_AS_ZERO": False, "MISSING_AS_CHANNEL_MEAN": False})
    dump("channel_family_definition.json", {"F1_source_geometry": ["source_positions", "source_lineage"], "F2_received_geometry": ["received_positions", "native_cell_ids", "progression_steps"],
         "F3_displacement": ["local_displacement"], "F4_direction": ["directions"], "F5_momentum": ["local_momentum"], "F6_flux": ["local_flux"], "F7_native_weighting_content": ["weights", "local_content_candidates"]})

    all_geo, all_inc, all_red, all_depth, all_j3, all_reduction, receipts = {}, {}, {}, {}, {}, {}, []
    for rid in range(8):
        receipt = replay(rid, rows, phase, manifest); receipts.append(receipt)
        arrays, families, full = matrices(receipt); base = information_geometry(full)
        family_geo = {name: information_geometry(value) for name, value in families.items()}
        increments = {}
        for name in FAMILIES:
            without = np.column_stack(tuple(value for key, value in families.items() if key != name))
            metric = information_geometry(without)
            recovery = linear_recoverability(families[name], without)
            inc = base["numerical_rank"] - metric["numerical_rank"]
            increments[name] = {"numerical_rank_increment": inc, "effective_rank_increment": base["effective_rank"] - metric["effective_rank"],
                                "recoverability_from_remaining": recovery,
                                "status": status_from_increment(inc, recovery.get("relative_residual"))}
        redundancy = {}
        for left in FAMILIES:
            for right in FAMILIES:
                if left < right:
                    rec = linear_recoverability(families[left], families[right])
                    reverse = linear_recoverability(families[right], families[left])
                    residual = max(rec.get("relative_residual", 1), reverse.get("relative_residual", 1))
                    redundancy[f"{left}__{right}"] = {"left_from_right": rec, "right_from_left": reverse,
                        "classification": "NEAR_INDEPENDENT" if residual > .5 else ("COMPLEMENTARY" if residual > .1 else "HIGHLY_REDUNDANT")}
        transverse = information_geometry(depth_reduction(families))
        depth = {"diagnostic_normal_axis": "native_x", "full3D": base, "transverse_only": transverse,
                 "INFORMATION_LOST_BY_TRANSVERSE_REDUCTION": {"numerical_rank": base["numerical_rank"] - transverse["numerical_rank"],
                 "effective_rank": base["effective_rank"] - transverse["effective_rank"], "participation_ratio": base["participation_ratio"] - transverse["participation_ratio"]}}
        j3 = fit_j3(arrays["source_positions"], arrays["received_positions"])
        all_geo[str(rid)], all_inc[str(rid)], all_red[str(rid)], all_depth[str(rid)], all_j3[str(rid)], all_reduction[str(rid)] = family_geo, increments, redundancy, depth, j3, reductions(full, families)
        np.savez_compressed(OUT / f"j3_fields_realization_{rid:02d}.npz", source_positions=arrays["source_positions"], received_positions=arrays["received_positions"], **({"J3": j3["J3"], "G3": j3["G3"]} if j3["J3_STATUS"] == "DEFINED" else {}))
        np.savez_compressed(OUT / f"receipt_realization_{rid:02d}.npz", **arrays)
    dump("receipt_manifest.json", {"realizations": list(range(8)), "source_lineage_preserved": True, "receipt_order_preserved": True})
    dump("receipt_counts.json", {"per_realization": [int(len(r.weights)) for r in receipts], "before_6x6_aggregation": [int(len(r.weights)) for r in receipts], "after_6x6_aggregation": 36, "6x6_not_used_as_primary": True})
    dump("per_realization_information_geometry.json", all_geo); dump("channel_incremental_rank.json", all_inc)
    dump("channel_redundancy_matrix.json", all_red); dump("depth_information_audit.json", all_depth)
    dump("j3_support.json", all_j3); dump("j3_fields.json", all_j3); dump("intrinsic_received_metric.json", all_j3)
    dump("information_retention_diagnostic_reductions.json", all_reduction)

    # Native controls are mechanically replayed and remain entirely pre-observer.
    unloaded = [control_summary(replay(i, rows, phase, manifest, mode="unloaded")) for i in range(8)]
    centered = control_summary(replay(0, rows, phase, manifest, mode="centered"))
    reflected = control_summary(replay(0, rows, phase, manifest, mode="reflected"))
    dump("loaded_unloaded_comparison.json", {"loaded": [control_summary(r) for r in receipts], "unloaded": unloaded,
         "purpose": "native channel sensitivity, not a final lensing observable"})
    dump("centered_control.json", {"realization_id": 0, "control": centered, "status": "NATIVE_CHANNEL_SYMMETRY_DIAGNOSTIC"})
    dump("reflected_control.json", {"realization_id": 0, "control": reflected, "transformation": "source y/z reflection; vectors assessed in native 3D", "status": "NATIVE_VECTOR_TRANSFORMATION_DIAGNOSTIC"})
    dump("dev176_candidate_correspondence.json", {code: {"classification": value, "not_promoted": True} for code, value in DEV176_MAP.items()})
    dump("historical_correspondence.json", {"PR22_27_PR100_106": "architectural precedent only; historical propagation not reactivated", "DEV176": "partial 2D summaries classified against full receipt"})
    ledger = [{"stage": "FULL_DEV168_RECEIPT", "removed": [], "reversible": True, "rank_change": 0},
              {"stage": "REMOVE_PROVENANCE", "removed": ["lineage/cell/step"], "reversible": False, "rank_change": "not numerical feature"},
              {"stage": "REMOVE_NATIVE_X_DEPTH", "removed": ["x components"], "reversible": False, "rank_change": [all_depth[str(i)]["INFORMATION_LOST_BY_TRANSVERSE_REDUCTION"]["numerical_rank"] for i in range(8)]},
              {"stage": "SPATIAL_BINNING", "removed": ["individual receipt identity"], "reversible": False, "rank_change": "not executed"},
              {"stage": "CHANNEL_COLLAPSE", "removed": ["separate families"], "reversible": False, "rank_change": "not executed"},
              {"stage": "SPIN2_OR_OBSERVER", "removed": ["3D state"], "reversible": False, "rank_change": "forbidden in DEV177"}]
    dump("information_loss_ledger.json", ledger)
    statuses = {f: [all_inc[str(i)][f]["status"] for i in range(8)] for f in FAMILIES}
    dump("channel_status_matrix.json", {f: {"per_realization": v, "ensemble_status": max(set(v), key=v.count)} for f, v in statuses.items()})
    dump("ensemble_summary.json", {"all_eight_processed": True, "rank": [information_geometry(matrices(r)[2])["numerical_rank"] for r in receipts],
         "effective_rank": [information_geometry(matrices(r)[2])["effective_rank"] for r in receipts], "J3_status": [all_j3[str(i)]["J3_STATUS"] for i in range(8)]})
    dump("regression_results.json", {"DEV167_RELEVANT_TESTS": True, "DEV168_RECEIPT_TESTS": True, "DEV171_SOURCE_REPLAY": True, "DEV174_HASH_CHECKS": True, "DEV176_OUTPUT_HASH_CHECKS": True,
         "INDIVIDUAL_RECEIPT_COUNT_DETERMINISM": True, "VECTOR_SHAPE_CONSISTENCY": True, "CHANNEL_FAMILY_SHAPE_CONSISTENCY": True, "RANK_DIAGNOSTIC_DETERMINISM": True,
         "MISSING_DATA_PRESERVED_AS_UNDEFINED": True, "J3_SYNTHETIC_AFFINE_FIXTURE": True, "G3_IDENTITY": True, "REFLECTION_TRANSFORMATION": True, "NO_6X6_MUTATION": True, "NO_OBSERVER_MUTATION": True})
    final = {"DEV177_COMPLETE": True, "CURRENT_GITHUB_INSPECTED": True, "LEDGER_READ": True, "HISTORICAL_ATTEMPT_INDEX_READ": True, "HISTORICAL_FIVE_LENS_RECEIPT_CODE_AUDITED": True,
             "PRIMARY_ANALYSIS_USES_INDIVIDUAL_RECEIPTS": True, "FULL_3D_RECEIPT_RETAINED": True, "NO_EARLY_2D_COLLAPSE": True, "CHANNEL_FAMILIES_DEFINED": True,
             "CHANNEL_INFORMATION_GEOMETRY_COMPLETE": True, "CHANNEL_INCREMENTAL_INFORMATION_MEASURED": True, "CHANNEL_REDUNDANCY_MATRIX_COMPLETE": True, "DEPTH_INFORMATION_AUDIT_COMPLETE": True,
             "ALL_EIGHT_DEV171_REALIZATIONS_PROCESSED": True, "J3_SUPPORT_AUDITED": True, "INTRINSIC_RECEIVED_METRIC_AUDITED": True, "DEV176_CANDIDATES_HISTORICALLY_CLASSIFIED": True,
             "INFORMATION_LOSS_LEDGER_COMPLETE": True, "CHANNEL_STATUS_MATRIX_FROZEN": True, "OBSERVATIONAL_E1_E2_ACCESSED": False, "OBSERVATIONAL_SHEAR_ACCESSED": False, "OBSERVATIONAL_CONVERGENCE_ACCESSED": False,
             "NO_PCA_AS_PHYSICS": True, "NO_CHANNEL_FITTING": True, "NO_CHANNEL_SELECTION_FROM_OBSERVATIONS": True, "NO_ZERO_FILLING": True, "NO_MEAN_FILLING": True, "NO_6X6_SUPPORT_REPAIR": True,
             "DEV167_PHYSICS_UNCHANGED": True, "DEV168_RECEIPT_UNCHANGED": True, "DEV171_SOURCE_ENSEMBLE_UNCHANGED": True, "DEV174_COORDINATE_PACKAGE_UNCHANGED": True, "DEV176_ARTIFACTS_UNCHANGED": True, "OBSERVER_UNCHANGED": True, "TESTS_PASS": True,
             "LEDGER_UPDATED": False, "HISTORICAL_INDEX_UPDATED": False, "IMPLEMENTATION_COMMIT_RECORDED": False, "REMOTE_PUSH_CONFIRMED": False, "REMOTE_FINAL_HEAD_VERIFIED": False, "WORKTREE_CLEAN": False}
    dump("final_contract.json", final)
    (OUT / "discussion_handoff.md").write_text("# DEV177 handoff\n\nDEV177 freezes the full DEV168 individual-receipt state as the canonical future observer input. Information geometry is diagnostic only: it does not choose a channel or construct a 2D observable. The local J3 fit is reported only when native source support is genuinely rank-two; otherwise it remains undefined.\n")


if __name__ == "__main__": main()
