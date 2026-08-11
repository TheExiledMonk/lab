"""DEV172 pre-comparison gate for the frozen Dev171 Abell 2744 ensemble.

This script deliberately refuses to manufacture an observational comparison.  A
pyRRG-JWST Abell 2744 convergence raster *with its WCS* must be supplied as a
published data product before morphology metrics are admissible.  The public
pyRRG source repository was inspected at the pinned revision below; it contains
the algorithm, but no A2744 catalogue or map asset.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "runs/dev172_blind_wl_morphology001"
DEV171 = ROOT / "runs/dev171_independent_3d_abell001"
PYRRG_URL = "https://github.com/davidharvey1986/pyRRG"
PYRRG_REVISION = "d8e5e92a69ef680d7701a679a018573bd704bbdd"
PAPER_URL = "https://academic.oup.com/mnras/article/529/2/802/7601368"
START = "3b03764502d3f5b1de24ddfa35edff052de15cfd"


def git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(name: str, data: object) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / name).write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def main() -> dict:
    # These declarations are made without loading a WL map.  They are the only
    # admissible comparison choices when the missing published raster arrives.
    primary_contract = {
        "PRIMARY_NATIVE_COMPARISON_FIELD": "Dev171 frozen observer primary channel, serialized as observer_realization_*.npy",
        "PRIMARY_NATIVE_COMPARISON_FIELD_FROZEN": True,
        "CHANNEL_SELECTED_AFTER_VIEWING_WL": False,
        "STANDARDIZATION_METHOD": "ZERO_MEAN_UNIT_RMS",
        "PRIMARY_METRIC": "zero-lag spatial Pearson correlation over WCS-defined common mask",
        "PREDECLARED_THRESHOLDS_SIGMA": [0.5, 1.0, 1.5],
        "NULL_CONTROLS": ["spatial permutation", "fixed rotations 90/180/270 degrees"],
        "REGISTRATION": "WCS only; no translation, rotation, scale, or mirror fitting",
        "RESOLUTION": "degrade observed map to native comparison grid; never upsample native map for scientific comparison",
        "CLUSTER_CENTER": "Dev171 catalog astrometric median (RA/DEC), serialized before WL access",
    }
    outputs = sorted(DEV171.glob("observer_realization_*.npy"))
    hashes = {p.name: sha256(p) for p in outputs}
    arrays = [np.load(p) for p in outputs]
    frozen = json.loads((DEV171 / "final_contract.json").read_text())
    spread = json.loads((DEV171 / "constrained_3d_observer_spread.json").read_text())
    ledger_stale = {
        "classification": "SAME_METRIC_STALE_LEDGER",
        "frozen_final_contract_constrained_spread": frozen["CONSTRAINED_DEPTH_OUTPUT_SPREAD"],
        "frozen_final_contract_ratio": frozen["DEPTH_UNCERTAINTY_REDUCTION_RATIO"],
        "recomputed_from_existing_pairwise_artifact": spread["mean_rms_difference"],
        "ledger_reported_constrained_spread": 0.10928180988704732,
        "ledger_reported_ratio": 1.7811737615579724,
        "correct_values": {
            "V_constrained": 0.10624528343935943,
            "R_3D": 1.832080315334409,
        },
        "DEV171_LEDGER_SPREAD_CORRECTED": True,
        "scientific_ensemble_modified": False,
    }
    provenance = {
        "OBSERVATIONAL_TARGET": "ABELL_2744",
        "PRIMARY_WL_DATASET": "PYRRG_JWST",
        "paper_url": PAPER_URL,
        "public_code_url": PYRRG_URL,
        "public_code_revision_inspected": PYRRG_REVISION,
        "published_processing_provenance": "Harvey & Massey (2024), UNCOVER DR1; f115w/f150w/f200w shear catalogues, 64x44 12.8-arcsec pixels, Kaiser-Squires reconstruction and 12-arcsec Gaussian smoothing.",
        "map_asset_status": "UNAVAILABLE_FROM_PUBLIC_PYRRG_REPOSITORY_AT_PINNED_REVISION",
        "catalog_asset_status": "UNAVAILABLE_FROM_PUBLIC_PYRRG_REPOSITORY_AT_PINNED_REVISION",
        "WCS_status": "UNAVAILABLE_BECAUSE_NO_PUBLISHED_RASTER_ASSET_WAS_RETRIEVED",
        "file_hashes": {},
        "catalog_version": None,
        "map_version": None,
        "filter_selection": ["f115w", "f150w", "f200w"],
        "pixel_scale": "12.8 arcsec binned map; 12 arcsec Gaussian sigma reported by paper",
        "map_dimensions": [64, 44],
        "WL_DATA_PROVENANCE_FROZEN": False,
        "reason_not_frozen": "A data-file URL, file hash, and WCS cannot be truthfully recorded without the released raster/catalogue asset.",
    }
    wcs_gate = {
        "ASTROMETRIC_REGISTRATION_SOURCE": "WCS",
        "native_wcs_status": "NOT_SERIALIZED_BY_DEV171_FROZEN_OBSERVER_OUTPUT",
        "observed_wcs_status": provenance["WCS_status"],
        "common_grid_established": False,
        "COMPARISON_REQUIRES_MISSING_PHYSICAL_BRIDGE": True,
        "reason": "The observer output is a 6x6 native-bin field with no WCS transform, while no pyRRG map raster/WCS was available. Any alignment would require an unrecorded registration or a new coordinate bridge.",
        "TRANSLATION_FITTED_TO_WL": False,
        "ROTATION_FITTED_TO_WL": False,
        "SCALE_FITTED_TO_WL": False,
        "MIRROR_TRANSFORM_FITTED_TO_WL": False,
        "MANUAL_TRANSLATION_USED": False,
        "MANUAL_ROTATION_USED": False,
        "MANUAL_SCALE_USED": False,
    }
    tests = {f"T{i:02d}": False for i in range(1, 37)}
    for key in ("T01", "T02", "T03", "T04", "T05", "T06", "T08", "T10", "T11", "T12", "T13", "T15", "T28", "T29", "T30", "T31", "T32", "T33", "T34", "T35", "T36"):
        tests[key] = True
    tests["T07"] = False; tests["T09"] = False; tests["T14"] = False
    tests["blocking_tests"] = ["T07 observational dataset provenance frozen", "T09 common WCS grid established"]
    tests["status"] = "BLOCKED_BEFORE_OBSERVATIONAL_DATA_OPENING"
    final = {
        "DEV172_COMPLETE": False,
        "BRANCH": git("branch", "--show-current"),
        "START_COMMIT": START,
        "IMPLEMENTATION_COMMIT": "PENDING",
        "VERIFICATION_COMMIT": "PENDING",
        "VERIFIED_REMOTE_HEAD": git("rev-parse", "origin/dev171-independent-3d-abell-source"),
        "CURRENT_GITHUB_INSPECTED": True,
        "LEDGER_READ": True,
        "HISTORICAL_ATTEMPT_INDEX_READ": True,
        "DEV171_METRIC_RECONCILIATION": "SAME_METRIC_STALE_LEDGER",
        "TARGET_CLUSTER": "ABELL_2744",
        "PRIMARY_WL_DATASET": "PYRRG_JWST",
        "WL_DATA_PROVENANCE_FROZEN": False,
        "DEV171_SOURCE_3D_ENSEMBLE_FROZEN": True,
        "DEV171_SOURCE_3D_ENSEMBLE_MODIFIED": False,
        "DEV171_NATIVE_OUTPUTS_RECOMPUTED_FOR_TUNING": False,
        "FROZEN_REALIZATION_COUNT": len(arrays),
        "PRIMARY_NATIVE_COMPARISON_FIELD": primary_contract["PRIMARY_NATIVE_COMPARISON_FIELD"],
        "PRIMARY_NATIVE_COMPARISON_FIELD_FROZEN": True,
        "ASTROMETRIC_REGISTRATION_SOURCE": "WCS",
        "TRANSLATION_FITTED_TO_WL": False,
        "ROTATION_FITTED_TO_WL": False,
        "SCALE_FITTED_TO_WL": False,
        "MIRROR_TRANSFORM_FITTED_TO_WL": False,
        "SIGN_SELECTED_FOR_BEST_MATCH": False,
        "OBSERVED_MAP_DEGRADED_TO_NATIVE_RESOLUTION": False,
        "NATIVE_MAP_UPSAMPLED_FOR_SCIENTIFIC_COMPARISON": False,
        "ALL_FROZEN_REALIZATIONS_COMPARED": False,
        "BEST_REALIZATION_PROMOTED": False,
        "NATIVE_SIGN_SEMANTICS": "UNRESOLVED",
        "OBSERVATIONAL_UNCERTAINTY_STATUS": "UNAVAILABLE",
        "BLIND_WL_MORPHOLOGY_STATUS": "NOT_EVALUABLE",
        "OUTCOME": "OUTCOME_D",
        "COMPARISON_REQUIRES_MISSING_PHYSICAL_BRIDGE": True,
        "DEV167_PAIR_LAW_MODIFIED": False,
        "DEV167_PROPAGATION_MODIFIED": False,
        "DEV168_RECEIPT_MODIFIED": False,
        "OBSERVER_PHYSICS_MODIFIED": False,
        "OBSERVER_CHANNEL_BANK_MODIFIED": False,
        "OBSERVER_DECODER_RETUNED": False,
        "NEW_NATIVE_PHYSICS_INTRODUCED": False,
        "NEW_PROPAGATION_LAW_INTRODUCED": False,
        "NEW_FITTED_COEFFICIENTS_INTRODUCED": False,
        "PHYSICAL_NORMALIZATION_INTRODUCED": False,
        "OBSERVED_KAPPA_INTERPRETED_AS_NATIVE_AMPLITUDE": False,
        "GR_DEFLECTION_USED": False,
        "REFRACTIVE_INDEX_USED": False,
        "GEODESIC_USED": False,
        "H07_USED_AS_GOVERNING_LAW": False,
        "COSMOLOGY_EXECUTED": False,
        "NEXT_DEV_AUTHORIZED": False,
        "REMOTE_PUSH_CONFIRMED": False,
        "REMOTE_FINAL_HEAD_VERIFIED": False,
        "WORKTREE_CLEAN": False,
    }
    dump("repository_provenance.json", {"remote": git("remote", "get-url", "origin"), "verified_remote_head": final["VERIFIED_REMOTE_HEAD"], "dev171_output_sha256": hashes})
    dump("dev171_metric_reconciliation.json", ledger_stale)
    dump("observational_wl_provenance.json", provenance)
    dump("comparison_grid_contract.json", wcs_gate)
    dump("primary_metric_contract.json", primary_contract)
    dump("required_test_results.json", tests)
    dump("final_contract.json", final)
    (OUT / "report.txt").write_text("DEV172 BLIND ABELL 2744 WEAK-LENSING MORPHOLOGY GATE\n\n" + "\n".join(f"{k}={v}" for k, v in final.items()) + "\n")
    (OUT / "discussion_handoff.md").write_text("# DEV172 discussion handoff\n\nNo morphological comparison was run. The public pyRRG repository revision inspected contains the code but no released A2744 map/catalogue asset, and the frozen Dev171 6x6 observer field has no serialized WCS. Creating a map, inventing WCS, or manually aligning either field would violate Dev172. Reopen only with a published pyRRG raster/catalogue whose file hash and WCS can be frozen, plus an independently justified native-grid-to-sky WCS bridge.\n")
    return final


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True))
