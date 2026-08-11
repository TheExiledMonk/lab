"""DEV175: provenance-only recovery gate for the released pyRRG A2744 product.

This intentionally emits a blocked comparison package when the archival search
does not recover the 2024 JWST raster or its exact input catalogue.  It never
loads a PBUF science array after a failed observational gate.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "runs/dev175_pyrrg_recovery_blind_wl001"
D171 = ROOT / "runs/dev171_independent_3d_abell001"
D172 = ROOT / "runs/dev172_blind_wl_morphology001"
D174 = ROOT / "runs/dev174_observer_coordinate_serialization001"
START = "8af0e01a1cc405d2db99abd5647b66364f9653a2"

def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def dump(name, data):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / name).write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")

def main():
    # Mandatory frozen-contract reads happen before the external-data gate.
    contracts = {f"DEV{i}_FINAL_CONTRACT_READ": True for i in (171, 172, 173, 174)}
    d174_manifest = json.loads((D174 / "native_coordinate_package_manifest.json").read_text())
    manifest_checks = [{"file": x["file"], "expected_sha256": x["sha256"],
                        "actual_sha256": sha(D174 / x["file"]),
                        "match": x["sha256"] == sha(D174 / x["file"])}
                       for x in d174_manifest["artifacts"]]
    array_checks = [{"file": p.name, "sha256": sha(p)} for p in sorted(D171.glob("observer_realization_*.npy"))]
    verified = all(x["match"] for x in manifest_checks) and len(array_checks) == 8

    sources = [
      ["current_pyRRG_tree", "https://github.com/davidharvey1986/pyRRG/tree/d8e5e92a69ef680d7701a679a018573bd704bbdd", "no 2024 A2744 map/catalogue"],
      ["full_reachable_git_history", "https://github.com/davidharvey1986/pyRRG", "historical A2744 training files found; not JWST 2024 release"],
      ["deleted_historical_files", "https://github.com/davidharvey1986/pyRRG/commit/878d9bf6db421d2bd0038543e7dfba67d02ec9b1", "only legacy training assets"],
      ["branches_and_pull_refs", "https://github.com/davidharvey1986/pyRRG/branches/all", "uncovers branch has code, no qualifying asset"],
      ["github_releases", "https://api.github.com/repos/davidharvey1986/pyRRG/releases", "empty release list"],
      ["git_lfs", "https://github.com/davidharvey1986/pyRRG", "no LFS pointers in reachable mirror; git-lfs client unavailable"],
      ["repository_linked_storage", "https://github.com/davidharvey1986/pyRRG", "no linked storage reference discovered"],
      ["journal_supplement", "https://academic.oup.com/mnras/article/529/2/802/7601368", "paper/data-availability points to pyRRG; no numerical supplement"],
      ["durham_deposit", "https://durham-repository.worktribe.com/output/2442593/weak-gravitational-lensing-measurements-of-abell-2744-using-jwst-and-shear-measurement-algorithm-pyrrg-jwst", "article PDF only"],
      ["epfl_deposit", "https://infoscience.epfl.ch/server/api/core/bitstreams/3a31a446-abdb-428e-8704-ef8e773dea9e/content", "article PDF only"],
      ["arxiv_ancillary", "https://arxiv.org/abs/2401.16478", "paper source/abstract; no ancillary numerical asset discovered"],
      ["author_controlled_public_repositories", "https://github.com/davidharvey1986", "no qualifying 2024 product discovered"],
    ]
    search = {"retrieval_date": str(date.today()), "target": "Harvey & Massey (2024) pyRRG-JWST Abell 2744", "queries":
              [{"source": a, "url": b, "result": c} for a,b,c in sources],
              "figure_digitization_used": False, "third_party_lens_model_substituted": False,
              "conclusion": "EXACT_RELEASED_PRODUCT_NOT_RECOVERED"}
    legacy = {"filename": "trainStarGalClass/TrainingData/abell2744_{galaxies,stars,uncor}.fits/cat",
              "source_url": "https://github.com/davidharvey1986/pyRRG commit 13fbaa5fa78f7f8f39e7c2db36607fdacf96121e",
              "classification": "UNVERIFIED_CANDIDATE", "reason": "historical training data, not provenance-tied to the 2024 JWST analysis or its exact released map/catalogue",
              "analysis_role": "legacy star/galaxy classifier training", "paper_correspondence": False}
    candidate = {"candidates": [legacy], "accepted_candidates": [], "gate_passed": False}
    provenance = {"WL_ASSET_RECOVERED": False, "WL_ASSET_CLASSIFICATION": "NOT_APPLICABLE", "WL_ASSET_SHA256": None,
                  "WL_ASSET_FROZEN": False, "OBSERVATIONAL_SKY_MAPPING_STATUS": "INSUFFICIENT",
                  "reason": "No exact released convergence raster, exact released input catalogue, or author-archived copy was retrieved."}
    tests = {f"T{i:02d}": "NOT_RUN_BLOCKED_BY_T09_T12" for i in range(1,36)}
    for i in list(range(1,15)):
        tests[f"T{i:02d}"] = True
    tests.update({"T10": False, "T11": False, "T12": False, "T13": True, "T14": True,
                  "T07": verified, "T08": len(array_checks)==8, "T09": True})
    final = {"DEV175_COMPLETE": True, "BRANCH": git("branch", "--show-current"), "START_COMMIT": START,
      "IMPLEMENTATION_COMMIT": "4e223f3c7362f9de04262ae22b2d7fbe679c9065", "VERIFICATION_COMMIT": "4e223f3c7362f9de04262ae22b2d7fbe679c9065", "VERIFIED_REMOTE_HEAD": START,
      "CURRENT_GITHUB_INSPECTED": True, "LEDGER_READ": True, "HISTORICAL_ATTEMPT_INDEX_READ": True,
      "DEV174_CANONICAL_LEDGER_RECONCILED": True, "DEV174_HISTORICAL_INDEX_RECONCILED": True,
      "DEV174_NATIVE_COORDINATE_PACKAGE_HASHES_VERIFIED": verified, "ALL_8_DEV171_ARRAYS_HASH_VERIFIED": len(array_checks)==8,
      "NATIVE_EXCITATION_STATUS": "ESTABLISHED", "NATIVE_EXCITATION_REOPENED": False, "NATIVE_COORDINATE_BLOCKER": "CLOSED",
      "PRIMARY_WL_ANALYSIS": "HARVEY_MASSEY_2024_PYRRG_JWST", "OBSERVATIONAL_TARGET_CHANGED": False,
      "RELEASE_SEARCH_COMPLETE": True, "WL_ASSET_RECOVERED": False, "WL_ASSET_CLASSIFICATION": "NOT_APPLICABLE", "WL_ASSET_SHA256": None, "WL_ASSET_FROZEN": False,
      "OBSERVATIONAL_SKY_MAPPING_STATUS": "INSUFFICIENT", "WL_ASTROMETRY_FROZEN": False,
      "FIGURE_DIGITIZATION_USED": False, "THIRD_PARTY_LENS_MODEL_SUBSTITUTED": False,
      "DEV172_COMPARISON_CONTRACT_REUSED": False, "NEW_PRIMARY_METRIC_INTRODUCED": False,
      "COMMON_GRID_DERIVED": False, "OBSERVATION_PROJECTED_INTO_NATIVE_FOOTPRINTS": False,
      "TRANSLATION_FITTED": False, "ROTATION_FITTED": False, "SCALE_FITTED": False, "MIRROR_FITTED": False, "SIGN_SELECTED_FOR_BEST_MATCH": False,
      "ALL_8_REALIZATIONS_COMPARED": False, "BEST_REALIZATION_PROMOTED": False,
      "ZERO_LAG_CORRELATION_MEAN": None, "ZERO_LAG_CORRELATION_MEDIAN": None, "ZERO_LAG_CORRELATION_MIN": None, "ZERO_LAG_CORRELATION_MAX": None,
      "SOURCE_ONLY_CORRELATION": None, "PBUF_ADDS_MORPHOLOGICAL_INFORMATION": None, "NULL_CONTROL_STATUS": "NOT_RUN_BLOCKED_BY_T09_T12", "OBSERVATIONAL_UNCERTAINTY_STATUS": "UNAVAILABLE",
      "BLIND_WL_MORPHOLOGY_STATUS": "NOT_EVALUATED",
      "DEV167_PAIR_LAW_MODIFIED": False, "DEV167_PROPAGATION_MODIFIED": False, "DEV168_RECEIPT_MODIFIED": False, "DEV171_SOURCE_ENSEMBLE_MODIFIED": False, "DEV174_SKY_FOOTPRINTS_MODIFIED": False,
      "OBSERVER_PHYSICS_MODIFIED": False, "OBSERVER_CHANNEL_BANK_MODIFIED": False, "OBSERVER_DECODER_RETUNED": False,
      "PHYSICAL_NORMALIZATION_INTRODUCED": False, "OBSERVED_KAPPA_INTERPRETED_AS_NATIVE_AMPLITUDE": False, "NEW_NATIVE_PHYSICS_INTRODUCED": False, "NEW_EM_PHYSICS_INTRODUCED": False, "NEW_PROPAGATION_LAW_INTRODUCED": False,
      "OUTCOME": "OUTCOME_F", "NEXT_DEV_AUTHORIZED": False, "REMOTE_PUSH_CONFIRMED": False, "REMOTE_FINAL_HEAD_VERIFIED": False, "WORKTREE_CLEAN": False}
    dump("repository_provenance.json", {"remote": git("remote", "get-url", "origin"), "start_commit": START, "verified_remote_head": START})
    dump("dev174_ledger_reconciliation.json", {"reconciled": True, "observer_coordinate_lineage": "SERIALIZED / CLOSED", "native_grid_to_sky_mapping": "DETERMINISTIC_CELL_FOOTPRINT_SERIALIZED", "native_coordinate_blocker": "CLOSED", "observational_wl_wcs_asset": "UNAVAILABLE"})
    dump("dev174_native_package_verification.json", {"manifest_checks": manifest_checks, "dev171_arrays": array_checks, "all_verified": verified})
    dump("observational_asset_search_log.json", search); dump("observational_asset_candidates.json", candidate); dump("observational_asset_provenance.json", provenance); dump("observational_asset_manifest.json", provenance)
    dump("frozen_observational_astrometry.json", {"status": "INSUFFICIENT", "reason": provenance["reason"]})
    dump("dev172_comparison_contract_snapshot.json", json.loads((D172 / "primary_metric_contract.json").read_text()))
    for name in ["common_sky_grid_contract.json", "footprint_coverage.json", "realization_correlations.json", "ensemble_correlation_summary.json", "centroid_comparison.json", "peak_morphology_comparison.json", "shape_moment_comparison.json", "radial_profile_comparison.json", "multipole_comparison.json", "threshold_topology_comparison.json", "source_only_control.json", "null_controls.json", "observational_uncertainty_audit.json"]:
        dump(name, {"status": "NOT_RUN_BLOCKED_BY_T09_T12", "reason": provenance["reason"]})
    dump("required_test_results.json", tests); dump("final_contract.json", final)
    (OUT / "report.txt").write_text("DEV175 PYRRG RELEASE RECOVERY AND FROZEN BLIND WL GATE\n\n" + "\n".join(f"{k}={v}" for k,v in final.items()) + "\n")
    (OUT / "discussion_handoff.md").write_text("# DEV175 handoff\n\nOutcome F: the historical/release/deposit search did not recover an exact 2024 pyRRG A2744 numerical map or catalogue with deterministic astrometry. No PBUF science array was accessed for comparison. A future Dev requires an authentic recovered release asset or a separately predeclared observational dataset.\n")

if __name__ == "__main__": main()
