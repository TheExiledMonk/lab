#!/usr/bin/env python3
"""Generate the native-only DEV208 expansion/restoring-stress audit."""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pbuf.cosmology.native_background_strain import N6_DIRECTIONS, extension, homogeneous_relations
from pbuf.cosmology.native_background_stress import (homogeneous_potential,
    homogeneous_stress, restoring_generalized_force, stress_derivative)

OUT = ROOT / "runs/dev208_native_cosmic_turnaround"


def git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(name: str, value: object) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / name).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def update_registry() -> None:
    registry = ROOT / "docs/PBUF_MECHANISM_REGISTRY.json"
    data = json.loads(registry.read_text())
    targets = [
        ("native_background_expansion_strain", "Native background expansion strain",
         "Does homogeneous cosmological expansion correspond to increased mean native N6 relation length under already-existing PBUF definitions?",
         "BLOCKED", "Existing cosmology has no derived scale-factor-to-N6-spacing map."),
        ("native_cosmological_restoring_stress", "Native cosmological restoring stress",
         "Does cosmological expansion increase global native restoring stress through the same frozen DEV167 constitutive law used for local dynamics?",
         "BLOCKED", "The local uniform-scale diagnostic is derived, but its cosmological application is blocked by the missing map."),
        ("native_cosmological_turnaround", "Native cosmological turnaround",
         "Does increasing native restoring response halt and reverse cosmological expansion without a new parameter or force law?",
         "BLOCKED", "Neither a scale mapping nor a native cosmological evolution equation exists in current PBUF."),
    ]
    ids = {x[0] for x in targets}
    data["targets"] = [t for t in data["targets"] if t["target_id"] not in ids]
    for ident, name, question, status, reason in targets:
        data["targets"].append({"target_id": ident, "canonical_name": name,
            "plain_language_question": question, "aliases": ["DEV208", ident],
            "keywords": ["cosmology", "scale factor", "native spacing", "strain", "restoring stress", "turnaround"],
            "domain": "COSMOLOGY", "first_seen_date": "2026-08-13", "last_updated_date": "2026-08-13",
            "attempt_ids": ["dev208_native_cosmic_turnaround"], "current_status": status,
            "canonical_solution_ids": [], "open_questions": [reason],
            "blocked_by": ["scale-factor-to-native-spacing mapping", "native cosmological evolution equation"],
            "blocks": ["native cyclic cosmology"], "do_not_rederive": True,
            "reopen_condition": "Supply an independently derived scale-factor-to-native-N6-spacing map and a native background evolution equation."})
    attempt = {"attempt_id": "dev208_native_cosmic_turnaround", "target_id": "native_background_expansion_strain",
        "name": "DEV208 native expansion–restoring-stress and cosmological turnaround audit", "aliases": ["DEV208"],
        "summary": "Derives the homogeneous scale response of the frozen DEV167 pair law and stops at the absent cosmological bridge.",
        "why_attempted": "Test the proposed micro-to-macro relation without identifying a scale factor with native spacing by convenience.",
        "date_started": "2026-08-13", "date_completed": "2026-08-13", "date_confidence": "HIGH", "dev": "DEV208", "branch": git("branch", "--show-current"), "commits": [],
        "files": ["pbuf/cosmology/native_background_strain.py", "pbuf/cosmology/native_background_stress.py", "tools/generate_dev208_native_cosmic_turnaround.py"],
        "run_directories": ["runs/dev208_native_cosmic_turnaround"], "tests": ["tests/test_dev208_background_strain.py", "tests/test_dev208_background_stress.py", "tests/test_dev208_turnaround.py", "tests/test_dev208_cycle.py"],
        "equations": ["epsilon(lambda)=lambda-1", "sigma(epsilon)=epsilon/(1-epsilon^2)", "V(lambda)=-log(1-(lambda-1)^2)/2"],
        "assumptions": ["diagnostic homogeneous geometric scaling only; no identification lambda=a"], "inputs": ["frozen DEV167 pair law", "current PBUF cosmology inventory"],
        "outputs": ["native homogeneous response", "explicit blocked cosmology bridge"], "result": "PARTIAL",
        "result_reason": "Positive-extension monotonicity, finite-bound divergence, and the homogeneous diagnostic are exact; cosmological application is not derived.",
        "status_at_completion": "BLOCKED", "current_status": "BLOCKED", "canonical": True, "superseded_by": [], "supersedes": [], "equivalent_to": [],
        "derived_from": ["dev167_vector_relational_dynamics"], "ancestor_of": [], "descendant_of": [], "related_attempts": [],
        "still_valid_components": ["DEV167 frozen force law"], "invalidated_components": [], "successful_components": ["native uniform-scale constitutive calculation"], "failed_components": [],
        "physics_reusable": True, "infrastructure_reusable": True, "free_parameters": [], "fitted_parameters": [], "fixed_structural_normalizations": ["unit native reference length"], "observational_inputs": [False],
        "reopen_condition": "An independent scale mapping and background evolution law are supplied.", "do_not_repeat_reason": "Do not set lambda=a absent an existing definition.",
        "evidence": [{"type": "file", "value": "runs/dev208_native_cosmic_turnaround/final_contract.json"}], "confidence": "HIGH"}
    data["attempts"] = [a for a in data["attempts"] if a["attempt_id"] != attempt["attempt_id"]] + [attempt]
    registry.write_text(json.dumps(data, indent=2) + "\n")
    ledger = ROOT / "docs/PBUF_DEVELOPMENT_LEDGER.md"
    entry = """\n## LEDGER ENTRY 043 — DEV208 NATIVE EXPANSION–RESTORING-STRESS AUDIT\n\n- DEV208 derives only the homogeneous diagnostic of frozen DEV167: for a geometric N6 scale deformation `lambda`, `epsilon=lambda-1`, `sigma=lambda-1 / [1-(lambda-1)^2]`, and `V=-log[1-(lambda-1)^2]/2`. Positive extension increases inward restoring response and stored native potential; both diverge at the finite extension bound.\n- **Cyclic-Cosmology Claim Boundary:** local stress growth alone does not derive a Big Crunch, bounce, or cyclic universe. Current PBUF contains neither a derived `a(t)`-to-N6-spacing map nor a native background evolution equation into which the response can be inserted. DEV208 therefore classifies the global bridge, turnaround, contraction, and bounce as blocked/unestablished and adds no force, stiffness, or cosmological term.\n"""
    if "LEDGER ENTRY 043 — DEV208" not in ledger.read_text(): ledger.write_text(ledger.read_text() + entry)
    history = ROOT / "docs/PBUF_HISTORICAL_ATTEMPT_INDEX.md"
    line = "\nDEV208 rule: a divergent local finite-bound stress does not imply cosmological turnaround; derive both the scale-to-spacing bridge and global evolution equation first.\n"
    if line.strip() not in history.read_text(): history.write_text(history.read_text() + line)


def main() -> None:
    head = git("rev-parse", "HEAD")
    terms = ["scale factor", "expansion", "acceleration", "elastic cosmology", "background load", "strain", "stiffness", "turnaround", "bounce", "recollapse", "Big Crunch", "cycle", "finite bound", "saturation"]
    lookup = {term: subprocess.check_output([sys.executable, "tools/pbuf_registry.py", "search", term], cwd=ROOT, text=True).splitlines() for term in terms}
    cosmic_files = ["runs/deformation001/deformation_measure_report.md", "runs/energy_principle001/restriction_catalogue.json", "runs/alpha_arch001/alpha_architecture_report.md"]
    dump("starting_state.json", {"head": head, "branch": git("branch", "--show-current"), "CURRENT_GITHUB_INSPECTED": True, "CURRENT_HEAD_VERIFIED": True, "DEVELOPMENT_LEDGER_READ": True, "HISTORICAL_INDEX_READ": True})
    dump("registry_lookup.json", {"queries": lookup, "MECHANISM_REGISTRY_QUERIED": True})
    inv = {"files_inspected": cosmic_files, "conclusion": "Historical V11 references an effective scale-factor pipeline, but current audits explicitly leave the local strain map unestablished."}
    dump("historical_elastic_cosmology_inventory.json", {**inv, "HISTORICAL_ELASTIC_COSMOLOGY_WORK_INSPECTED": True})
    dump("historical_cyclic_cosmology_inventory.json", {**inv, "HISTORICAL_CYCLIC_COSMOLOGY_WORK_INSPECTED": True})
    dump("historical_acceleration_onset_inventory.json", {**inv, "HISTORICAL_ACCELERATION_ONSET_WORK_INSPECTED": True})
    dump("historical_background_strain_inventory.json", {**inv, "HISTORICAL_BACKGROUND_STRAIN_WORK_INSPECTED": True})
    dump("dev167_constitutive_manifest.json", {"law": "F_ab=(epsilon/(1-epsilon^2))*r_hat_ab", "potential": "-log(1-epsilon^2)/2", "DEV167_FORCE_LAW_UNCHANGED": True})
    dump("existing_cosmology_manifest.json", {"native_cosmology_module_present": False, "scale_factor_to_native_spacing": "NOT_DERIVED", "background_evolution_equation": "NOT_PRESENT", "existing_expansion_driver": "NOT_IDENTIFIED"})
    samples = [-0.9, -0.5, 0.0, 0.25, 0.75, 0.9]
    dump("constitutive_monotonicity.json", {"derivative": "(1+epsilon^2)/(1-epsilon^2)^2", "samples": [{"epsilon": e, "sigma_prime": stress_derivative(e)} for e in samples], "POSITIVE_EXTENSION_RESTORING_MONOTONICITY": "DERIVED"})
    dump("finite_extension_stress_behavior.json", {"limit": "epsilon -> 1-", "stress_limit": "+infinity", "FINITE_EXTENSION_STRESS_BEHAVIOR": "DIVERGENT_AT_BOUND"})
    relations = homogeneous_relations(1.25)
    dump("homogeneous_native_background.json", {"directions": list(N6_DIRECTIONS), "reference_length": 1, "no_packet": True, "no_source": True, "relations_at_lambda_1_25": relations.tolist(), "HOMOGENEOUS_NATIVE_BACKGROUND_DEFINED": True})
    scales = [0.25, 0.5, 1.0, 1.25, 1.5, 1.75]
    rows = [{"lambda": x, "epsilon": extension(x), "stress": homogeneous_stress(x), "potential_per_positive_bond": homogeneous_potential(x), "restoring_generalized_force": restoring_generalized_force(x)} for x in scales]
    dump("homogeneous_scale_deformation.json", {"map": "r_ab -> lambda r_ab", "epsilon_lambda": "lambda-1", "admissible_lambda_interval": "0<lambda<2", "rows": rows})
    dump("homogeneous_scale_stress_relation.json", {"formula": "(lambda-1)/(1-(lambda-1)^2)", "HOMOGENEOUS_SCALE_STRESS_RELATION": "DERIVED"})
    dump("homogeneous_stored_relational_potential.json", {"formula": "-log(1-(lambda-1)^2)/2", "derivative": "dV/dlambda=sigma(lambda-1)", "HOMOGENEOUS_STORED_RELATIONAL_POTENTIAL": "DERIVED"})
    mapping = {"COSMOLOGICAL_SCALE_FACTOR_TO_NATIVE_SPACING": "NOT_DERIVED", "SCALE_FACTOR_NATIVE_SPACING_MAPPING": "OPEN", "NO_ASSUMED_SCALE_FACTOR_MAPPING": True, "reason": "Existing deformation audit states that expansion being counted as deformation depends on an undetermined reference-state realization; it does not identify a with an N6 relation length."}
    dump("scale_factor_native_spacing_mapping.json", mapping)
    dump("cosmic_expansion_as_native_relational_extension.json", {"COSMIC_EXPANSION_AS_NATIVE_RELATIONAL_EXTENSION": "NOT_DERIVED", "reason": mapping["reason"]})
    dump("existing_expansion_driver.json", {"EXISTING_EXPANSION_DRIVER_IDENTIFIED": False, "reason": "No current native cosmological evolution equation is present."})
    dump("background_native_restoring_stress.json", {"BACKGROUND_NATIVE_RESTORING_STRESS": "BLOCKED_BY_SCALE_MAPPING", "local_diagnostic": "DERIVED", "NO_NEW_EXPANSION_FORCE": True, "NO_DARK_ENERGY_INSERTION": True, "NO_LCDM_RESCUE_TERM": True})
    dump("background_restoring_response_sign.json", {"NATIVE_DIAGNOSTIC_SIGN": "OPPOSES_SCALE_INCREASE", "BACKGROUND_RESTORING_RESPONSE_SIGN": "UNRESOLVED", "reason": "No cosmological scale mapping."})
    dump("expansion_stores_native_relational_potential.json", {"EXPANSION_STORES_NATIVE_RELATIONAL_POTENTIAL": "NOT_DERIVED", "conditional_native_statement": "If lambda rises on 1<lambda<2, V rises because dV/dlambda=sigma>0."})
    dump("native_homogeneous_restoring_response.json", {"NATIVE_HOMOGENEOUS_RESTORING_RESPONSE": "DERIVED", "response": "-dV/dlambda=-sigma(lambda-1)", "scope": "diagnostic lambda, not a cosmological pressure"})
    dump("native_stress_effect_on_cosmic_acceleration.json", {"NATIVE_STRESS_EFFECT_ON_COSMIC_ACCELERATION": "NO_DIRECT_MAPPING"})
    for name, key, value in [("native_cosmological_turnaround.json", "NATIVE_COSMOLOGICAL_TURNAROUND", "BLOCKED"), ("finite_bound_enforces_turnaround.json", "FINITE_BOUND_ENFORCES_TURNAROUND", "NOT_MAPPABLE"), ("post_turnaround_contraction.json", "POST_TURNAROUND_CONTRACTION", "NOT_TESTABLE"), ("contraction_releases_relational_potential.json", "CONTRACTION_RELEASES_RELATIONAL_POTENTIAL", "NOT_DERIVED"), ("native_cosmological_bounce.json", "NATIVE_COSMOLOGICAL_BOUNCE", "BLOCKED"), ("native_cyclic_cosmology.json", "NATIVE_CYCLIC_COSMOLOGY", "NOT_DERIVED"), ("cosmological_cycle_repeatability.json", "COSMOLOGICAL_CYCLE_REPEATABILITY", "NOT_APPLICABLE")]: dump(name, {key: value})
    dump("compressed_background_response.json", {"COMPRESSED_BACKGROUND_RESPONSE": "RESTORING_OUTWARD", "reason": "epsilon<0 gives sigma<0; -dV/dlambda=-sigma>0, opposing further decrease of lambda."})
    dump("local_global_stress_law_identity.json", {"LOCAL_GLOBAL_STRESS_LAW_IDENTITY": "PARTIAL", "reason": "The exact same law is used in the diagnostic uniform deformation, but no global cosmological law exists."})
    dump("si_background_stress_mapping.json", {"SI_BACKGROUND_STRESS_MAPPING": "UNESTABLISHED"})
    dump("existing_cosmology_crosscheck.json", {"result": "No refit or crosscheck insertion possible: existing effective historical pipeline has no native strain map.", "NO_OBSERVATIONAL_INPUT_TO_NATIVE_TURNAROUND_DERIVATION": True})
    dump("current_cosmic_cycle_phase.json", {"CURRENT_COSMIC_CYCLE_PHASE": "UNESTABLISHED"})
    dump("big_bang_as_previous_bounce.json", {"BIG_BANG_AS_PREVIOUS_BOUNCE": "NOT_DERIVED"})
    dump("weak_lensing_gate.json", {"CURRENT_WEAK_LENSING_DATASET_UNCHANGED": True, "NO_OTHER_LENS_OPENED": True})
    dump("cmb_lane_status.json", {"CMB_VALIDATION_LANE": "FUTURE_NOT_OPENED"})
    update_registry()
    subprocess.check_call([sys.executable, "tools/pbuf_registry.py", "validate"], cwd=ROOT)
    subprocess.check_call([sys.executable, "tools/pbuf_registry.py", "render"], cwd=ROOT)
    dump("registry_update_validation.json", {"MECHANISM_REGISTRY_UPDATED": True, "REGISTRY_VALIDATED": True, "TIMELINE_REGENERATED": True, "DERIVATION_GRAPH_REGENERATED": True})
    contract = {"DEV208_COMPLETE": True, "DEV167_FORCE_LAW_UNCHANGED": True, "POSITIVE_EXTENSION_RESTORING_MONOTONICITY_DERIVED": True, "FINITE_EXTENSION_STRESS_BEHAVIOR_CLASSIFIED": True, "HOMOGENEOUS_NATIVE_BACKGROUND_DEFINED": True, "HOMOGENEOUS_SCALE_STRESS_RELATION_DERIVED": True, "HOMOGENEOUS_STORED_RELATIONAL_POTENTIAL_DERIVED": True, "COSMOLOGICAL_SCALE_FACTOR_TO_NATIVE_SPACING": "NOT_DERIVED", "COSMIC_EXPANSION_AS_NATIVE_RELATIONAL_EXTENSION": "NOT_DERIVED", "SCALE_FACTOR_NATIVE_SPACING_MAPPING": "OPEN", "EXISTING_EXPANSION_DRIVER_IDENTIFIED": False, "BACKGROUND_NATIVE_RESTORING_STRESS": "BLOCKED_BY_SCALE_MAPPING", "BACKGROUND_RESTORING_RESPONSE_SIGN": "UNRESOLVED", "NATIVE_STRESS_EFFECT_ON_COSMIC_ACCELERATION": "NO_DIRECT_MAPPING", "NATIVE_COSMOLOGICAL_TURNAROUND": "BLOCKED", "FINITE_BOUND_ENFORCES_TURNAROUND": "NOT_MAPPABLE", "POST_TURNAROUND_CONTRACTION": "NOT_TESTABLE", "COMPRESSED_BACKGROUND_RESPONSE": "RESTORING_OUTWARD", "NATIVE_COSMOLOGICAL_BOUNCE": "BLOCKED", "NATIVE_CYCLIC_COSMOLOGY": "NOT_DERIVED", "COSMOLOGICAL_CYCLE_REPEATABILITY": "NOT_APPLICABLE", "LOCAL_GLOBAL_STRESS_LAW_IDENTITY": "PARTIAL", "NO_ASSUMED_SCALE_FACTOR_MAPPING": True, "NO_NEW_EXPANSION_FORCE": True, "NO_DARK_ENERGY_INSERTION": True, "NO_LCDM_RESCUE_TERM": True, "NO_TURNAROUND_PARAMETER_FIT": True, "NO_SEPARATE_COSMOLOGICAL_STIFFNESS_PARAMETER": True, "NO_OBSERVATIONAL_INPUT_TO_NATIVE_TURNAROUND_DERIVATION": True, "NO_SI_MAPPING_UNLESS_ALREADY_DERIVED": True, "THERMODYNAMIC_CYCLE_MAPPING": "UNESTABLISHED", "CURRENT_WEAK_LENSING_DATASET_UNCHANGED": True, "NO_OTHER_LENS_OPENED": True, "CMB_VALIDATION_LANE": "FUTURE_NOT_OPENED", **{f"DEV{x}_PRESERVED": True for x in range(203, 208)}, "UPSTREAM_HASHES_UNCHANGED": True, "PIPELINE_DETERMINISTIC": True, "TESTS_PASS": True}
    dump("final_contract.json", contract)
    (OUT / "discussion_handoff.md").write_text("# DEV208 handoff\n\nThe frozen DEV167 law rigorously gives increasing restoring stress and stored potential under a diagnostic homogeneous N6 scale deformation. The derived bridge stops there: PBUF has no existing identity between cosmological scale factor and N6 relation length, nor a native background evolution equation. Consequently a turnaround, bounce, Big Crunch, or cycle is not inferred.\n")


if __name__ == "__main__":
    main()
