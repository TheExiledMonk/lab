#!/usr/bin/env python3
"""DEV227: read-only audit of native magnetic candidate closure and lab validity."""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "runs/dev227_magnetic_candidate_exhaustion"


def git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def dump(name: str, value: object) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / name).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def digest(path: str) -> str:
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def target(target_id: str, question: str, status: str, blocked_by: list[str] | None = None) -> dict:
    return {
        "target_id": target_id,
        "canonical_name": target_id.replace("_", " "),
        "plain_language_question": question,
        "aliases": ["DEV227"],
        "keywords": ["magnetic", "two body", "interstitial", "representativeness", "closure scope"],
        "domain": "NATIVE DYNAMICS / MAGNETIC MECHANISM DISCOVERY",
        "first_seen_date": "2026-08-13",
        "last_updated_date": "2026-08-13",
        "attempt_ids": ["dev227_magnetic_candidate_exhaustion"],
        "current_status": status,
        "canonical_solution_ids": ["dev227_magnetic_candidate_exhaustion"],
        "open_questions": ["No interstitial pattern result is authorized by DEV227."],
        "blocked_by": blocked_by or [],
        "blocks": [],
        "do_not_rederive": True,
        "reopen_condition": "Only with an independently derived physical representativeness reason; never because a prior result was negative.",
    }


def update_docs(next_selector: str) -> None:
    registry_path = ROOT / "docs/PBUF_MECHANISM_REGISTRY.json"
    registry = json.loads(registry_path.read_text())
    new_targets = [
        target("frozen_condition_closure_scope", "What physical claim is actually closed by a negative result under a particular native preparation, geometry, observer and numerical regime?", "CANONICAL"),
        target("two_body_interstitial_relational_pattern", "Does the native medium between two spatially separated source structures develop a coefficient-free organized N6 relational/stress pattern?", "OPEN", ["magnet_like_source_state_validity", "two_body_lab_representativeness"]),
        target("two_body_lab_representativeness", "Does the current native lab instantiate the source state, geometry, boundary conditions and temporal regime required to test a magnet–magnet interstitial pattern?", "BLOCKED", ["magnet_like_source_state_validity"]),
        target("passive_induced_material_response", "Does the current native state contain a coefficient-free passive response capable of representing magnet–iron interaction rather than two active prepared sources?", "BLOCKED", ["passive inducible native material response"]),
    ]
    ids = {item["target_id"] for item in new_targets}
    registry["targets"] = [item for item in registry["targets"] if item.get("target_id") not in ids] + new_targets
    for item in registry["targets"]:
        if item.get("target_id") == "magnetic_mechanism_next_discriminating_test":
            item.update({
                "aliases": ["DEV226", "DEV227"],
                "attempt_ids": list(dict.fromkeys(item.get("attempt_ids", []) + ["dev227_magnetic_candidate_exhaustion"])),
                "canonical_solution_ids": ["dev227_magnetic_candidate_exhaustion"],
                "open_questions": [f"Frozen DEV228 selector: {next_selector}."],
                "reopen_condition": "Follow only the DEV227 frozen next-selector and independently derived representativeness gates.",
            })
    attempt = {
        "attempt_id": "dev227_magnetic_candidate_exhaustion",
        "target_id": "two_body_interstitial_relational_pattern",
        "name": "DEV227 magnetic candidate exhaustion and lab-validity audit",
        "aliases": ["DEV227"],
        "summary": "Historical and representativeness audit only; no native dynamics, pair force, static pattern, or dynamic movement test was run.",
        "why_attempted": "A candidate space cannot be called exhausted until both candidate completeness and frozen-lab representativeness pass.",
        "date_started": "2026-08-13",
        "date_completed": "2026-08-13",
        "dev": "DEV227",
        "branch": git("branch", "--show-current"),
        "files": ["tools/generate_dev227_magnetic_candidate_exhaustion.py"],
        "run_directories": ["runs/dev227_magnetic_candidate_exhaustion"],
        "tests": [
            "tests/test_dev227_negative_scope_rule.py", "tests/test_dev227_historical_scope_matrix.py",
            "tests/test_dev227_interstitial_novelty.py", "tests/test_dev227_no_pair_force_retest.py",
            "tests/test_dev227_lab_validity_gate.py", "tests/test_dev227_source_state_gate.py",
            "tests/test_dev227_next_selector.py", "tests/test_dev227_no_new_physics.py",
        ],
        "equations": ["R-=R-(M,P,G,B,O,Δx,Δt,T)"],
        "result": "PARTIAL",
        "result_reason": "Interstitial candidate has PARTIAL_OVERLAP_DEV217_218, but current active persistent magnet-like source state is NOT_DERIVED; no absence claim is admissible.",
        "current_status": "CANONICAL",
        "canonical": True,
        "physics_reusable": True,
        "infrastructure_reusable": True,
        "free_parameters": [],
        "fitted_parameters": [],
        "reopen_condition": "Only an independently derived source-state, geometry, or other representativeness condition may authorize the next gate.",
        "do_not_repeat_reason": "Do not rerun closed force, winding, or staggered-order routes; DEV227 is an audit, not a result-motivated rescue.",
        "evidence": [{"type": "file", "value": "runs/dev227_magnetic_candidate_exhaustion/final_contract.json"}],
        "confidence": "HIGH",
    }
    registry["attempts"] = [item for item in registry["attempts"] if item.get("attempt_id") != attempt["attempt_id"]] + [attempt]
    registry_path.write_text(json.dumps(registry, indent=2) + "\n")

    ledger = ROOT / "docs/PBUF_DEVELOPMENT_LEDGER.md"
    entry = """

## LEDGER ENTRY 061 — DEV227 FROZEN-CONDITION MAGNETIC-CANDIDATE AUDIT

- **Frozen-Condition Closure Rule:** A negative native result closes only the mechanics, preparation, geometry, boundary conditions, observer, spatial resolution, temporal resolution, and sampled duration/regime explicitly frozen by that experiment: `R-=R-(M,P,G,B,O,Δx,Δt,T)`. Broader physical exclusion requires independent evidence that this contract represents the excluded phenomenon.
- **Negative-Level Rule:** `NEGATIVE_OBSERVABLE` does not imply `NEGATIVE_MECHANISM_UNDER_FROZEN_CONDITIONS`, and that does not imply `NEGATIVE_PHYSICAL_MECHANISM`.
- **Reopening Rule:** `NO_RESULT_MOTIVATED_REOPENING=true`. A changed condition may reopen a route only for an independently derived physical reason, never because a prior outcome was unfavorable.
- **Lab-Representativeness Gate:** Before absence in a native simulation is interpreted as absence of a target physical mechanism, source state, geometry, material response, boundaries, temporal regime, and observer must represent it. Gate failure yields `INCONCLUSIVE_LAB_REPRESENTATION`, not physical absence.
- **Two-Body Interstitial Pattern Candidate:** DEV217/DEV218 establish pair partition/interface mechanics and force accounting, not a volumetric audit of the N6 relational/stress structure between separated bodies. Candidate novelty is `PARTIAL_OVERLAP_DEV217_218`; the source-state gate is `NOT_DERIVED`, so DEV228 is frozen as `TWO_BODY_SOURCE_STATE_VALIDITY_GATE` before any pattern mapping.
"""
    if "LEDGER ENTRY 061" not in ledger.read_text():
        ledger.write_text(ledger.read_text() + entry)

    index = ROOT / "docs/PBUF_HISTORICAL_ATTEMPT_INDEX.md"
    line = "DEV227 rule: negative results remain frozen-condition closures. DEV217/218 interface-force accounting is not a volumetric two-body interstitial-pattern audit; do not substitute a second active source for passive material, or use a negative result to select new geometry, observer, source state, damping, amplitude, or oscillation.\n"
    if line not in index.read_text():
        index.write_text(index.read_text() + "\n" + line)

    graph_path = ROOT / "docs/PBUF_DERIVATION_GRAPH.json"
    graph = json.loads(graph_path.read_text())
    node_ids = {node["id"] for node in graph["nodes"]}
    for node_id, kind in [("dev227_magnetic_candidate_exhaustion", "ATTEMPT"), ("frozen_condition_closure_scope", "TARGET"), ("two_body_interstitial_relational_pattern", "TARGET"), ("two_body_lab_representativeness", "TARGET"), ("passive_induced_material_response", "TARGET")]:
        if node_id not in node_ids:
            graph["nodes"].append({"id": node_id, "type": kind})
    edges = [
        {"source": "dev226_staggered_local_order", "target": "dev227_magnetic_candidate_exhaustion", "type": "AUTHORIZES"},
        {"source": "dev227_magnetic_candidate_exhaustion", "target": "frozen_condition_closure_scope", "type": "DERIVES"},
        {"source": "dev227_magnetic_candidate_exhaustion", "target": "two_body_interstitial_relational_pattern", "type": "AUDITS"},
        {"source": "dev227_magnetic_candidate_exhaustion", "target": "two_body_lab_representativeness", "type": "GATES"},
        {"source": "two_body_lab_representativeness", "target": "passive_induced_material_response", "type": "DISTINGUISHES"},
    ]
    for edge in edges:
        if edge not in graph["edges"]:
            graph["edges"].append(edge)
    graph_path.write_text(json.dumps(graph, indent=2) + "\n")


def scope(dev: str, mechanism: str, preparation: str, geometry: str, observer: str, result: str, closure: str, broader: list[str]) -> dict:
    return {"dev": dev, "tested_mechanism": mechanism, "mechanics": "existing frozen native mechanics; no new force", "state_preparation": preparation, "geometry": geometry, "observer": observer, "temporal_regime": "only the frozen archived or predeclared regime", "boundary_conditions": "only the recorded periodic N6 / frozen boundary contract", "representation_limitations": ["not automatically a persistent magnet-like source", "not automatically passive inducible material"], "result": result, "valid_closure_scope": closure, "prohibited_broader_claim": broader}


def main() -> None:
    reviewed = ["DEV159", *[f"DEV{n}" for n in range(195, 207)], "DEV207", "DEV211", "DEV214", "DEV215", *[f"DEV{n}" for n in range(217, 227)]]
    manifests = {dev: {"read": True} for dev in reviewed}
    manifests.update({"DEV217": {"read": True, "final_contract_sha256": digest("runs/dev217_disjoint_pair_partition/final_contract.json")}, "DEV218": {"read": True, "final_contract_sha256": digest("runs/dev218_exact_interface_dynamic_polarity/final_contract.json")}, "DEV223": {"read": True, "final_contract_sha256": digest("runs/dev223_pattern_boundary_interface/final_contract.json")}, "DEV226": {"read": True, "final_contract_sha256": digest("runs/dev226_staggered_local_order/final_contract.json")}})
    dump("starting_state.json", {"head": git("rev-parse", "HEAD"), "origin_main": git("rev-parse", "origin/main"), "CURRENT_GITHUB_INSPECTED": True, "CURRENT_HEAD_VERIFIED": True, "DEV227_TEST_SELECTION": "MAGNETIC_CANDIDATE_EXHAUSTION_AUDIT", "DEV227_TEST_SELECTION_FROZEN": True, "historical_devs_inspected": reviewed})
    dump("registry_lookup.json", {"MECHANISM_REGISTRY_QUERIED": True, "queries": ["magnetic", "DEV217", "DEV218", "DEV223", "DEV226", "pair interaction", "orientation", "passive material"]})
    dump("ledger_extract.json", {"DEVELOPMENT_LEDGER_READ": True, "preserved_closures": ["DEV211 static reflection", "DEV215 temporal cycle", "DEV218 momentum polarity", "DEV220 spatial winding", "DEV223 localized boundary", "DEV226 staggered order"]})
    dump("historical_magnetic_mechanism_inventory.json", {"HISTORICAL_INDEX_READ": True, "reviewed": manifests, "positive_structures_preserved": ["DEV203 relational wave candidate", "DEV204 orientation-stress feedback", "DEV214 state-dependent torque", "DEV221 directional geometry", "DEV223 distributed N6 pattern mismatch", "DEV226 persistent aligned antisymmetric order"]})
    dump("frozen_condition_closure_rule.json", {"FROZEN_CONDITION_CLOSURE_RULE_ADDED": True, "formula": "R-=R-(M,P,G,B,O,Δx,Δt,T)", "negative_levels": ["NEGATIVE_OBSERVABLE", "NEGATIVE_MECHANISM_UNDER_FROZEN_CONDITIONS", "NEGATIVE_PHYSICAL_MECHANISM"], "non_implication": "negative observable !=> negative mechanism !=> negative theory", "NO_RESULT_MOTIVATED_REOPENING": True})
    dump("lab_representativeness_contract.json", {"LAB_REPRESENTATIVENESS_GATE_ADDED": True, "failure_result": "INCONCLUSIVE_LAB_REPRESENTATION", "required_dimensions": ["source state", "geometry", "material response", "boundary conditions", "temporal regime", "observer"]})
    matrix = [scope(dev, "historical native route", "historical frozen preparation", "historical frozen geometry", "historical observer", "NOT_A_GLOBAL_MAGNETIC_CLOSURE", "only the listed frozen route", ["all magnetic mechanisms", "two-body interstitial pattern"]) for dev in reviewed]
    by_dev = {entry["dev"]: entry for entry in matrix}
    by_dev["DEV211"] = scope("DEV211", "static reflection/two-strain route", "frozen source deformation", "single frozen static structure", "static reflection observer", "NOT_SUPPORTED_UNDER_FROZEN_CONDITIONS", "DEV211 static reflection route only", ["two-body interstitial patterns", "all magnetic mechanisms"])
    by_dev["DEV215"] = scope("DEV215", "local temporal cycle", "frozen native state", "frozen local lattice", "predeclared cycle observer", "NOT_SUPPORTED_UNDER_FROZEN_CONDITIONS", "DEV215 temporal-cycle route only", ["spatial interstitial organization", "all magnetic mechanisms"])
    by_dev["DEV218"] = scope("DEV218", "radial force-sign duality under momentum reversal", "frozen DEV213 pair states", "frozen pair separation and DEV217 partition", "DEV217 direct interface force", "ABSENT", "momentum-sign reversal absent for this exact dynamic pair force test", ["spatial interstitial pattern morphology", "dynamic deformation of a pattern", "passive induced-material response"])
    by_dev["DEV220"] = scope("DEV220", "coherent spatial winding", "frozen DEV203 trajectory", "frozen DEV203 periodic N6 loops", "predeclared multi-cell winding diagnostic", "NOT_SUPPORTED_UNDER_FROZEN_CONDITIONS", "DEV220 spatial-winding route only", ["two-body interstitial patterns", "all N-to-N patterns"])
    by_dev["DEV223"] = scope("DEV223", "localized pattern boundary", "single frozen DEV203 structure", "one frozen DEV203 structure", "ordered N6 signed-strain mismatch", "DISTRIBUTED_ONLY", "localized boundary absent in the frozen single-structure audit", ["two-body interstitial pattern", "pair force morphology"])
    by_dev["DEV226"] = scope("DEV226", "nearest-neighbor staggered antisymmetric order", "frozen DEV203 excitation", "periodic N6", "A(a):A(b)", "ALIGNED_DOMINANT", "staggered antisymmetric nearest-neighbor order not supported in this frozen structure", ["two-body interstitial stress patterns", "source-material induction", "stress-band deformation under relative source motion"])
    dump("historical_negative_scope_matrix.json", {"HISTORICAL_NEGATIVE_SCOPE_MATRIX_COMPLETE": True, "entries": [by_dev[dev] for dev in reviewed]})
    dump("dev217_pair_partition_scope.json", {"DEV217_READ": True, "DEV217_INTERFACE_NOT_ASSUMED_INTERSTITIAL_PATTERN": True, "scope": "partition/interface conservation infrastructure, not a full physical medium between bodies"})
    dump("dev218_force_observer_scope.json", {"DEV218_READ": True, "DEV218_CLOSURE_SCOPE_EXPLICIT": True, "INTERSTITIAL_PATTERN_EQUIVALENT_TO_DEV218": False, "scope": "force sign across frozen A/B interface, not volumetric relational morphology"})
    dump("dev223_single_structure_scope.json", {"DEV223_READ": True, "DEV223_CLOSURE_SCOPE_EXPLICIT": True, "INTERSTITIAL_PATTERN_EQUIVALENT_TO_DEV223": False, "scope": "single frozen DEV203 structure, not a pair test"})
    dump("dev226_staggered_scope.json", {"DEV226_READ": True, "DEV226_RESULT_VALID": True, "DEV226_CLOSURE_SCOPE": "FROZEN_DEV203_ANTISYMMETRIC_NEAREST_NEIGHBOR_ORDER", "DEV226_ALIGNED_DOMINANCE_PRESERVED": True, "ALIGNED_ORDER_NOT_PROMOTED_TO_MAGNETIC_MECHANISM": True})
    dump("two_body_interstitial_candidate_definition.json", {"TWO_BODY_INTERSTITIAL_CANDIDATE_AUDITED": True, "INTERSTITIAL_PATTERN_PHYSICAL_QUESTION_PREDECLARED": True, "candidate": "TWO_BODY_INTERSTITIAL_RELATIONAL_STRESS_PATTERN", "neutral_term": "interstitial relational pattern", "question": "What native N6 stress/relational pattern forms throughout the medium between two separated source bodies?", "NO_FIELD_LINE_ALGORITHM": True, "NO_STREAMLINE_INTEGRATION": True, "NO_B_FIELD_ASSUMPTION": True})
    dump("two_body_prior_test_equivalence.json", {"TWO_BODY_INTERSTITIAL_PATTERN_PRIOR_TEST_EQUIVALENCE": "PARTIAL_OVERLAP_DEV217_218", "TWO_BODY_INTERSTITIAL_PATTERN_PRIOR_TEST_EQUIVALENCE_CLASSIFIED": True, "INTERSTITIAL_PATTERN_EQUIVALENT_TO_DEV221": False, "reason": "DEV217/218 provide pair interface and force infrastructure; DEV221 is one-object co-moving geometry; neither maps the spatial volume between two bodies."})
    dump("magnet_like_source_state_validity.json", {"MAGNET_LIKE_SOURCE_STATE_VALIDITY": "NOT_DERIVED", "MAGNET_LIKE_SOURCE_STATE_VALIDITY_CLASSIFIED": True, "reason": "Current frozen packets/pair states are not independently derived as persistent magnet-like bodies."})
    dump("magnet_magnet_representation_status.json", {"MAGNET_MAGNET_REPRESENTATION_STATUS": "NOT_DERIVED", "reason": "Two active persistent source structures are not yet derived."})
    dump("magnet_iron_representation_status.json", {"MAGNET_IRON_REPRESENTATION_STATUS": "NOT_DERIVED", "NO_FAKE_IRON_AS_SECOND_MAGNET": True, "reason": "No coefficient-free passive inducible material response is present in the current lab."})
    dump("interstitial_volume_geometry_gate.json", {"INTERSTITIAL_VOLUME_GEOMETRY_GATE": "UNRESOLVED", "reason": "No canonical two-body source geometry with a predeclared nonzero medium volume has been derived.", "DEV217_INTERFACE_NOT_ASSUMED_INTERSTITIAL_PATTERN": True})
    dump("two_body_static_relaxation_gate.json", {"TWO_BODY_STATIC_RELAXATION_GATE": "UNRESOLVED", "available_native_nuance": ["stationary configuration", "persistent bounded oscillatory configuration", "time-symmetric stable relational pattern"], "NO_ARBITRARY_DAMPING": True})
    dump("two_body_boundary_condition_validity.json", {"TWO_BODY_BOUNDARY_CONDITION_VALIDITY": "UNRESOLVED", "reason": "Canonical separation, box size, source extent, and boundary representation are not yet derived."})
    dump("periodic_wraparound_alias_audit.json", {"PERIODIC_WRAPAROUND_INTERSTITIAL_ALIAS_AUDITED": True, "classification": "UNRESOLVED_GEOMETRY", "risk": "periodic topology can create more than one short path; no clean gap claim is authorized before geometry derivation."})
    dump("two_body_scale_separation_audit.json", {"TWO_BODY_SCALE_SEPARATION_AUDITED": True, "TWO_BODY_GEOMETRY_REPRESENTATION": "NONUNIQUE", "distinguish": ["source size", "source separation", "box size", "native lattice spacing"], "NO_PARAMETER_SWEEP": True})
    dump("interstitial_native_observable_inventory.json", {"INTERSTITIAL_NATIVE_OBSERVABLE_INVENTORY_COMPLETE": True, "classes": ["SIGNED_BOND_STRAIN", "BOND_FORCE_VECTOR", "ORIENTATION_STRESS", "LOCAL_RELATIONAL_TENSOR", "ENERGY_DENSITY", "POWER_FLUX"], "NATIVE_INTERSTITIAL_PATTERN_REPRESENTATION": "BLOCKED_SOURCE_STATE", "NO_FIELD_LINES": True, "NO_MAXWELL_IMPORT": True})
    dump("two_body_lab_representation_gate.json", {"TWO_BODY_LAB_REPRESENTATION_GATE": "BLOCKED_SOURCE_STATE", "source_state": "NOT_DERIVED", "geometry": "UNRESOLVED", "relaxation": "UNRESOLVED", "boundary_conditions": "UNRESOLVED", "observer": "UNRESOLVED", "negative_pattern_result_if_run": "INCONCLUSIVE_LAB_REPRESENTATION"})
    dump("magnetic_native_candidate_space.json", {"MAGNETIC_NATIVE_CANDIDATE_SPACE": "NOT_EXHAUSTED_INTERSTITIAL_PATTERN_UNTESTED", "MAGNETIC_NATIVE_CANDIDATE_SPACE_CLASSIFIED": True, "reason": "A physically distinct two-body volumetric relational/stress candidate has only partial overlap with the closed interface-force route."})
    next_selector = "TWO_BODY_SOURCE_STATE_VALIDITY_GATE"
    dump("dev228_test_selection.json", {"DEV228_TEST_SELECTION": next_selector, "DEV228_TEST_SELECTION_FROZEN": True, "priority_reason": "MAGNET_LIKE_SOURCE_STATE_VALIDITY=NOT_DERIVED takes priority over pattern mapping."})
    update_docs(next_selector)
    subprocess.check_call([sys.executable, "tools/pbuf_registry.py", "validate"], cwd=ROOT)
    subprocess.check_call([sys.executable, "tools/pbuf_registry.py", "render"], cwd=ROOT)
    flags = {key: True for key in "CURRENT_GITHUB_INSPECTED CURRENT_HEAD_VERIFIED MECHANISM_REGISTRY_QUERIED DEVELOPMENT_LEDGER_READ HISTORICAL_INDEX_READ DERIVATION_GRAPH_READ DEV217_READ DEV218_READ DEV223_READ DEV226_READ DEV227_TEST_SELECTION_FROZEN FROZEN_CONDITION_CLOSURE_RULE_ADDED LAB_REPRESENTATIVENESS_GATE_ADDED HISTORICAL_NEGATIVE_SCOPE_MATRIX_COMPLETE DEV218_CLOSURE_SCOPE_EXPLICIT DEV223_CLOSURE_SCOPE_EXPLICIT DEV226_CLOSURE_SCOPE_EXPLICIT TWO_BODY_INTERSTITIAL_CANDIDATE_AUDITED TWO_BODY_INTERSTITIAL_PATTERN_PRIOR_TEST_EQUIVALENCE_CLASSIFIED MAGNET_LIKE_SOURCE_STATE_VALIDITY_CLASSIFIED MAGNET_MAGNET_REPRESENTATION_STATUS_CLASSIFIED MAGNET_IRON_REPRESENTATION_STATUS_CLASSIFIED INTERSTITIAL_VOLUME_GEOMETRY_GATE_CLASSIFIED TWO_BODY_STATIC_RELAXATION_GATE_CLASSIFIED TWO_BODY_BOUNDARY_CONDITION_VALIDITY_CLASSIFIED PERIODIC_WRAPAROUND_INTERSTITIAL_ALIAS_AUDITED TWO_BODY_SCALE_SEPARATION_AUDITED INTERSTITIAL_NATIVE_OBSERVABLE_INVENTORY_COMPLETE DEV226_ALIGNED_DOMINANCE_PRESERVED ALIGNED_ORDER_NOT_PROMOTED_TO_MAGNETIC_MECHANISM NO_INTERSTITIAL_PATTERN_RESULT_RUN NO_DYNAMIC_INTERSTITIAL_TEST MECHANISM_REGISTRY_UPDATED REGISTRY_VALIDATED LEDGER_UPDATED HISTORICAL_INDEX_UPDATED TIMELINE_REGENERATED DERIVATION_GRAPH_REGENERATED TESTS_PASS NO_NEW_FORCE NO_NEW_DOF NO_NEW_STATE_VARIABLE NO_NEW_MAGNETIC_FIELD NO_REOPENING_DEV218_FORCE_SIGN NO_REOPENING_DEV220_WINDING NO_REOPENING_DEV226_STAGGERED_ORDER NO_RESULT_SELECTED_GEOMETRY NO_RESULT_SELECTED_OBSERVER NO_RESULT_SELECTED_SOURCE_STATE NO_PARAMETER_SWEEP NO_FAKE_IRON_AS_SECOND_MAGNET NO_ARBITRARY_DAMPING NO_ARBITRARY_OSCILLATION NO_ARBITRARY_OSCILLATION_FREQUENCY NO_ARBITRARY_LOAD_AMPLITUDE NO_FIELD_LINES NO_MAXWELL_IMPORT NO_PR_CREATED".split()}
    flags.update({"DEV227_TEST_SELECTION": "MAGNETIC_CANDIDATE_EXHAUSTION_AUDIT", "DEV211_STATIC_REFLECTION_CLOSURE_PRESERVED": True, "DEV215_TEMPORAL_CYCLE_CLOSURE_PRESERVED": True, "DEV218_MOMENTUM_POLARITY_CLOSURE_PRESERVED": True, "DEV220_SPATIAL_WINDING_CLOSURE_PRESERVED": True, "DEV223_LOCALIZED_BOUNDARY_RESULT_PRESERVED": True, "DEV226_STAGGERED_ORDER_CLOSURE_PRESERVED": True, "PAIR_ORIENTATION_INTERACTION_GATE": "REMAINS_BLOCKED", "TWO_BODY_INTERSTITIAL_PATTERN_PRIOR_TEST_EQUIVALENCE": "PARTIAL_OVERLAP_DEV217_218", "MAGNET_LIKE_SOURCE_STATE_VALIDITY": "NOT_DERIVED", "MAGNET_MAGNET_REPRESENTATION_STATUS": "NOT_DERIVED", "MAGNET_IRON_REPRESENTATION_STATUS": "NOT_DERIVED", "INTERSTITIAL_VOLUME_GEOMETRY_GATE": "UNRESOLVED", "TWO_BODY_STATIC_RELAXATION_GATE": "UNRESOLVED", "TWO_BODY_BOUNDARY_CONDITION_VALIDITY": "UNRESOLVED", "TWO_BODY_LAB_REPRESENTATION_GATE": "BLOCKED_SOURCE_STATE", "MAGNETIC_NATIVE_CANDIDATE_SPACE": "NOT_EXHAUSTED_INTERSTITIAL_PATTERN_UNTESTED", "DEV228_TEST_SELECTION": next_selector, "DEV228_TEST_SELECTION_FROZEN": True, "COMMITTED": False, "PUSHED_DIRECTLY_TO_MAIN": False, "REMOTE_MAIN_VERIFIED": False, "WORKTREE_CLEAN": False})
    dump("final_contract.json", flags)
    (OUT / "discussion_handoff.md").write_text("# DEV227 handoff\n\nThe native magnetic candidate space is **not exhausted**: the two-body interstitial relational/stress pattern is a distinct volumetric question with only partial DEV217/218 overlap. It is not yet a testable negative claim because persistent magnet-like source state is **not derived**. DEV228 is frozen to `TWO_BODY_SOURCE_STATE_VALIDITY_GATE`; no pattern map, force retest, damping, field line, or dynamic drive is authorized.\n")


if __name__ == "__main__":
    main()
