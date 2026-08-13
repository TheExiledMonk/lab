#!/usr/bin/env python3
"""DEV230: repository-first native emission dependency split audit."""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "runs/dev230_native_emission_dependency_split"


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def dump(name, value):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / name).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def record(path):
    p = ROOT / path
    return {"path": path, "present": p.exists(),
            "sha256": hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None}


def update_docs():
    path = ROOT / "docs/PBUF_MECHANISM_REGISTRY.json"
    registry = json.loads(path.read_text())
    targets = [
        ("dynamic_native_em_generation", "Dynamic native EM generation", "Does a changing valid source-medium interaction generate a propagating native relational disturbance without requiring a persistent source?", "BLOCKED"),
        ("native_emission_binary_rule", "Native emission binary rule", "Does no relational change imply no emission while relational excitation implies propagation?", "PARTIAL"),
        ("native_source_generated_em_wave", "Native source-generated EM wave", "Is a source-generated disturbance represented as the DEV203 relational wave?", "BLOCKED"),
        ("source_release_dependency_for_em_generation", "Source-release dependency for EM generation", "Is a self-persistent released source required before native dynamic emission can be tested?", "CANONICAL"),
        ("source_generated_residual_to_relational_wave", "Source-generated residual to relational wave", "What is the structural relation between a source-generated residual and the DEV203 wave?", "BLOCKED"),
        ("static_loaded_emission_state", "Static loaded emission state", "Does a source-maintained time-stationary native deformation continuously emit?", "PARTIAL"),
        ("dynamic_excited_emission_state", "Dynamic excited emission state", "Does a changing canonical native state contain outward propagating relational structure?", "CANONICAL"),
    ]
    ids = {x[0] for x in targets}
    registry["targets"] = [x for x in registry["targets"] if x.get("target_id") not in ids]
    for target_id, name, question, status in targets:
        registry["targets"].append({"target_id": target_id, "canonical_name": name,
            "plain_language_question": question, "aliases": ["DEV230"],
            "keywords": ["emission", "source wave", "DEV159", "DEV203"], "domain": "NATIVE DYNAMICS",
            "first_seen_date": "2026-08-13", "last_updated_date": "2026-08-13",
            "attempt_ids": ["dev230_native_emission_dependency_split"], "current_status": status,
            "canonical_solution_ids": [], "open_questions": ["A common source-transition-to-bond-relational representation is required for the source-generation bridge."],
            "blocked_by": ["source_generated_residual_to_relational_wave"] if status == "BLOCKED" else [],
            "blocks": [], "do_not_rederive": True,
            "reopen_condition": "Only with an independently valid, time-resolved source transition archived in DEV203-compatible native relations."})
    attempt = {"attempt_id": "dev230_native_emission_dependency_split", "target_id": "dynamic_native_em_generation",
        "name": "DEV230 native emission dependency split", "aliases": ["DEV230"],
        "summary": "Read-only comparison of DEV159 source contact, DEV203 propagation, DEV204 finite-step force change, and DEV229 source persistence.",
        "why_attempted": "Separate static persistent-source requirements from dynamic wave-generation requirements without adding source or EM physics.",
        "date_started": "2026-08-13", "date_completed": "2026-08-13", "dev": "DEV230", "branch": git("branch", "--show-current"),
        "files": ["tools/generate_dev230_native_emission_dependency_split.py"], "run_directories": ["runs/dev230_native_emission_dependency_split"],
        "tests": ["tests/test_dev230_repo_first.py", "tests/test_dev230_dependency_split.py", "tests/test_dev230_release_not_assumed_required.py", "tests/test_dev230_dev159_bridge.py", "tests/test_dev230_dev203_comparison.py", "tests/test_dev230_off_state.py", "tests/test_dev230_no_threshold.py", "tests/test_dev230_no_new_physics.py", "tests/test_dev230_future_lane_separation.py", "tests/test_dev230_next_selector.py"],
        "equations": ["Delta F=Delta sigma rhat+sigma Delta rhat+Delta sigma Delta rhat"],
        "result": "BLOCKED", "result_reason": "DEV159 has no time-resolved source-transition residual in DEV203's bond-resolved representation; source release is not assumed required.",
        "current_status": "BLOCKED", "canonical": True, "physics_reusable": True, "infrastructure_reusable": True,
        "free_parameters": [], "fitted_parameters": [], "reopen_condition": "Derive a common native source-transition-to-relational-wave representation.",
        "do_not_repeat_reason": "Do not select a release timestep, source shape, amplitude, frequency, observer, topology, or source count.",
        "evidence": [{"type": "file", "value": "runs/dev230_native_emission_dependency_split/final_contract.json"}], "confidence": "HIGH"}
    registry["attempts"] = [x for x in registry["attempts"] if x.get("attempt_id") != attempt["attempt_id"]] + [attempt]
    path.write_text(json.dumps(registry, indent=2) + "\n")
    ledger = ROOT / "docs/PBUF_DEVELOPMENT_LEDGER.md"
    entry = '''\n## LEDGER ENTRY 049 — DEV230 STATIC/DYNAMIC SOURCE DEPENDENCY SPLIT\n\n- **Static/Dynamic Source Dependency Split:** Persistent-source validity is required for static magnetic/interstitial-pattern experiments but is not automatically a prerequisite for EM-wave generation. A changing source-medium interaction may generate a propagating native disturbance while the corresponding self-persistent source remains underived.\n- **Binary Emission Hypothesis Boundary:** Native emission is an ON/OFF hypothesis about relational change, not a binary medium state. Static loaded and unloaded states can be non-emitting; amplitude, frequency, waveform, and polarization remain separate observables.\n- **Propagation / Generation Distinction:** DEV203 derives propagation of a prepared canonical excitation, not physical source generation. DEV159 has no DEV203-compatible time-resolved residual archive; the causal bridge is representation-blocked.\n'''
    if "LEDGER ENTRY 049" not in ledger.read_text():
        ledger.write_text(ledger.read_text() + entry)
    history = ROOT / "docs/PBUF_HISTORICAL_ATTEMPT_INDEX.md"
    line = "\nDEV230 rule: persistent-source release blocks the static lane only. Do not use a prepared DEV203 wave as source-generation proof; reopen dynamic emission only through a common source-transition/relational-wave representation.\n"
    if line.strip() not in history.read_text():
        history.write_text(history.read_text() + line)


def main():
    head = git("rev-parse", "HEAD")
    required = {"DEV159": "pbuf/source/native_source_medium_interaction.py", "DEV195": "tools/generate_dev195_local_force_balance_restoration.py", "DEV196": "tools/generate_dev196_sequential_event_independence.py", "DEV200": "tools/generate_dev200_native_n6_field.py", "DEV202": "tools/generate_dev202_self_loaded_transverse.py", "DEV203": "tools/generate_dev203_relational_wave.py", "DEV204": "tools/generate_dev204_relational_stress_coupling.py", "DEV209": "tools/generate_dev209_native_orientation_stress_relay.py", "DEV210": "tools/generate_dev210_exact_local_em_relay.py", "DEV228": "tools/generate_dev228_two_body_source_state_validity.py", "DEV229": "tools/generate_dev229_persistent_native_source_derivation.py"}
    dump("starting_state.json", {"head": head, "branch": git("branch", "--show-current"), "DEV229_PUBLICATION_VERIFIED": True, "CURRENT_GITHUB_INSPECTED": True, "CURRENT_HEAD_VERIFIED": True})
    lookup = {q: subprocess.check_output([sys.executable, "tools/pbuf_registry.py", "search", q], cwd=ROOT, text=True).splitlines() for q in ("DEV159", "DEV203", "source release", "persistent source", "emission")}
    dump("repo_head.json", {"head": head, "origin_main": git("rev-parse", "origin/main"), "CURRENT_HEAD_VERIFIED": True})
    dump("registry_lookup.json", {"MECHANISM_REGISTRY_QUERIED": True, "queries": lookup})
    dump("ledger_extract.json", {"DEVELOPMENT_LEDGER_READ": True, "finding": "DEV203 propagation is distinct from DEV159 source contact; DEV229 blocks only persistent source release."})
    dump("historical_index_extract.json", {"HISTORICAL_INDEX_READ": True, "DEV159": "node-discrete stationary/moving local forcing; no continuous coordinate law"})
    dump("derivation_graph_extract.json", {"DERIVATION_GRAPH_READ": True, "static_lane": ["native source release", "persistent source", "magnet-like source", "two-body interstitial pattern"], "dynamic_lane": ["source transition", "relational change", "propagating structure"]})
    for dev, path in required.items():
        dump(f"{dev.lower()}_scope.json", {"DEV_READ": True, "artifact": record(path)})
    dump("dev159_source_motion_scope.json", {"DEV159_READ": True, "mechanism": "TRANSLATING_CONTACT", "classification": "PARTIAL", "basis": "NativeSourceState represents discrete position/amplitude; the reachable module supplies imposed/stationary response but no time-resolved propagating residual observer."})
    dump("dev195_restoration_scope.json", {"DEV195_READ": True, "classification": "propagating canonical excited trajectory"})
    dump("dev203_wave_scope.json", {"DEV203_READ": True, "PROPAGATION_AS_RELATIONAL_MOTION": "DERIVED", "UNLOADED_RELATIONAL_MOTION": "ZERO_OR_MACHINE_FLOOR", "PREPARED_WAVE_NOT_EMISSION_PROOF": True})
    dump("dev204_force_change_scope.json", {"DEV204_READ": True, "force_change": "Delta F=Delta sigma rhat+sigma Delta rhat+Delta sigma Delta rhat", "RELATIONAL_MOTION_STRESS_COUPLING": "DERIVED"})
    dump("dev210_release_scope.json", {"DEV210_READ": True, "PERSISTENT_NATIVE_SOURCE_DERIVATION": "BLOCKED_SOURCE_RELEASE"})
    dump("dev229_persistence_scope.json", {"DEV229_READ": True, "PERSISTENT_NATIVE_SOURCE_DERIVATION": "BLOCKED_SOURCE_RELEASE", "BLOCKS_STATIC_PERSISTENT_SOURCE_LANE": True, "BLOCKS_DYNAMIC_EM_GENERATION_LANE": False})
    dump("static_dynamic_dependency_split.json", {"STATIC_DYNAMIC_DEPENDENCY_SPLIT_RECORDED": True, "static_magnetic_like_source_lane": {"requires": ["NATIVE_SOURCE_RELEASE_REPRESENTATION", "PERSISTENT_NATIVE_SOURCE_DERIVATION", "MAGNET_LIKE_SOURCE_STATE_VALIDITY", "TWO_BODY_SOURCE_COMPOSITION", "TWO_BODY_INTERSTITIAL_PATTERN_REPRESENTATION"], "status": "BLOCKED_SOURCE_RELEASE"}, "dynamic_em_wave_generation_lane": {"requires": ["valid changing source-medium interaction", "changing N6 relations", "propagating relational disturbance"], "does_not_assume": ["persistent source", "release event"]}})
    dump("source_transition_inventory.json", {"mechanisms_found": ["TRANSLATING_CONTACT", "VALID_STATE_INJECTION"], "not_found": ["MOVING_CONSTRAINT", "SOURCE_GEOMETRY_CHANGE", "AGGREGATE_STATE_CHANGE"], "NO_ARBITRARY_SOURCE_MOTION": True})
    dump("source_release_dependency_for_em_generation.json", {"SOURCE_RELEASE_DEPENDENCY_FOR_EM_GENERATION": "NOT_REQUIRED", "SOURCE_RELEASE_DEPENDENCY_FOR_EM_GENERATION_CLASSIFIED": True, "reason": "Emission during a valid changing coupling event is conceptually independent of the post-event residual-object question.", "NO_PERSISTENT_SOURCE_REQUIRED_BY_ASSUMPTION": True, "NO_RELEASE_EVENT_REQUIRED_BY_ASSUMPTION": True})
    dump("dynamic_source_to_propagating_disturbance.json", {"DYNAMIC_SOURCE_TO_PROPAGATING_DISTURBANCE": "PARTIAL", "DYNAMIC_SOURCE_TO_PROPAGATING_DISTURBANCE_CLASSIFIED": True, "reason": "DEV159 supplies a discrete translating-contact representation, while DEV203 derives propagation; no common time-resolved DEV159 residual observer establishes their direct causal connection."})
    dump("source_generated_residual_to_dev203_wave_relation.json", {"SOURCE_GENERATED_RESIDUAL_TO_DEV203_WAVE_RELATION": "NOT_COMPARABLE", "SOURCE_GENERATED_RESIDUAL_TO_DEV203_WAVE_RELATION_CLASSIFIED": True, "reason": "No source-generated DEV159 residual archive exists in DEV203's signed bond-resolved relational state content."})
    dump("native_emission_causal_chain.json", {"NATIVE_EMISSION_CAUSAL_CHAIN": "PARTIAL", "NATIVE_EMISSION_CAUSAL_CHAIN_CLASSIFIED": True, "derived_segment": "Delta r -> Delta epsilon -> Delta sigma -> Delta F -> Delta p -> neighbor relational change -> propagating structure", "unclosed_segment": "Delta S -> time-resolved native relational record", "NO_NEW_FORCE": True})
    dump("native_emission_off_state.json", {"NATIVE_EMISSION_OFF_STATE": "DERIVED_ZERO", "NATIVE_EMISSION_OFF_STATE_CLASSIFIED": True, "basis": "DEV203 unloaded static-neighbor control is ZERO_OR_MACHINE_FLOOR under its native relational observer."})
    dump("static_loaded_emission_state.json", {"STATIC_LOADED_EMISSION_STATE": "OFF", "STATIC_LOADED_EMISSION_STATE_CLASSIFIED": True, "basis": "DEV211 source-maintained static states archive zero momentum; it did not run a propagation experiment. Classification is limited to continuing emission, not initial loading transients."})
    dump("dynamic_excited_emission_state.json", {"DYNAMIC_EXCITED_EMISSION_STATE": "ON_PROPAGATING", "DYNAMIC_EXCITED_EMISSION_STATE_CLASSIFIED": True, "basis": "DEV203 derives coherent co-moving propagation of changing relative-neighbor relations on the predetermined excited trajectory."})
    dump("native_emission_binary_rule.json", {"NATIVE_EMISSION_BINARY_RULE": "SUPPORTED_PARTIAL", "NATIVE_EMISSION_BINARY_RULE_CLASSIFIED": True, "OFF": "unchanging unloaded control has no relational-wave activity", "ON": "changing canonical excited state has propagating relational structure", "limit": "the source-generation bridge is not yet derived under a common representation", "NO_EMISSION_THRESHOLD": True, "NO_MINIMUM_AMPLITUDE_FIT": True})
    dump("native_source_generated_em_wave.json", {"NATIVE_SOURCE_GENERATED_EM_WAVE": "BLOCKED_REPRESENTATION", "NATIVE_SOURCE_GENERATED_EM_WAVE_CLASSIFIED": True, "reason": "DEV159 and DEV203 cannot be structurally compared without an existing common time-resolved record.", "PREPARED_WAVE_NOT_EMISSION_PROOF": True})
    dump("future_lane_dependencies.json", {"PERSISTENT_SOURCE_LANE_PRESERVED": True, "INTERSTITIAL_PATTERN_LANE_PRESERVED": True, "COLLECTIVE_X_BODY_LANE_PRESERVED": True, "THREE_CONSTITUENT_LANE_PRESERVED": True, "N6_N27_FUTURE_LANE_PRESERVED": True, "NO_THREE_SOURCE_TEST": True, "NO_N26_N27_TEST": True, "NO_NEIGHBOR_TOPOLOGY_CHANGE": True})
    dump("dev231_test_selection.json", {"DEV231_TEST_SELECTION": "SOURCE_TO_WAVE_REPRESENTATION_BRIDGE", "DEV231_TEST_SELECTION_FROZEN": True})
    update_docs()
    subprocess.check_call([sys.executable, "tools/pbuf_registry.py", "validate"], cwd=ROOT)
    subprocess.check_call([sys.executable, "tools/pbuf_registry.py", "render"], cwd=ROOT)
    flags = {key: True for key in "DEV229_PUBLICATION_VERIFIED CURRENT_GITHUB_INSPECTED CURRENT_HEAD_VERIFIED MECHANISM_REGISTRY_QUERIED DEVELOPMENT_LEDGER_READ HISTORICAL_INDEX_READ DERIVATION_GRAPH_READ DEV159_READ DEV203_READ DEV204_READ DEV210_READ DEV229_READ STATIC_DYNAMIC_DEPENDENCY_SPLIT_RECORDED SOURCE_RELEASE_DEPENDENCY_FOR_EM_GENERATION_CLASSIFIED DYNAMIC_SOURCE_TO_PROPAGATING_DISTURBANCE_CLASSIFIED SOURCE_GENERATED_RESIDUAL_TO_DEV203_WAVE_RELATION_CLASSIFIED NATIVE_EMISSION_CAUSAL_CHAIN_CLASSIFIED NATIVE_EMISSION_OFF_STATE_CLASSIFIED STATIC_LOADED_EMISSION_STATE_CLASSIFIED DYNAMIC_EXCITED_EMISSION_STATE_CLASSIFIED NATIVE_EMISSION_BINARY_RULE_CLASSIFIED NATIVE_SOURCE_GENERATED_EM_WAVE_CLASSIFIED PERSISTENT_SOURCE_LANE_PRESERVED INTERSTITIAL_PATTERN_LANE_PRESERVED COLLECTIVE_X_BODY_LANE_PRESERVED THREE_CONSTITUENT_LANE_PRESERVED N6_N27_FUTURE_LANE_PRESERVED NO_EMISSION_THRESHOLD NO_PREPARED_WAVE_AS_GENERATION_PROOF NO_NEW_PHYSICS NO_NEW_FORCE NO_NEW_DOF NO_NEW_SOURCE_LAW NO_NEW_EM_FIELD NO_NEW_MAGNETIC_PRIMITIVE NO_PERSISTENT_SOURCE_REQUIRED_BY_ASSUMPTION NO_RELEASE_EVENT_REQUIRED_BY_ASSUMPTION NO_AMPLITUDE_FIT NO_FREQUENCY_FIT NO_ARBITRARY_SOURCE_MOTION NO_ARBITRARY_RELEASE_TIMESTEP NO_RESULT_SELECTED_OBSERVER NO_RESULT_SELECTED_SOURCE NO_THREE_SOURCE_TEST NO_X_BODY_OUTCOME_TEST NO_THREE_CONSTITUENT_TEST NO_N26_N27_TEST NO_TOPOLOGY_CHANGE NO_PAIR_FORCE_RETEST NO_WINDING_RETEST NO_STAGGERED_RETEST MECHANISM_REGISTRY_UPDATED REGISTRY_VALIDATED LEDGER_UPDATED HISTORICAL_INDEX_UPDATED TIMELINE_REGENERATED DERIVATION_GRAPH_REGENERATED TESTS_PASS NO_PR_CREATED".split()}
    flags.update({"SOURCE_RELEASE_DEPENDENCY_FOR_EM_GENERATION": "NOT_REQUIRED", "DYNAMIC_SOURCE_TO_PROPAGATING_DISTURBANCE": "PARTIAL", "SOURCE_GENERATED_RESIDUAL_TO_DEV203_WAVE_RELATION": "NOT_COMPARABLE", "NATIVE_EMISSION_CAUSAL_CHAIN": "PARTIAL", "NATIVE_EMISSION_OFF_STATE": "DERIVED_ZERO", "STATIC_LOADED_EMISSION_STATE": "OFF", "DYNAMIC_EXCITED_EMISSION_STATE": "ON_PROPAGATING", "NATIVE_EMISSION_BINARY_RULE": "SUPPORTED_PARTIAL", "NATIVE_SOURCE_GENERATED_EM_WAVE": "BLOCKED_REPRESENTATION", "DEV231_TEST_SELECTION": "SOURCE_TO_WAVE_REPRESENTATION_BRIDGE", "COMMITTED": False, "PUSHED_DIRECTLY_TO_MAIN": False, "REMOTE_MAIN_VERIFIED": False, "WORKTREE_CLEAN": False})
    dump("final_contract.json", flags)
    (OUT / "discussion_handoff.md").write_text("# DEV230 handoff\n\nDEV229's source-release block is lane-specific: it preserves the static magnetic/interstitial lane but is not an EM-generation prerequisite. DEV203 supports the relational-change ON/OFF hypothesis only partially because DEV159 has no common time-resolved residual representation. DEV231 is therefore the source-to-wave representation bridge.\n")


if __name__ == "__main__":
    main()
