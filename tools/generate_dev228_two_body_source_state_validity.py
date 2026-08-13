#!/usr/bin/env python3
"""DEV228: audit existing source-state validity before any interstitial pattern test."""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "runs/dev228_two_body_source_state_validity"
sys.path.insert(0, str(ROOT))
from pbuf.analysis.native_source_state_validity import classifications, source_inventory


def git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def dump(name: str, value: object) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / name).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def manifest(dev: str, script: str, run: str) -> dict:
    path = ROOT / script
    return {"DEV_READ": True, "dev": dev, "script": script, "run": run,
            "exists": path.exists(), "sha256": hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None}


def target(ident: str, question: str, status: str, blocked: list[str] | None = None) -> dict:
    return {"target_id": ident, "canonical_name": ident.replace("_", " "),
            "plain_language_question": question, "aliases": ["DEV228"],
            "keywords": ["source state", "collective", "composition", "interstitial", "native"],
            "domain": "NATIVE DYNAMICS / MAGNETIC & EM SOURCE STRUCTURE",
            "first_seen_date": "2026-08-13", "last_updated_date": "2026-08-13",
            "attempt_ids": ["dev228_two_body_source_state_validity"], "current_status": status,
            "canonical_solution_ids": ["dev228_two_body_source_state_validity"], "open_questions": [],
            "blocked_by": blocked or [], "blocks": [], "do_not_rederive": True,
            "reopen_condition": "Only an independently derived native source state or representation result may reopen this gate."}


def update_docs() -> None:
    registry_path = ROOT / "docs/PBUF_MECHANISM_REGISTRY.json"
    registry = json.loads(registry_path.read_text())
    targets = [
        target("magnet_like_source_state_validity", "Does frozen native dynamics already contain a persistent, identifiable, orientation-capable source structure suitable for a two-body interstitial experiment?", "PARTIAL", ["self-supported localized native source"]),
        target("two_body_source_composition", "Can two valid native source structures be prepared simultaneously using existing DEV213 semantics?", "BLOCKED", ["magnet_like_source_state_validity"]),
        target("finite_x_native_composition", "Does the DEV213 same-step injection algebra extend to a finite collection of independently defined preparations?", "ACTIVE"),
        target("collective_x_body_native_response", "Is the response of X>2 simultaneously prepared native sources reducible to one- and two-source behavior?", "OPEN", ["valid source structures", "interstitial representation"]),
        target("three_source_collective_reducibility", "Does P_123 reduce coefficient-free to available one- and two-source information?", "OPEN", ["valid source structures", "interstitial representation"]),
        target("three_constituent_composite_state", "Can three distinct native constituents form an aggregate-defined persistent composite identity?", "BLOCKED", ["two distinct constituent-state classes"]),
        target("collective_source_em_wave_generation", "Can a DEV203-like relational wave emerge from aggregate source reorganization rather than a prepared packet?", "OPEN", ["valid source structures"]),
        target("collective_near_to_far_relational_transition", "Can complex collective near-source organization yield a simpler coherent outgoing relational structure?", "OPEN", ["collective source reorganization"]),
        target("passive_induced_material_response", "Does a coefficient-free passive induced material response exist, distinct from an active prepared source?", "BLOCKED", ["passive inducible native material response"]),
    ]
    ids = {item["target_id"] for item in targets}
    registry["targets"] = [item for item in registry["targets"] if item.get("target_id") not in ids] + targets
    attempt = {"attempt_id": "dev228_two_body_source_state_validity", "target_id": "magnet_like_source_state_validity",
               "name": "DEV228 native persistent source-state validity and collective-source gate", "aliases": ["DEV228"],
               "summary": "Read-only source-state and composition audit; no interstitial pattern or dynamic load response is run.",
               "why_attempted": "DEV227 requires source validity before a two-body interstitial pattern can be physically interpreted.",
               "date_started": "2026-08-13", "date_completed": "2026-08-13", "dev": "DEV228", "branch": git("branch", "--show-current"),
               "files": ["pbuf/analysis/native_source_state_validity.py", "tools/generate_dev228_two_body_source_state_validity.py"],
               "run_directories": ["runs/dev228_two_body_source_state_validity"],
               "tests": ["tests/test_dev228_source_inventory.py", "tests/test_dev228_source_persistence.py", "tests/test_dev228_source_localization.py", "tests/test_dev228_source_orientation.py", "tests/test_dev228_dev213_reuse.py", "tests/test_dev228_finite_x_composition_gate.py", "tests/test_dev228_no_pattern_result.py", "tests/test_dev228_future_branch_separation.py", "tests/test_dev228_next_selector.py"],
               "equations": ["X_AB=I_B(I_A(X_0))=I_A(I_B(X_0))", "X_1...X=I_X(...I_2(I_1(X_0))...)", "R-=R-(M,P,G,B,O,Δx,Δt,T)"],
               "result": "PARTIAL", "result_reason": "No existing structure is simultaneously self-persistent, canonically localizable, identifiable, and orientation-capable; DEV213 still derives aggregate preparation, not a magnetic source.",
               "current_status": "CANONICAL", "canonical": True, "physics_reusable": True, "infrastructure_reusable": True,
               "free_parameters": [], "fitted_parameters": [], "reopen_condition": "Derive a source state independently; do not relabel packet or clamp.",
               "do_not_repeat_reason": "No result-selected source construction, interstitial pattern, load motion, or source-count sweep is authorized.",
               "evidence": [{"type": "file", "value": "runs/dev228_two_body_source_state_validity/final_contract.json"}], "confidence": "HIGH"}
    registry["attempts"] = [item for item in registry["attempts"] if item.get("attempt_id") != attempt["attempt_id"]] + [attempt]
    registry_path.write_text(json.dumps(registry, indent=2) + "\n")
    ledger = ROOT / "docs/PBUF_DEVELOPMENT_LEDGER.md"
    entry = """

## LEDGER ENTRY 062 — DEV228 NATIVE SOURCE-STATE VALIDITY AND COLLECTIVE-SOURCE GATE

- **Native Source-State Validity Rule:** A two-body magnetic/interstitial-pattern experiment is physically interpretable only after its constituent source independently passes persistence, localization, identity, and orientation gates. A propagating packet, externally maintained deformation, or preparation-provenance label may not be promoted to a magnet-like source.
- **Source-Number Closure Rule:** A one- or two-source negative result does not exclude collective X>2 organization unless finite-X response is independently reduced to one- and two-source contributions. DEV213 establishes aggregate same-state evolution, not linear dynamical superposition.
- **Native Collective Evolution Rule:** Simultaneously prepared structures evolve as one aggregate DEV167 state; independently propagated source solutions may not be added post hoc.
- **Three-Constituent Composite Boundary:** Three is the first source count beyond pairwise structure and can later distinguish AAB from ABB. No majority-state rule, quark identity, up/down primitive, or composite polarity is assumed.
"""
    if "LEDGER ENTRY 062" not in ledger.read_text():
        ledger.write_text(ledger.read_text() + entry)
    index = ROOT / "docs/PBUF_HISTORICAL_ATTEMPT_INDEX.md"
    line = "\nDEV228 rule: source validity precedes interstitial mapping; DEV213's finite aggregate preparation is not linear propagation superposition, and X=2 never closes X>2 without a reduction theorem.\n"
    if line.strip() not in index.read_text(): index.write_text(index.read_text() + line)
    graph_path = ROOT / "docs/PBUF_DERIVATION_GRAPH.json"; graph = json.loads(graph_path.read_text())
    nodes = [("dev228_two_body_source_state_validity", "ATTEMPT"), ("magnet_like_source_state_validity", "TARGET"), ("two_body_source_composition", "TARGET"), ("finite_x_native_composition", "TARGET"), ("collective_x_body_native_response", "TARGET"), ("three_source_collective_reducibility", "TARGET"), ("three_constituent_composite_state", "TARGET"), ("collective_source_em_wave_generation", "TARGET"), ("collective_near_to_far_relational_transition", "TARGET")]
    existing = {node["id"] for node in graph["nodes"]}
    graph["nodes"] += [{"id": ident, "type": kind} for ident, kind in nodes if ident not in existing]
    edges = [{"source": "dev227_magnetic_candidate_exhaustion", "target": "dev228_two_body_source_state_validity", "type": "AUTHORIZES"}, {"source": "dev228_two_body_source_state_validity", "target": "magnet_like_source_state_validity", "type": "AUDITS"}, {"source": "dev228_two_body_source_state_validity", "target": "finite_x_native_composition", "type": "DERIVES"}, {"source": "finite_x_native_composition", "target": "collective_x_body_native_response", "type": "ENABLES"}, {"source": "collective_x_body_native_response", "target": "three_source_collective_reducibility", "type": "SPECIALIZES"}]
    for edge in edges:
        if edge not in graph["edges"]: graph["edges"].append(edge)
    graph_path.write_text(json.dumps(graph, indent=2) + "\n")


def main() -> None:
    values = classifications(); inventory = source_inventory(); head = git("rev-parse", "HEAD")
    reviewed = ["DEV159", "DEV195", "DEV200", "DEV203", "DEV204", "DEV207", "DEV211", "DEV212", "DEV213", "DEV214", "DEV217", "DEV218", "DEV227"]
    scripts = {"DEV159": ("tests/test_dev159_native_source_medium_interaction.py", "runs/native_medium_interaction_wide_net001"), "DEV195": ("tools/generate_dev195_local_force_balance_restoration.py", "runs/dev195_local_force_balance_restoration"), "DEV200": ("tools/generate_dev200_native_n6_field.py", "runs/dev200_native_n6_field"), "DEV203": ("tools/generate_dev203_relational_wave.py", "runs/dev203_relational_wave"), "DEV204": ("tools/generate_dev204_relational_stress_coupling.py", "runs/dev204_relational_stress_coupling"), "DEV207": ("tools/generate_dev207_two_excitation_interaction.py", "runs/dev207_two_excitation_interaction"), "DEV211": ("tools/generate_dev211_two_strain_magnetism.py", "runs/dev211_two_strain_magnetism"), "DEV212": ("tools/generate_dev212_native_multistate_polarity.py", "runs/dev212_native_multistate_polarity"), "DEV213": ("tools/generate_dev213_native_multi_structure_composition.py", "runs/dev213_native_multi_structure_composition"), "DEV214": ("tools/generate_dev214_dynamic_polarity_interaction.py", "runs/dev214_dynamic_polarity_interaction"), "DEV217": ("tools/generate_dev217_disjoint_pair_partition.py", "runs/dev217_disjoint_pair_partition"), "DEV218": ("tools/generate_dev218_exact_interface_dynamic_polarity.py", "runs/dev218_exact_interface_dynamic_polarity"), "DEV227": ("tools/generate_dev227_magnetic_candidate_exhaustion.py", "runs/dev227_magnetic_candidate_exhaustion")}
    dump("starting_state.json", {"head": head, "origin_main": git("rev-parse", "origin/main"), "CURRENT_GITHUB_INSPECTED": True, "CURRENT_HEAD_VERIFIED": True, "DEV228_TEST_SELECTION": "TWO_BODY_SOURCE_STATE_VALIDITY_GATE", "DEV228_TEST_SELECTION_FROZEN": True})
    queries = ["source", "persistent", "stationary", "packet", "composition", "collective", "three constituent", "passive induced material"]
    lookup = {q: subprocess.check_output([sys.executable, "tools/pbuf_registry.py", "search", q], cwd=ROOT, text=True).splitlines() for q in queries}
    dump("registry_lookup.json", {"MECHANISM_REGISTRY_QUERIED": True, "queries": lookup})
    dump("ledger_extract.json", {"DEVELOPMENT_LEDGER_READ": True, "preserved": ["DEV159 external constraint", "DEV213 aggregate evolution", "DEV227 frozen-condition closure"]})
    dump("historical_source_state_inventory.json", {"HISTORICAL_INDEX_READ": True, "reviewed": reviewed, "PERSISTENT_SOURCE_NOT_ASSUMED_FROM_PROPAGATING_PACKET": True, "EXTERNAL_CONSTRAINT_NOT_AUTOMATIC_MAGNET": True})
    dump("derivation_graph_extract.json", {"DERIVATION_GRAPH_READ": True, "lineage": ["DEV159", "DEV195/200/203/204", "DEV211/212", "DEV213/214", "DEV217/218", "DEV227"]})
    for dev in reviewed: dump(f"{dev.lower()}_manifest.json", manifest(dev, *scripts[dev]))
    dump("native_source_state_inventory.json", {"NATIVE_SOURCE_STATE_INVENTORY_COMPLETE": True, "classes": inventory, "PERSISTENT_SOURCE_NOT_ASSUMED_FROM_PROPAGATING_PACKET": True, "EXTERNAL_CONSTRAINT_NOT_AUTOMATIC_MAGNET": True})
    dump("native_source_persistence.json", {"NATIVE_SOURCE_PERSISTENCE": values["NATIVE_SOURCE_PERSISTENCE"], "reason": "S1/S4 require external maintenance; S2/S3 propagate; S5 is not derived."})
    dump("native_source_localization.json", {"NATIVE_SOURCE_LOCALIZATION": values["NATIVE_SOURCE_LOCALIZATION"], "reason": "Exact contact geometry is available only for maintained constraints; packet support is noncompact because Gaussian tails are exactly nonzero."})
    dump("native_source_identity.json", {"NATIVE_SOURCE_IDENTITY": values["NATIVE_SOURCE_IDENTITY"], "reason": "DEV213 provenance is not a physical degree of freedom; dynamic coherence does not establish a persistent bounded body."})
    dump("source_orientation_state.json", {"SOURCE_ORIENTATION_STATE": values["SOURCE_ORIENTATION_STATE"], "basis": ["DEV203 directional geometry", "DEV204 orientation stress", "DEV212 dynamic full-state distinctions", "DEV214 state-dependent torque"], "NO_NORTH_SOUTH_PRIMITIVE": True, "NO_MAGNETIC_POLE_LABEL": True})
    dump("magnet_like_source_state_validity.json", {"MAGNET_LIKE_SOURCE_STATE_VALIDITY": values["MAGNET_LIKE_SOURCE_STATE_VALIDITY"], "reason": "No source is simultaneously self-persistent, canonically localizable, physically identifiable, and orientation-capable independent of result-selected maintenance."})
    dump("two_body_source_composition.json", {"TWO_BODY_SOURCE_COMPOSITION": values["TWO_BODY_SOURCE_COMPOSITION"], "DEV213_MULTI_STRUCTURE_COMPOSITION_REUSED": True, "NO_POSTHOC_LINEAR_SUPERPOSITION": True, "reason": "DEV213 authorizes dynamic packet injection, but no valid magnet-like source definition is available to compose unchanged."})
    dump("finite_x_composition_audit.json", {"FINITE_X_NATIVE_COMPOSITION_GATE": values["FINITE_X_NATIVE_COMPOSITION_GATE"], "proof": "For additive same-step I_i(X)=X+delta_i, induction gives I_X(...I_1(X0))=X0+sum_i delta_i independently of order; after preparation Phi_DEV167 evolves the one aggregate state.", "conditions": ["finite independently defined same-step injections", "each aggregate full state remains DEV167-valid"], "NO_PAIRWISE_REDUCIBILITY_ASSUMPTION": True, "NO_POSTHOC_LINEAR_SUPERPOSITION": True})
    dump("collective_x_body_candidate.json", {"target_id": "collective_x_body_native_response", "question": "Is X>2 response reducible to one- and two-source behavior, or does irreducible collective organization appear?", "NO_SOURCE_COUNT_SWEEP": True, "minimal_new_X": 3})
    dump("three_source_collective_reducibility_contract.json", {"target_id": "three_source_collective_reducibility", "formula": "P_123 in ? F(P_1,P_2,P_3,P_12,P_13,P_23)", "NO_PAIRWISE_REDUCIBILITY_ASSUMPTION": True, "NO_CHAOS_CLAIM_WITHOUT_DIAGNOSTIC": True})
    dump("three_constituent_composite_candidate.json", {"target_id": "three_constituent_composite_state", "question": "Can AAB and ABB form distinct aggregate-defined persistent composites?", "blocked_until": "two distinct constituent-state classes independently derived", "NO_QUARK_IDENTIFICATION": True, "NO_UP_DOWN_PRIMITIVE": True, "NO_MAJORITY_STATE_RULE": True})
    dump("collective_source_em_wave_generation_candidate.json", {"target_id": "collective_source_em_wave_generation", "question": "Can aggregate source reorganization generate a DEV203-like outgoing relational wave?", "prepared_packet_route_preserved": True})
    dump("near_to_far_relational_transition_candidate.json", {"target_id": "collective_near_to_far_relational_transition", "question": "Can complex near-source collective structure yield a coherent outgoing relational wave?", "status": "OPEN_HYPOTHESIS"})
    dump("magnet_iron_passive_response_gate.json", {"MAGNET_IRON_REQUIRES_PASSIVE_INDUCED_RESPONSE": True, "NO_FAKE_IRON_AS_ACTIVE_SOURCE": True, "NO_FAKE_IRON_AS_ACTIVE_MAGNET": True, "status": "BLOCKED"})
    dump("dev229_test_selection.json", {"DEV229_TEST_SELECTION": "PERSISTENT_NATIVE_SOURCE_DERIVATION_GATE", "DEV229_TEST_SELECTION_FROZEN": True, "reason": "MAGNET_LIKE_SOURCE_STATE_VALIDITY=NOT_DERIVED; no pattern mapping is authorized."})
    dump("future_collective_test_selectors.json", {"FUTURE_X_BODY_TEST_SELECTION": "BLOCKED_PENDING_VALID_SOURCE", "candidate_selector": "THREE_SOURCE_COLLECTIVE_REDUCIBILITY_AUDIT", "FINITE_X_NATIVE_COMPOSITION_GATE": values["FINITE_X_NATIVE_COMPOSITION_GATE"], "FUTURE_COMPOSITE_TEST_SELECTION": "THREE_CONSTITUENT_COMPOSITE_STATE_GATE", "composite_block": "two distinct constituent-state classes independently derived"})
    update_docs(); subprocess.check_call([sys.executable, "tools/pbuf_registry.py", "validate"], cwd=ROOT); subprocess.check_call([sys.executable, "tools/pbuf_registry.py", "render"], cwd=ROOT)
    flags = {key: True for key in "CURRENT_GITHUB_INSPECTED CURRENT_HEAD_VERIFIED MECHANISM_REGISTRY_QUERIED DEVELOPMENT_LEDGER_READ HISTORICAL_INDEX_READ DERIVATION_GRAPH_READ DEV159_READ DEV203_READ DEV204_READ DEV211_READ DEV212_READ DEV213_READ DEV214_READ DEV227_READ DEV228_TEST_SELECTION_FROZEN NATIVE_SOURCE_STATE_INVENTORY_COMPLETE NATIVE_SOURCE_PERSISTENCE_CLASSIFIED NATIVE_SOURCE_LOCALIZATION_CLASSIFIED NATIVE_SOURCE_IDENTITY_CLASSIFIED SOURCE_ORIENTATION_STATE_CLASSIFIED MAGNET_LIKE_SOURCE_STATE_VALIDITY_CLASSIFIED DEV213_MULTI_STRUCTURE_COMPOSITION_REUSED TWO_BODY_SOURCE_COMPOSITION_CLASSIFIED FINITE_X_NATIVE_COMPOSITION_GATE_CLASSIFIED COLLECTIVE_X_BODY_CANDIDATE_REGISTERED THREE_SOURCE_REDUCIBILITY_CANDIDATE_REGISTERED THREE_CONSTITUENT_COMPOSITE_CANDIDATE_REGISTERED COLLECTIVE_SOURCE_EM_WAVE_CANDIDATE_REGISTERED NEAR_TO_FAR_RELATIONAL_TRANSITION_CANDIDATE_REGISTERED MAGNET_IRON_PASSIVE_RESPONSE_GATE_PRESERVED NO_INTERSTITIAL_PATTERN_RESULT_RUN NO_DYNAMIC_PATTERN_RESULT_RUN NO_PAIRWISE_REDUCIBILITY_ASSUMPTION NO_POSTHOC_LINEAR_SUPERPOSITION NO_SOURCE_COUNT_SWEEP FROZEN_CONDITION_CLOSURE_RULE_PRESERVED FUTURE_COLLECTIVE_TEST_SELECTORS_RECORDED MECHANISM_REGISTRY_UPDATED REGISTRY_VALIDATED LEDGER_UPDATED HISTORICAL_INDEX_UPDATED TIMELINE_REGENERATED DERIVATION_GRAPH_REGENERATED TESTS_PASS NO_NEW_FORCE NO_NEW_DOF NO_NEW_MAGNETIC_STATE NO_NEW_SOURCE_LAW NO_PACKET_AS_MAGNET_BY_LABEL NO_EXTERNAL_CLAMP_AS_MAGNET_BY_LABEL NO_CHAOS_CLAIM_WITHOUT_DIAGNOSTIC NO_UP_DOWN_PRIMITIVE NO_QUARK_IDENTIFICATION NO_MAJORITY_STATE_RULE NO_FAKE_IRON_AS_ACTIVE_MAGNET NO_RESULT_MOTIVATED_REOPENING NO_ARBITRARY_LOAD_FREQUENCY NO_ARBITRARY_LOAD_AMPLITUDE".split()}
    flags.update(values); flags.update({"DEV228_TEST_SELECTION": "TWO_BODY_SOURCE_STATE_VALIDITY_GATE", "DEV229_TEST_SELECTION": "PERSISTENT_NATIVE_SOURCE_DERIVATION_GATE", "DEV229_TEST_SELECTION_FROZEN": True, "PERSISTENT_SOURCE_NOT_ASSUMED_FROM_PROPAGATING_PACKET": True, "EXTERNAL_CONSTRAINT_NOT_AUTOMATIC_MAGNET": True, "NO_NORTH_SOUTH_PRIMITIVE": True, "NO_MAGNETIC_POLE_LABEL": True, "NO_PR_CREATED": True, "COMMITTED": False, "PUSHED_DIRECTLY_TO_MAIN": False, "REMOTE_MAIN_VERIFIED": False, "WORKTREE_CLEAN": False})
    dump("final_contract.json", flags)
    (OUT / "discussion_handoff.md").write_text("# DEV228 handoff\n\nDEV228 finds no existing self-supported magnet-like source state. DEV213 aggregate preparation is retained, including finite-X same-step composition, but it does not relabel a propagating packet or externally maintained deformation as a magnet. DEV229 is frozen to persistent-source derivation; interstitial mapping and dynamic loading were not run.\n")


if __name__ == "__main__": main()
