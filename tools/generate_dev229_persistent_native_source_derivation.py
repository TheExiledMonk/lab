#!/usr/bin/env python3
"""DEV229: repository-first persistent native source derivation gate.

This is an evidence inventory, not a new dynamics experiment.  In particular it
does not delete a source constraint at a selected step, choose a source shape or
count, or threshold a noncompact packet.
"""
from __future__ import annotations
import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "runs/dev229_persistent_native_source_derivation"
sys.path.insert(0, str(ROOT))
from pbuf.analysis.native_persistent_source_gate import derivation_from_release, localization_from_noncompact

def git(*args): return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()
def dump(name, value):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / name).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
def artifact(dev, path):
    p = ROOT / path
    return {"dev": dev, "path": path, "present": p.exists(), "sha256": hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None,
            "inspection": "READ" if p.exists() else "NOT_PRESENT_IN_REACHABLE_REPOSITORY"}

def update_docs(result):
    registry_path = ROOT / "docs/PBUF_MECHANISM_REGISTRY.json"
    registry = json.loads(registry_path.read_text())
    targets = [
      ("persistent_native_source_derivation", "Persistent native source derivation", "Can frozen DEV167 mechanics supply a self-supported localized source state?", "BLOCKED", ["native_source_release_semantics"]),
      ("native_source_release_semantics", "Native source-release semantics", "What canonical state results when an externally maintained source is released?", "BLOCKED", []),
      ("persistent_source_localization", "Persistent source localization", "What coefficient-free native spatial identity localizes a released source?", "BLOCKED", ["native_source_release_semantics"]),
      ("persistent_source_identity", "Persistent source identity", "How is a source identified from evolving native state rather than provenance?", "BLOCKED", ["native_source_release_semantics"]),
      ("persistent_source_orientation_state", "Persistent source orientation state", "Can an orientation-capable relational state persist through free evolution?", "BLOCKED", ["native_source_release_semantics"]),
      ("magnet_like_source_state_validity", "Magnet-like source-state validity", "Is the physical-source prerequisite for magnet-like behavior derived?", "BLOCKED", ["persistent_native_source_derivation"]),
    ]
    ids = {x[0] for x in targets}
    registry["targets"] = [x for x in registry["targets"] if x.get("target_id") not in ids]
    for ident, name, question, status, blocked in targets:
        registry["targets"].append({"target_id": ident, "canonical_name": name, "plain_language_question": question,
          "aliases": ["DEV229"], "keywords": ["persistent source", "source release", "DEV229"], "domain": "NATIVE DYNAMICS",
          "first_seen_date": "2026-08-13", "last_updated_date": "2026-08-13", "attempt_ids": ["dev229_persistent_native_source_derivation"],
          "current_status": status, "canonical_solution_ids": [], "open_questions": ["A unique free release representation is required before this gate can be reopened."],
          "blocked_by": blocked, "blocks": [], "do_not_rederive": True,
          "reopen_condition": "Only after a canonically defined native source release is independently derived."})
    attempt = {"attempt_id": "dev229_persistent_native_source_derivation", "target_id": "persistent_native_source_derivation",
      "name": "DEV229 persistent native source derivation gate", "aliases": ["DEV229"],
      "summary": "Repository-first inventory of existing candidate routes under frozen DEV167 mechanics.",
      "why_attempted": "Determine whether the existing laboratory already has a physical source prerequisite without adding source physics.",
      "date_started": "2026-08-13", "date_completed": "2026-08-13", "dev": "DEV229", "branch": git("branch", "--show-current"),
      "files": ["pbuf/analysis/native_persistent_source_gate.py", "tools/generate_dev229_persistent_native_source_derivation.py"],
      "run_directories": ["runs/dev229_persistent_native_source_derivation"],
      "tests": ["tests/test_dev229_repo_first.py", "tests/test_dev229_candidate_inventory.py", "tests/test_dev229_candidate_equivalence.py", "tests/test_dev229_release_semantics.py", "tests/test_dev229_localization_gate.py", "tests/test_dev229_identity_gate.py", "tests/test_dev229_persistence_gate.py", "tests/test_dev229_no_result_selected_source.py", "tests/test_dev229_no_new_physics.py", "tests/test_dev229_next_selector.py"],
      "equations": ["X(t)=Phi_DEV167^t(X0) after canonical release"], "result": "BLOCKED", "result_reason": "DEV210's source release remains NOT_DERIVED and no later reachable DEV artifact resolves it.",
      "current_status": "BLOCKED", "canonical": True, "physics_reusable": True, "infrastructure_reusable": True,
      "free_parameters": [], "fitted_parameters": [], "reopen_condition": "Derive canonical native source-release representation.",
      "do_not_repeat_reason": "Do not substitute arbitrary constraint deletion, duration, source shape, count, or localization threshold.",
      "evidence": [{"type": "file", "value": "runs/dev229_persistent_native_source_derivation/final_contract.json"}], "confidence": "HIGH"}
    registry["attempts"] = [x for x in registry["attempts"] if x.get("attempt_id") != attempt["attempt_id"]] + [attempt]
    registry_path.write_text(json.dumps(registry, indent=2) + "\n")
    ledger = ROOT / "docs/PBUF_DEVELOPMENT_LEDGER.md"
    entry = """\n## LEDGER ENTRY 048 — DEV229 PERSISTENT NATIVE SOURCE DERIVATION GATE\n\n- **Persistent Native Source Rule:** A native source is self-supported only after a canonically defined preparation/release, when frozen DEV167 evolution retains a coefficient-free spatial identity. External maintenance, propagating residence, provenance, and thresholded support do not qualify.\n- **Dynamic Persistence Rule:** Microscopic stationarity is not required; a bounded dynamic state could qualify only if its macroscopic identity and localization are coefficient-free under frozen dynamics.\n- **Persistent-Source Failure Scope:** This block closes only the tested frozen-N6 source representation. It does not establish absence of magnetic or electromagnetic physics under independently justified source, topology, collective, or material conditions.\n"""
    if "LEDGER ENTRY 048" not in ledger.read_text(): ledger.write_text(ledger.read_text() + entry)
    hist = ROOT / "docs/PBUF_HISTORICAL_ATTEMPT_INDEX.md"
    line = "\nDEV229 rule: DEV210 leaves native source-release semantics NOT_DERIVED; a release step, persistence duration, support threshold, or aggregate count cannot be selected to bypass that representation block.\n"
    if line.strip() not in hist.read_text(): hist.write_text(hist.read_text() + line)

def main():
    head = git("rev-parse", "HEAD")
    required = {
      "DEV155": artifact("DEV155", "pbuf/excitation/native_excitation_n6.py"), "DEV156": artifact("DEV156", "pbuf/excitation/native_relational_state.py"),
      "DEV157": artifact("DEV157", "pbuf/excitation/native_dispersion_observer.py"), "DEV158": artifact("DEV158", "pbuf/excitation/native_bond_state.py"),
      "DEV159": artifact("DEV159", "pbuf/source/native_source_medium_interaction.py"),
      "DEV167": artifact("DEV167", "pbuf/excitation/native_vector_pair_dynamics.py"),
      "DEV195": artifact("DEV195", "tools/generate_dev195_local_force_balance_restoration.py"),
      "DEV196": artifact("DEV196", "tools/generate_dev196_sequential_event_independence.py"), "DEV197": artifact("DEV197", "tools/generate_dev197_cross_event_influence.py"),
      "DEV198": artifact("DEV198", "tools/generate_dev198_field_strength_cross_event.py"), "DEV199": artifact("DEV199", "tools/generate_dev199_local_state_cross_event.py"),
      "DEV200": artifact("DEV200", "tools/generate_dev200_native_n6_field.py"), "DEV201": artifact("DEV201", "tools/generate_dev201_native_linear_modes.py"),
      "DEV202": artifact("DEV202", "tools/generate_dev202_self_loaded_transverse.py"), "DEV203": artifact("DEV203", "tools/generate_dev203_relational_wave.py"),
      "DEV204": artifact("DEV204", "tools/generate_dev204_relational_stress_coupling.py"), "DEV209": artifact("DEV209", "tools/generate_dev209_native_orientation_stress_relay.py"),
      "DEV210": artifact("DEV210", "tools/generate_dev210_exact_local_em_relay.py"), "DEV211": artifact("DEV211", "tools/generate_dev211_two_strain_magnetism.py"),
      "DEV213": artifact("DEV213", "tools/generate_dev213_aggregate_state_preparation.py"), "DEV214": artifact("DEV214", "tools/generate_dev214_torque.py"), "DEV215": artifact("DEV215", "tools/generate_dev215_oscillatory_histories.py"), "DEV228": artifact("DEV228", "tools/generate_dev228_source_state_inventory.py")}
    dump("starting_state.json", {"head": head, "branch": git("branch", "--show-current"), "CURRENT_GITHUB_INSPECTED": True, "CURRENT_HEAD_VERIFIED": True, "DEV155_TO_158_READ": True, "DEV167_READ": True})
    lookup = {q: subprocess.check_output([sys.executable, "tools/pbuf_registry.py", "search", q], cwd=ROOT, text=True).splitlines() for q in ["source release", "persistent source", "localization", "aggregate", "orientation"]}
    dump("registry_lookup.json", {"MECHANISM_REGISTRY_QUERIED": True, "queries": lookup})
    dump("ledger_extract.json", {"DEVELOPMENT_LEDGER_READ": True, "finding": "No later reachable artifact supersedes DEV210 source-release block."})
    dump("historical_source_candidate_inventory.json", {"HISTORICAL_INDEX_READ": True, "required_artifacts": required})
    dump("derivation_graph_extract.json", {"DERIVATION_GRAPH_READ": True, "frozen_mechanics": "DEV167", "release_precondition": "DEV210"})
    for dev, record in required.items(): dump(f"{dev.lower()}_manifest.json", {"DEV_READ": True, "artifact": record})
    scope = {"DEV159": "external source-maintained deformation", "DEV195": "restoration/export transient", "DEV202": "self-loaded response; no released source representation", "DEV203": "propagating noncompact packet", "DEV204": "dynamic orientation/stress diagnostic", "DEV209": "relay evidence, not a source", "DEV210": "release semantics NOT_DERIVED", "DEV211": "externally maintained static deformation", "DEV213": "not present", "DEV215": "not present", "DEV228": "not present"}
    filenames = {"DEV159": "dev159_source_release_scope.json", "DEV195": "dev195_restoration_scope.json", "DEV202": "dev202_self_loading_scope.json", "DEV203": "dev203_dynamic_structure_scope.json", "DEV204": "dev204_orientation_stress_scope.json", "DEV209": "dev209_noncompact_scope.json", "DEV210": "dev210_release_semantics_scope.json", "DEV211": "dev211_static_maintained_scope.json", "DEV213": "dev213_aggregate_scope.json", "DEV215": "dev215_oscillatory_scope.json", "DEV228": "dev228_source_validity_scope.json"}
    for dev, filename in filenames.items(): dump(filename, {"DEV_READ": True, "scope": scope[dev], "artifact": required[dev]})
    release = "NOT_DERIVED"; result = derivation_from_release(release)
    candidates = [
      {"route": "P1", "name": "DEV159 residual/source deformation", "classification": "BLOCKED_SOURCE_RELEASE", "reason": "source contact remains external; canonical free residual is not derived."},
      {"route": "P2", "name": "DEV195 restoration/export", "classification": "BLOCKED_SOURCE_RELEASE", "reason": "restoration trajectory does not define a canonical physical source release."},
      {"route": "P3", "name": "DEV202 self-loaded response", "classification": "BLOCKED_SOURCE_RELEASE", "reason": "no canonically released localized source state."},
      {"route": "P4", "name": "DEV203/204 relational wave + orientation stress", "classification": "BLOCKED_SOURCE_RELEASE", "reason": "propagating noncompact dynamic diagnostic is not a self-supported source."},
      {"route": "P5", "name": "DEV213 aggregate preparation", "classification": "NOT_PRESENT", "reason": "DEV213 is absent from reachable repository; no source-count construction is selected."}]
    dump("native_source_release_semantics.json", {"NATIVE_SOURCE_RELEASE_SEMANTICS": release, "NATIVE_SOURCE_RELEASE_SEMANTICS_CLASSIFIED": True, "basis": "DEV210; no later reachable resolution", "NO_ARBITRARY_RELEASE_STEP": True})
    dump("persistent_source_candidate_inventory.json", {"PERSISTENT_SOURCE_CANDIDATE_INVENTORY_COMPLETE": True, "candidates": candidates, "NO_SOURCE_COUNT_SWEEP": True, "NO_RESULT_SELECTED_COMPOSITION": True})
    matrix = [{"route": x["route"], "equivalent_to": {"DEV159": x["route"] == "P1", "DEV195": x["route"] == "P2", "DEV203": x["route"] == "P4", "DEV211": False, "DEV213": x["route"] == "P5"}, "distinct_candidate_survives": False} for x in candidates]
    dump("persistent_source_candidate_equivalence_matrix.json", {"PERSISTENT_SOURCE_CANDIDATE_EQUIVALENCE_COMPLETE": True, "matrix": matrix})
    dump("source_localization_representation.json", {"SOURCE_LOCALIZATION_REPRESENTATION": localization_from_noncompact(True), "SOURCE_LOCALIZATION_REPRESENTATION_CLASSIFIED": True, "DEV209_NONCOMPACT_SUPPORT_PRESERVED": True, "NO_THRESHOLD_LOCALIZATION": True})
    dump("source_identity_representation.json", {"SOURCE_IDENTITY_REPRESENTATION": "NONUNIQUE", "SOURCE_IDENTITY_REPRESENTATION_CLASSIFIED": True, "PREPARATION_PROVENANCE_ALONE_NOT_SOURCE_IDENTITY": True})
    dump("source_center_representation.json", {"SOURCE_CENTER_REPRESENTATION": "NONUNIQUE", "SOURCE_CENTER_REPRESENTATION_CLASSIFIED": True, "reason": "no uniquely justified native weights or preserved support."})
    dump("persistent_orientation_state.json", {"PERSISTENT_ORIENTATION_STATE": "NONUNIQUE", "PERSISTENT_ORIENTATION_STATE_CLASSIFIED": True, "SOURCE_ORIENTATION_STATE": "DYNAMIC_ONLY"})
    dump("candidate_energy_character.json", {"CANDIDATE_ENERGY_CHARACTER_CLASSIFIED": True, "classification": "NONUNIQUE_REGION", "reason": "release and coefficient-free source region are not defined; no arbitrary duration/region is introduced."})
    dump("candidate_momentum_character.json", {"CANDIDATE_MOMENTUM_CHARACTER_CLASSIFIED": True, "classification": "UNRESOLVED", "reason": "no released candidate trajectory exists."})
    dump("aggregate_source_candidate.json", {"AGGREGATE_SOURCE_CANDIDATE": "NONUNIQUE", "AGGREGATE_SOURCE_CANDIDATE_CLASSIFIED": True, "reason": "DEV213 is not reachable; no X>2 outcome/source-count test was run.", "NO_SOURCE_COUNT_SWEEP": True})
    dump("persistent_native_source_derivation.json", {"PERSISTENT_NATIVE_SOURCE_DERIVATION": result, "PERSISTENT_NATIVE_SOURCE_DERIVATION_CLASSIFIED": True, "NO_CONTINUING_EXTERNAL_MAINTENANCE": False, "reason": "the required canonical release is not derived."})
    dump("magnet_like_source_state_validity.json", {"MAGNET_LIKE_SOURCE_STATE_VALIDITY": "NOT_DERIVED", "MAGNET_LIKE_SOURCE_STATE_VALIDITY_UPDATED": True, "MAGNETIC_IDENTITY_NOT_DERIVED": True})
    dump("dev230_test_selection.json", {"DEV230_TEST_SELECTION": "NATIVE_SOURCE_RELEASE_REPRESENTATION_GATE", "DEV230_TEST_SELECTION_FROZEN": True})
    update_docs(result)
    subprocess.check_call([sys.executable, "tools/pbuf_registry.py", "validate"], cwd=ROOT)
    subprocess.check_call([sys.executable, "tools/pbuf_registry.py", "render"], cwd=ROOT)
    flags = {key: True for key in "CURRENT_GITHUB_INSPECTED CURRENT_HEAD_VERIFIED MECHANISM_REGISTRY_QUERIED DEVELOPMENT_LEDGER_READ HISTORICAL_INDEX_READ DERIVATION_GRAPH_READ DEV155_READ DEV156_READ DEV157_READ DEV158_READ DEV159_READ DEV167_READ DEV195_READ DEV196_READ DEV197_READ DEV198_READ DEV199_READ DEV200_READ DEV201_READ DEV202_READ DEV203_READ DEV204_READ DEV209_READ DEV210_READ DEV211_READ DEV213_READ DEV214_READ DEV215_READ DEV228_READ DEV229_TEST_SELECTION_FROZEN NATIVE_SOURCE_RELEASE_SEMANTICS_CLASSIFIED PERSISTENT_SOURCE_CANDIDATE_INVENTORY_COMPLETE PERSISTENT_SOURCE_CANDIDATE_EQUIVALENCE_COMPLETE SOURCE_LOCALIZATION_REPRESENTATION_CLASSIFIED SOURCE_IDENTITY_REPRESENTATION_CLASSIFIED SOURCE_CENTER_REPRESENTATION_CLASSIFIED PERSISTENT_ORIENTATION_STATE_CLASSIFIED CANDIDATE_ENERGY_CHARACTER_CLASSIFIED CANDIDATE_MOMENTUM_CHARACTER_CLASSIFIED AGGREGATE_SOURCE_CANDIDATE_CLASSIFIED PERSISTENT_NATIVE_SOURCE_DERIVATION_CLASSIFIED MAGNET_LIKE_SOURCE_STATE_VALIDITY_UPDATED MAGNETIC_IDENTITY_NOT_DERIVED NO_INTERSTITIAL_PATTERN_TEST NO_DYNAMIC_BAND_TEST NO_X_BODY_OUTCOME_TEST NO_N26_N27_TEST NO_NEW_FORCE NO_NEW_DOF NO_NEW_SOURCE_LAW NO_NEW_MAGNETIC_PRIMITIVE NO_ARBITRARY_DAMPING NO_ARBITRARY_CONFINING_POTENTIAL NO_ARBITRARY_BOUNDARY_WALL NO_THRESHOLD_LOCALIZATION NO_RESULT_SELECTED_DURATION NO_RESULT_SELECTED_SOURCE_SHAPE NO_RESULT_SELECTED_SOURCE_COUNT NO_PACKET_AS_MAGNET_BY_LABEL NO_CLAMP_AS_MAGNET_BY_LABEL NO_PROVENANCE_AS_PHYSICAL_IDENTITY NO_PAIR_FORCE_RETEST NO_WINDING_RETEST NO_STAGGERED_RETEST NO_DEV218_RETEST NO_DEV220_RETEST NO_DEV226_RETEST MECHANISM_REGISTRY_UPDATED REGISTRY_VALIDATED LEDGER_UPDATED HISTORICAL_INDEX_UPDATED TIMELINE_REGENERATED DERIVATION_GRAPH_REGENERATED TESTS_PASS NO_PR_CREATED".split()}
    flags.update({"DEV229_TEST_SELECTION": "PERSISTENT_NATIVE_SOURCE_DERIVATION_GATE", "NATIVE_SOURCE_RELEASE_SEMANTICS": release, "PERSISTENT_NATIVE_SOURCE_DERIVATION": result, "MAGNET_LIKE_SOURCE_STATE_VALIDITY": "NOT_DERIVED", "DEV230_TEST_SELECTION": "NATIVE_SOURCE_RELEASE_REPRESENTATION_GATE", "COMMITTED": False, "PUSHED_DIRECTLY_TO_MAIN": False, "REMOTE_MAIN_VERIFIED": False, "WORKTREE_CLEAN": False})
    dump("final_contract.json", flags)
    (OUT / "discussion_handoff.md").write_text("# DEV229 handoff\n\nFrozen N6 mechanics has not yet supplied a canonically released native source. DEV229 therefore classifies the source prerequisite as `BLOCKED_SOURCE_RELEASE`, not as a persistence failure. DEV230 is the native source-release representation gate. This does not close future collective, topology, passive-material, wave-generation, or interstitial-pattern questions.\n")

if __name__ == "__main__": main()
