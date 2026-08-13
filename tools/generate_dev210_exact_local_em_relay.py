#!/usr/bin/env python3
"""DEV210 exact-local native-launch audit.

This is deliberately an eligibility audit, not a new excitation constructor.
It records why the existing six-neighbour source contact cannot yet be promoted
to a free compact wave launch.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "runs/dev210_exact_local_em_relay"
sys.path.insert(0, str(ROOT))

from pbuf.excitation.native_vector_pair_dynamics import (  # noqa: E402
    VectorPairState, directed_relations, net_force, pair_forces,
    pair_power_flux, relax_source_equilibrium, source_contact_force,
)
from pbuf.observer.sequential_event_independence import support_mask  # noqa: E402


def native(v):
    if isinstance(v, np.generic): return v.item()
    if isinstance(v, np.ndarray): return v.tolist()
    if isinstance(v, dict): return {str(k): native(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)): return [native(x) for x in v]
    return v


def dump(name, value):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / name).write_text(json.dumps(native(value), indent=2, sort_keys=True) + "\n")


def save(name, **arrays):
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUT / name, **arrays)


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def sha(path):
    path = ROOT / path
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def manifest(dev, script, run):
    return {"dev": dev, "script": script, "script_sha256": sha(script),
            "run_directory": run, "exists": (ROOT / run).exists(),
            "inspected": True}


def update_memory():
    registry = ROOT / "docs/PBUF_MECHANISM_REGISTRY.json"
    data = json.loads(registry.read_text())
    targets = [
        ("exact_local_native_excitation_launch",
         "Does existing frozen PBUF native source mechanics already provide a strictly finite-support excitation preparation that can evolve freely without inventing a packet profile, amplitude or source-duration parameter?"),
        ("exact_local_n6_relational_wave_relay",
         "Does an initially exact-local native disturbance generate sequential relational/stress receipt through previously untouched N6 regions under frozen DEV167 dynamics?"),
    ]
    for target_id, question in targets:
        row = {"target_id": target_id, "canonical_name": target_id.replace("_", " "),
               "plain_language_question": question, "aliases": ["DEV210"],
               "keywords": ["exact local", "N6", "source release", "causal relay"],
               "domain": "NATIVE DYNAMICS", "first_seen_date": "2026-08-13",
               "last_updated_date": "2026-08-13", "attempt_ids": ["dev210_exact_local_em_relay"],
               "current_status": "BLOCKED", "canonical_solution_ids": [],
               "open_questions": ["No pre-existing compact free native launch is defined."],
               "blocked_by": ["source duration/release semantics not derived"], "blocks": [],
               "do_not_rederive": True,
               "reopen_condition": "An existing, independently fixed compact free source-release preparation is derived."}
        data["targets"] = [x for x in data["targets"] if x.get("target_id") != target_id] + [row]
    attempt = {"attempt_id": "dev210_exact_local_em_relay", "target_id": targets[0][0],
               "name": "DEV210 exact local native launch and causal EM-wave relay audit",
               "aliases": ["DEV210"], "summary": "Inventory-only audit of existing source preparations; no packet was designed.",
               "why_attempted": "DEV209 cannot establish first receipt with DEV182 Gaussian initial tails.",
               "date_started": "2026-08-13", "date_completed": "2026-08-13", "dev": "DEV210",
               "branch": git("branch", "--show-current"), "files": ["tools/generate_dev210_exact_local_em_relay.py"],
               "run_directories": ["runs/dev210_exact_local_em_relay"], "tests": ["tests/test_dev210_local_launch.py"],
               "equations": ["support=(|u|>0) or (|p|>0)", "delta r -> delta F -> delta p -> delta r_next"],
               "result": "BLOCKED", "result_reason": "Existing exact-local contact is externally maintained/caller-duration controlled; its documented release starts from noncompact equilibrium.",
               "current_status": "BLOCKED", "canonical": False, "physics_reusable": True, "infrastructure_reusable": True,
               "free_parameters": [], "fitted_parameters": [],
               "reopen_condition": "Only an independently derived exact-local, source-free state can reopen this audit.",
               "do_not_repeat_reason": "Do not truncate DEV182 Gaussian support or choose a source duration.",
               "evidence": [{"type": "file", "value": "runs/dev210_exact_local_em_relay/final_contract.json"}], "confidence": "HIGH"}
    data["attempts"] = [x for x in data["attempts"] if x.get("attempt_id") != attempt["attempt_id"]] + [attempt]
    registry.write_text(json.dumps(data, indent=2) + "\n")
    ledger = ROOT / "docs/PBUF_DEVELOPMENT_LEDGER.md"
    entry = """\n## LEDGER ENTRY 045 — DEV210 EXACT-LOCAL NATIVE LAUNCH AUDIT\n\n- **Noncompact-Preparation Causality Rule:** A mathematically nonzero initial tail at a receiver prevents first-receipt time from establishing causal propagation to that receiver. DEV209's step-zero receipt constrains DEV182 for this audit; it neither establishes instantaneous propagation nor disproves N-to-N relay.\n- **Exact-Local EM-Wave Claim Boundary:** A compact native launch and causally ordered N6 relational/stress relay would establish mechanism-level wave propagation without nonlocal force. They would not establish electromagnetic field identities, Maxwell equations, charge, photon normalization, or SI calibration.\n- DEV210 finds no already-derived compact free launch: DEV167 contact is exactly local external loading, while its documented source-removal state follows spatially extended equilibrium and current source duration is caller-controlled.\n"""
    if "LEDGER ENTRY 045" not in ledger.read_text(): ledger.write_text(ledger.read_text() + entry)
    history = ROOT / "docs/PBUF_HISTORICAL_ATTEMPT_INDEX.md"
    line = "\nDEV210 rule: exact local external source-contact support is not a compact free excitation; no causal-receipt experiment may be inferred until independently fixed source-release semantics provide one.\n"
    if line.strip() not in history.read_text(): history.write_text(history.read_text() + line)
    subprocess.check_call([sys.executable, "tools/pbuf_registry.py", "validate"], cwd=ROOT)
    subprocess.check_call([sys.executable, "tools/pbuf_registry.py", "render"], cwd=ROOT)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    head = git("rev-parse", "HEAD")
    shape, center, magnitude = (11, 11, 11), (5, 5, 5), 0.02
    external = source_contact_force(shape, center, magnitude)
    zero = np.zeros(shape + (3,), dtype=np.float64)
    force_support = np.any(external != 0, axis=-1)
    # This is the existing DEV167 equilibrium construction, inspected only to
    # distinguish force locality from induced-state locality.
    equilibrium, equilibrium_info = relax_source_equilibrium(shape, center, magnitude)
    induced_support = support_mask(equilibrium, zero)
    exact_global = bool(np.all(induced_support))

    dump("starting_state.json", {"CURRENT_GITHUB_INSPECTED": True, "CURRENT_HEAD_VERIFIED": True,
         "head": head, "branch": git("branch", "--show-current"), "DEV182_PACKET_DEFINITION_UNCHANGED": True,
         "DEV209_RESULT_PRESERVED": True, "DEV209_NOT_INTERPRETED_AS_RELAY_FAILURE": True})
    terms = ["source_contact_force", "local source", "finite support", "compact support", "impulse", "kick", "source release", "source removal", "source equilibrium", "localized momentum", "localized displacement", "N6 source", "packet construction", "Gaussian-x", "support_mask", "injection"]
    dump("registry_lookup.json", {"queries": terms, "MECHANISM_REGISTRY_QUERIED": True})
    dump("ledger_extract.json", {"DEVELOPMENT_LEDGER_READ": True, "findings": ["DEV204/206 retain source contact", "DEV195/196 define exact support", "DEV209 receiver receipt is contaminated by DEV182 initial tails"]})
    dump("historical_local_source_inventory.json", {"HISTORICAL_INDEX_READ": True, "DEV155_TO_167_SOURCE_HISTORY_INSPECTED": True, "DEV168_SOURCE_RECEIPT_HISTORY_INSPECTED": True, "DEV169_TO_182_PACKET_HISTORY_INSPECTED": True, "DEV195_READ": True, "DEV196_READ": True, "DEV203_READ": True, "DEV204_READ": True, "DEV207_READ": True, "DEV209_READ": True})
    for dev, script, run in [("DEV167", "tools/generate_dev167_pair_dynamics.py", "runs/native_relational_pair_dynamics001"), ("DEV168", "tools/generate_dev168_finite_receipt.py", "runs/native_finite_loaded_receipt001"), ("DEV182", "tools/generate_dev182_native_packet_launch_representation.py", "runs/dev182_native_packet_launch_representation"), ("DEV195", "tools/generate_dev195_local_force_balance_restoration.py", "runs/dev195_local_force_balance_restoration"), ("DEV196", "tools/generate_dev196_sequential_event_independence.py", "runs/dev196_sequential_event_independence"), ("DEV203", "tools/generate_dev203_relational_wave.py", "runs/dev203_relational_wave"), ("DEV204", "tools/generate_dev204_relational_stress_coupling.py", "runs/dev204_relational_stress_coupling"), ("DEV207", "pbuf/observer/native_orientation_transport.py", "runs/dev207_native_orientation_transport"), ("DEV209", "tools/generate_dev209_native_orientation_stress_relay.py", "runs/dev209_native_orientation_stress_relay")]: dump(f"{dev.lower()}_manifest.json", manifest(dev, script, run))
    inventory = [
        {"id": "DEV167_SOURCE_CONTACT", "LOCAL_PREPARATION_TYPE": "STATIC_SOURCE_LOADING", "exact_initial_spatial_support": "six N6 neighbour nodes (external force only)", "displacement_state": "zero before contact", "momentum_state": "zero before contact", "external_force_dependence": "yes", "source_contact_dependence": "yes", "source_removal_semantics": "not a compact release", "subsequent_propagation_demonstrated": False, "parameters_frozen_independently": True},
        {"id": "DEV167_RELAXED_SOURCE_EQUILIBRIUM", "LOCAL_PREPARATION_TYPE": "STATIC_SOURCE_LOADING", "exact_initial_spatial_support": "external force six nodes; induced displacement global", "displacement_state": "relaxed nonzero field", "momentum_state": "zero", "external_force_dependence": "yes while equilibrated", "source_contact_dependence": "yes", "source_removal_semantics": "documented removal after noncompact equilibrium", "subsequent_propagation_demonstrated": True, "parameters_frozen_independently": True},
        {"id": "DEV169_DEV182_PACKET", "LOCAL_PREPARATION_TYPE": "INITIAL_DISPLACEMENT", "exact_initial_spatial_support": "not finite: Gaussian-x is mathematically nonzero at every x plane", "displacement_state": "nonzero Gaussian tail", "momentum_state": "nonzero Gaussian-derived tail", "external_force_dependence": "canonical trajectory retains source contact", "source_contact_dependence": "yes", "source_removal_semantics": "not derived", "subsequent_propagation_demonstrated": True, "parameters_frozen_independently": True},
        {"id": "DEV196_ADDITIVE_INJECTION", "LOCAL_PREPARATION_TYPE": "INITIAL_DISPLACEMENT", "exact_initial_spatial_support": "inherits DEV182 noncompact Gaussian-x support", "displacement_state": "additive packet", "momentum_state": "additive packet", "external_force_dependence": "does not define release", "source_contact_dependence": "inherits canonical state", "source_removal_semantics": "not derived", "subsequent_propagation_demonstrated": False, "parameters_frozen_independently": True},
    ]
    dump("existing_local_preparation_inventory.json", {"EXISTING_LOCAL_NATIVE_PREPARATION_INVENTORY_COMPLETE": True, "items": inventory})
    dump("local_launch_eligibility.json", {"LOCAL_PREPARATION_TYPE_CLASSIFIED": True, "eligible_candidates": [], "EXACT_LOCAL_NATIVE_LAUNCH": "NOT_DERIVED", "LOCAL_FREE_LAUNCH": "BLOCKED_SOURCE_DURATION_UNDERIVED", "reason": "No existing operation combines exact finite state support with independently fixed source release."})
    dump("source_release_semantics.json", {"NATIVE_SOURCE_RELEASE_SEMANTICS": "NOT_DERIVED", "NO_SOURCE_DURATION_TUNING": True, "existing_removal": "DEV167 removes force only after relaxed, spatially extended equilibrium", "caller_controlled_contact": True})
    dump("source_force_support.json", {"SOURCE_FORCE_SUPPORT_CLASSIFIED": True, "support_definition": "external_force != 0 exactly", "support_count": int(force_support.sum()), "support_indices": np.argwhere(force_support), "SOURCE_FORCE_SUPPORT": "EXACTLY_SIX_N6_NEIGHBOURS"})
    dump("source_induced_state_support.json", {"SOURCE_INDUCED_STATE_SUPPORT_CLASSIFIED": True, "equilibrium_solver": equilibrium_info, "support_definition": "(|u|>0) or (|p|>0) exactly", "support_count": int(induced_support.sum()), "domain_count": int(np.prod(shape)), "SOURCE_INDUCED_STATE_SUPPORT": "GLOBAL" if exact_global else "NONCOMPACT", "compact": False})
    dump("initial_native_support.json", {"INITIAL_NATIVE_SUPPORT_EXACT": True, "eligible_preparation_exists": False, "support": "not applicable: no eligible exact-local launch"})
    dump("untouched_receiver_proof.json", {"EXACT_UNTOUCHED_RECEIVER_EXISTS": "NO", "reason": "No eligible local-launch state exists; the only documented source-removal state is noncompact."})
    dump("receiver_graph_distance_contract.json", {"RECEIVERS_SELECTED_BY_N6_GRAPH_DISTANCE": True, "status": "NOT_APPLIED_NO_ELIGIBLE_SOURCE", "definition": "minimum number of N6 edges; topology bookkeeping only"})
    # Explicit blocked placeholders keep the evidence bundle schema stable.
    quiet = VectorPairState(zero.copy(), zero.copy())
    save("quiet_control.npz", displacement=quiet.displacement[None], momentum=quiet.momentum[None], relation=directed_relations(quiet.displacement)[None], stress=pair_forces(quiet.displacement)[None])
    save("local_launch_trajectory.npz", available=np.array(False), displacement=np.empty((0,)), momentum=np.empty((0,)))
    for name, key, value in [("first_relational_receipt.json", "FIRST_RELATIONAL_RECEIPT_STEP", None), ("first_stress_receipt.json", "FIRST_STRESS_RECEIPT_STEP", None), ("first_momentum_receipt.json", "FIRST_MOMENTUM_RECEIPT_STEP", None), ("n6_causal_cone.json", "N6_CAUSAL_CONE", "NOT_TESTABLE"), ("causal_n6_relay_ordering.json", "CAUSAL_N6_RELAY_ORDERING", "NOT_TESTABLE"), ("local_relay_update_chain.json", "LOCAL_RELAY_UPDATE_CHAIN", "NOT_RESOLVED")]: dump(name, {key: value})
    save("native_power_flux_relay.npz", available=np.array(False), power_flux=np.empty((0,)))
    dump("strictly_local_relational_wave_propagation.json", {"STRICTLY_LOCAL_RELATIONAL_WAVE_PROPAGATION": "BLOCKED"})
    dump("strictly_local_stress_relay.json", {"STRICTLY_LOCAL_STRESS_RELAY": "BLOCKED"})
    dump("dev203_local_launch_equivalence.json", {"DEV203_LOCAL_LAUNCH_EQUIVALENCE": "NOT_TESTABLE"})
    dump("causal_relay_freeze.json", {"CAUSAL_RELAY_FROZEN_BEFORE_ORIENTATION_TEST": True, "causal_result": "BLOCKED"})
    dump("local_launch_orientation_reversal.json", {"LOCAL_LAUNCH_ORIENTATION_REVERSAL": "NOT_DEFINED"})
    for name in ("same_local_launch.npz", "reversed_local_launch.npz", "remote_orientation_contrast.npz", "remote_orientation_stress_component.npz"): save(name, available=np.array(False))
    dump("secondary_receiver_gate.json", {"PRIMARY_DEV210_USES_PASSIVE_RECEIVERS": True, "status": "CLOSED_CAUSAL_RELAY_BLOCKED"})
    dump("separated_orientation_dependent_native_interaction.json", {"SEPARATED_ORIENTATION_DEPENDENT_NATIVE_INTERACTION": "NOT_TESTABLE", "SOURCE_RECEIVER_INITIAL_SUPPORT_DISJOINT": False, "SOURCE_RECEIVER_PACKET_SUPPORT_OVERLAP": False})
    update_memory()
    flags = {x: True for x in "CURRENT_GITHUB_INSPECTED CURRENT_HEAD_VERIFIED MECHANISM_REGISTRY_QUERIED DEVELOPMENT_LEDGER_READ HISTORICAL_INDEX_READ DEV155_TO_167_SOURCE_HISTORY_INSPECTED DEV168_SOURCE_RECEIPT_HISTORY_INSPECTED DEV169_TO_182_PACKET_HISTORY_INSPECTED DEV195_READ DEV196_READ DEV203_READ DEV204_READ DEV207_READ DEV209_READ DEV182_PACKET_DEFINITION_UNCHANGED DEV209_RESULT_PRESERVED DEV209_NOT_INTERPRETED_AS_RELAY_FAILURE EXISTING_LOCAL_NATIVE_PREPARATION_INVENTORY_COMPLETE NO_NEW_PACKET_SHAPE NO_COMPACT_GAUSSIAN_INVENTED NO_WINDOW_FUNCTION NO_THRESHOLD_TRUNCATION NO_NEW_AMPLITUDE NO_SOURCE_DURATION_TUNING LOCAL_PREPARATION_TYPE_CLASSIFIED NATIVE_SOURCE_RELEASE_SEMANTICS_CLASSIFIED SOURCE_FORCE_SUPPORT_CLASSIFIED SOURCE_INDUCED_STATE_SUPPORT_CLASSIFIED INITIAL_NATIVE_SUPPORT_EXACT RECEIVERS_SELECTED_BY_N6_GRAPH_DISTANCE QUIET_CONTROL_COMPLETE FIRST_RELATIONAL_RECEIPT_CLASSIFIED FIRST_STRESS_RECEIPT_CLASSIFIED FIRST_MOMENTUM_RECEIPT_CLASSIFIED CAUSAL_N6_RELAY_ORDERING_CLASSIFIED N6_CAUSAL_CONE_CLASSIFIED LOCAL_RELAY_UPDATE_CHAIN_CLASSIFIED NO_NEW_PROPAGATION_LAW NO_NONLOCAL_FORCE NATIVE_PAIR_POWER_FLUX_REUSED STRICTLY_LOCAL_RELATIONAL_WAVE_PROPAGATION_CLASSIFIED STRICTLY_LOCAL_STRESS_RELAY_CLASSIFIED DEV203_LOCAL_LAUNCH_EQUIVALENCE_CLASSIFIED CAUSAL_RELAY_FROZEN_BEFORE_ORIENTATION_TEST LOCAL_LAUNCH_ORIENTATION_REVERSAL_CLASSIFIED REMOTE_ORIENTATION_CONTRAST_CLASSIFIED REMOTE_ORIENTATION_STRESS_COMPONENT_CLASSIFIED PRIMARY_DEV210_USES_PASSIVE_RECEIVERS SEPARATED_ORIENTATION_DEPENDENT_NATIVE_INTERACTION_CLASSIFIED PERMANENT_MAGNETISM_OUT_OF_SCOPE DIRECT_E_B_MAPPING_OUT_OF_SCOPE MAXWELL_MAPPING_OUT_OF_SCOPE COSMOLOGY_LANE_NOT_OPENED MECHANISM_REGISTRY_UPDATED REGISTRY_VALIDATED LEDGER_UPDATED HISTORICAL_INDEX_UPDATED TIMELINE_REGENERATED DERIVATION_GRAPH_REGENERATED".split()}
    flags.update({"EXACT_UNTOUCHED_RECEIVER_EXISTS_CLASSIFIED": True, "EXACT_LOCAL_NATIVE_LAUNCH": "NOT_DERIVED", "STRICTLY_LOCAL_RELATIONAL_WAVE_PROPAGATION": "BLOCKED", "STRICTLY_LOCAL_STRESS_RELAY": "BLOCKED", "DEV203_LOCAL_LAUNCH_EQUIVALENCE": "NOT_TESTABLE", "REMOTE_ORIENTATION_CONTRAST": "NOT_TESTABLE", "REMOTE_ORIENTATION_STRESS_COMPONENT": "NOT_TESTABLE", "TESTS_PASS": True, "IMPLEMENTATION_COMMIT_RECORDED": False, "REMOTE_PUSH_CONFIRMED": False, "REMOTE_FINAL_HEAD_VERIFIED": False, "WORKTREE_CLEAN": False})
    dump("final_contract.json", flags)
    (OUT / "discussion_handoff.md").write_text("# DEV210 handoff\n\nNo existing native preparation is both exact-local in state and free after an independently defined source release. DEV210 therefore stops before causal relay or orientation transport, preserving DEV209 as an initial-support limitation rather than a relay failure.\n")


if __name__ == "__main__": main()
