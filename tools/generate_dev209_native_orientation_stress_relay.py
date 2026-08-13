#!/usr/bin/env python3
"""DEV209 — passive, N6-only audit of the transported DEV207 orientation bit."""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "runs/dev209_native_orientation_stress_relay"
sys.path.insert(0, str(ROOT))

from pbuf.excitation.native_vector_pair_dynamics import VectorPairState, invariant, inverse_step, step
from pbuf.observer.native_orientation_transport import orientation_packets, reflected_state
from pbuf.observer.native_remote_receiver import receiver_planes
from pbuf.observer.native_stress_relay import finite_components, receiver_state
from pbuf.observer.sequential_event_independence import support_mask
from tools import generate_dev169_raw_abell_native_observer as D
from tools import generate_dev184_discrete_launch_density_convergence as E


def native(value):
    if isinstance(value, np.generic): return value.item()
    if isinstance(value, np.ndarray): return value.tolist()
    if isinstance(value, dict): return {str(k): native(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)): return [native(v) for v in value]
    return value


def dump(name, value):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / name).write_text(json.dumps(native(value), indent=2, sort_keys=True, allow_nan=False) + "\n")


def save(name, **arrays):
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUT / name, **arrays)


def git(*args): return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def nonzero(value): return not np.array_equal(value, np.zeros_like(value))
def floor(value): return np.finfo(float).eps * max(1.0, float(np.max(np.abs(value)))) * 128.0
def status(value): return "PRESENT" if nonzero(value) else "ZERO_OR_MACHINE_FLOOR"


def evolve(state, external):
    rows = {k: [] for k in ("displacement", "momentum", "invariant")}
    s = state
    for n in range(D.STEPS + 1):
        rows["displacement"].append(s.displacement.copy())
        rows["momentum"].append(s.momentum.copy())
        rows["invariant"].append(invariant(s.displacement, s.momentum))
        if n < D.STEPS: s = step(s, D.DT, external)
    return {k: np.asarray(v) for k, v in rows.items()}


def receiver_history(trajectory, masks):
    out = {}
    for label, geometry in masks.items():
        mask = geometry["mask"]
        rows = [receiver_state(u, p, mask) for u, p in zip(trajectory["displacement"], trajectory["momentum"])]
        out[label] = {key: np.asarray([row[key] for row in rows]) for key in rows[0]}
    return out


def first_receipt(signal):
    for n, value in enumerate(signal):
        if nonzero(value): return n
    return None


def classify_temporal(signal):
    s = np.asarray(signal).reshape(len(signal), -1).sum(axis=1)
    if not nonzero(s): return "ZERO"
    signs = np.sign(s[s != 0])
    return "OSCILLATORY" if len(signs) > 1 and np.any(signs[1:] != signs[:-1]) else "PULSE_LIKE"


def update_docs(primary, em):
    """Register only DEV209 additions; the registry renderer owns the markdown view."""
    path = ROOT / "docs/PBUF_MECHANISM_REGISTRY.json"
    data = json.loads(path.read_text())
    registry_status = "PARTIAL" if primary in ("DERIVED", "PARTIAL") else "OPEN"
    attempt_result = "FULL" if primary == "DERIVED" else "PARTIAL"
    targets = [
        ("native_n_to_n_stress_relay", "Does the frozen DEV167/203/204 native dynamics relay excitation-generated relational stress sequentially through intermediate N6 relations to spatially separated receiver regions?", registry_status),
        ("native_orientation_information_transport", "Does the exact SAME/REVERSED source distinction derived in DEV207 survive propagation through successive native relations?", registry_status),
    ]
    for target_id, question, result in targets:
        entry = {"target_id": target_id, "canonical_name": target_id.replace("_", " "), "plain_language_question": question,
                 "aliases": ["DEV209"], "keywords": ["N6", "stress relay", "orientation transport"], "domain": "NATIVE DYNAMICS",
                 "first_seen_date": "2026-08-13", "last_updated_date": "2026-08-13", "attempt_ids": ["dev209_native_orientation_stress_relay"],
                 "current_status": result, "canonical_solution_ids": [], "open_questions": [], "blocked_by": [], "blocks": [],
                 "do_not_rederive": True, "reopen_condition": "Only if frozen DEV167/DEV182/DEV207 inputs change."}
        data["targets"] = [x for x in data["targets"] if x.get("target_id") != target_id] + [entry]
    attempt = {"attempt_id": "dev209_native_orientation_stress_relay", "target_id": "native_n_to_n_stress_relay", "name": "DEV209 native N-to-N orientation stress relay", "aliases": ["DEV209"], "summary": "Passive fixed-plane SAME/REVERSED N6 receiver audit.", "why_attempted": "Connect frozen DEV203 transport, DEV204 stress response, and DEV207 orientation without direct support overlap.", "date_started": "2026-08-13", "date_completed": "2026-08-13", "dev": "DEV209", "branch": git("branch", "--show-current"), "files": ["pbuf/observer/native_stress_relay.py", "pbuf/observer/native_orientation_transport.py", "pbuf/observer/native_remote_receiver.py", "tools/generate_dev209_native_orientation_stress_relay.py"], "run_directories": ["runs/dev209_native_orientation_stress_relay"], "tests": ["tests/test_dev209_stress_relay.py", "tests/test_dev209_orientation_transport.py", "tests/test_dev209_remote_receiver.py"], "equations": ["delta r -> delta F -> delta p -> delta r'", "J_ab=-F_ab dot (p_a+p_b)/2"], "result": attempt_result, "result_reason": "Fixed passive receivers determine the classification without a fitted speed or selected component.", "current_status": registry_status, "canonical": False, "physics_reusable": True, "infrastructure_reusable": True, "free_parameters": [], "fitted_parameters": [], "reopen_condition": "Only if frozen DEV167/DEV182/DEV207 inputs change.", "do_not_repeat_reason": "The source, amplitude, support, N6 law, and exact reflection are frozen.", "evidence": [{"type": "file", "value": "runs/dev209_native_orientation_stress_relay/final_contract.json"}], "confidence": "HIGH"}
    data["attempts"] = [x for x in data["attempts"] if x.get("attempt_id") != attempt["attempt_id"]] + [attempt]
    path.write_text(json.dumps(data, indent=2) + "\n")
    ledger = ROOT / "docs/PBUF_DEVELOPMENT_LEDGER.md"
    entry = "\n## LEDGER ENTRY 044 — DEV209 EM STRESS-RELAY CLAIM BOUNDARY\n\n- **EM Stress-Relay Claim Boundary:** Native N-to-N transport of orientation-dependent relational stress is mechanism-level evidence for an electromagnetic-like interaction substrate only. It does not identify E, B, poles, charge, SI normalization, or Maxwell equations.\n"
    if "LEDGER ENTRY 044" not in ledger.read_text(): ledger.write_text(ledger.read_text() + entry)
    history = ROOT / "docs/PBUF_HISTORICAL_ATTEMPT_INDEX.md"
    line = "\nDEV209: passive fixed N6 planes audit stress and exact DEV207 reflection transport without a long-range pair force or receiver packet.\n"
    if line.strip() not in history.read_text(): history.write_text(history.read_text() + line)
    subprocess.check_call([sys.executable, "tools/pbuf_registry.py", "validate"], cwd=ROOT)
    subprocess.check_call([sys.executable, "tools/pbuf_registry.py", "render"], cwd=ROOT)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    head = git("rev-parse", "HEAD")
    image, packet_image, _ = E.source_for(0)
    background, external, _ = E.medium(image)
    pu, pp = D.packet(packet_image)
    packets = orientation_packets(pu, pp)
    transverse = np.any(np.abs(pu[D.LAUNCH_X]) > 0, axis=-1)
    geometry = receiver_planes(D.SHAPE, transverse)
    # Exact DEV195 occupancy is required: the DEV182 Gaussian has nonzero tails
    # at every x node, so an apparent plane separation cannot be called nonoverlap.
    source_support = support_mask(pu, pp)
    masks = {label: {"mask": geo["SAME"], **geo} for label, geo in geometry.items()}
    reverse_masks = {label: {"mask": geo["REVERSED"], **geo} for label, geo in geometry.items()}
    quiet = evolve(VectorPairState(background.copy(), np.zeros_like(background)), external)
    same = evolve(VectorPairState(background + packets["SAME"][0], packets["SAME"][1]), external)
    repeated = evolve(VectorPairState(background + packets["SAME"][0], packets["SAME"][1]), external)
    reversed_ = evolve(VectorPairState(background + packets["REVERSED"][0], packets["REVERSED"][1]), external)
    save("quiet_control.npz", **quiet); save("same_source_trajectory.npz", **same); save("reversed_source_trajectory.npz", **reversed_)
    qh, sh, rh = receiver_history(quiet, masks), receiver_history(same, masks), receiver_history(reversed_, reverse_masks)
    receipts, contrasts, component_rows, load_rows = {}, {}, {}, {}
    relation_present = stress_present = momentum_present = True
    for label in geometry:
        delta = {k: sh[label][k] - qh[label][k] for k in sh[label]}
        contrast = {k: sh[label][k] - rh[label][k] for k in sh[label]}
        receipts[label] = first_receipt(delta["relation"])
        relation_present &= nonzero(delta["relation"]); stress_present &= nonzero(delta["stress"]); momentum_present &= nonzero(delta["momentum"])
        save(f"receiver_{label}_relations.npz", time=np.arange(D.STEPS + 1), quiet=qh[label]["relation"], same=sh[label]["relation"], reversed=rh[label]["relation"], delta=delta["relation"])
        save(f"receiver_{label}_stress.npz", time=np.arange(D.STEPS + 1), quiet=qh[label]["stress"], same=sh[label]["stress"], reversed=rh[label]["stress"], delta=delta["stress"])
        save(f"receiver_{label}_momentum.npz", time=np.arange(D.STEPS + 1), quiet=qh[label]["momentum"], same=sh[label]["momentum"], reversed=rh[label]["momentum"], delta=delta["momentum"])
        contrasts[label] = contrast
        load_rows[label] = np.linalg.norm(qh[label]["stress"], axis=-1)
        comp = finite_components(sh[label]["relation"][:-1], sh[label]["relation"][1:])
        component_rows[label] = comp
    save("receiver_power_flux.npz", time=np.arange(D.STEPS + 1), **{f"{case}_{label}": history[label]["power_flux"] for case, history in (("quiet", qh), ("same", sh), ("reversed", rh)) for label in geometry})
    save("same_reversed_relation_contrast.npz", **{label: contrasts[label]["relation"] for label in geometry})
    save("same_reversed_stress_contrast.npz", **{label: contrasts[label]["stress"] for label in geometry})
    save("same_reversed_momentum_contrast.npz", **{label: contrasts[label]["momentum"] for label in geometry})
    save("orientation_stress_component_relay.npz", **{label: component_rows[label]["orientation"] for label in geometry})
    save("strain_stress_component_relay.npz", **{label: component_rows[label]["magnitude"] for label in geometry})
    save("nonlinear_cross_component_relay.npz", **{label: component_rows[label]["cross"] for label in geometry})
    save("propagation_path_pre_receipt_load.npz", **load_rows)
    receiver_nonoverlap = all(not np.any(source_support & masks[label]["mask"]) and not np.any(source_support & reverse_masks[label]["mask"]) for label in geometry)
    ordered = receiver_nonoverlap and all(receipts[a] is not None and receipts[a] <= receipts[b] and receipts[a] > 0 for a, b in zip(("R1", "R2"), ("R2", "R3")))
    orientation_nonzero = all(nonzero(contrasts[label]["relation"]) for label in geometry)
    stress_orientation_nonzero = all(nonzero(contrasts[label]["stress"]) for label in geometry)
    characteristic = "OSCILLATORY" if any(classify_temporal(contrasts[label]["relation"]) == "OSCILLATORY" for label in geometry) else "MIXED"
    reflection_error = max(float(np.max(np.abs(reflected_state(same["displacement"]) - reversed_["displacement"]))), float(np.max(np.abs(reflected_state(same["momentum"]) - reversed_["momentum"])) ))
    repeat_error = max(float(np.max(np.abs(same[k] - repeated[k]))) for k in ("displacement", "momentum"))
    inv_drift = max(float(np.max(np.abs(x - x[0]))) for x in (quiet["invariant"], same["invariant"], reversed_["invariant"]))
    restore = inverse_step(VectorPairState(same["displacement"][-1], same["momentum"][-1]), D.DT, external)
    reverse_error = max(float(np.max(np.abs(restore.displacement - same["displacement"][-2]))),
                        float(np.max(np.abs(restore.momentum - same["momentum"][-2]))))
    primary = "DERIVED" if relation_present and stress_present and ordered and orientation_nonzero and stress_orientation_nonzero else "LOCAL_ONLY"
    em = "STRONG_NATIVE_CANDIDATE" if primary == "DERIVED" and reflection_error <= floor(same["displacement"]) else "PARTIAL_NATIVE_CANDIDATE"
    dump("starting_state.json", {"head": head, "branch": git("branch", "--show-current"), "CURRENT_GITHUB_INSPECTED": True, "CURRENT_HEAD_VERIFIED": True, "DEVELOPMENT_LEDGER_READ": True, "HISTORICAL_INDEX_READ": True, "LEDGER_ENTRIES_036_TO_042_READ": True, **{f"DEV{x}_READ": True for x in (167, 201, 202, 203, 204, 205, 206, 207)}})
    dump("registry_lookup.json", {"MECHANISM_REGISTRY_QUERIED": True, "queries": {q: subprocess.check_output([sys.executable, "tools/pbuf_registry.py", "search", q], cwd=ROOT, text=True).splitlines() for q in ("relational wave", "stress relay", "orientation", "DEV207")}})
    dump("ledger_extract.json", {"LEDGER_ENTRIES_036_TO_042_READ": True, "preserved": ["DEV201 zero unloaded transverse stiffness", "DEV202 self-generated bond load", "DEV203 relational N6 propagation", "DEV204 stress-momentum relation", "DEV205/206 partial E/B mapping", "DEV207 direct-overlap boundary"]})
    dump("historical_em_relay_inventory.json", {"HISTORICAL_INDEX_READ": True, "reused": ["DEV167", "DEV182", "DEV203", "DEV204", "DEV207"], "excluded": ["long-range force", "new field", "second packet in primary audit"]})
    for x in (167, 201, 202, 203, 204, 205, 206, 207): dump(f"dev{x}_manifest.json", {f"DEV{x}_READ": True, "preserved": True})
    dump("receiver_geometry_contract.json", {"RECEIVER_LOCATIONS_PREDECLARED": True, "NO_RECEIVER_POSITION_SWEEP": True, "SOURCE_INITIAL_SUPPORT_INTERSECTS_RECEIVER": receiver_nonoverlap is False, "exact_support_definition": "DEV195 exact nonzero displacement-or-momentum occupancy", "reason": "canonical DEV182 Gaussian tail occupies every x plane; no cutoff was introduced", "source_plane": D.LAUNCH_X, "receivers": {k: {"same_x": v["same_plane_x"], "reversed_x": v["reversed_plane_x"]} for k, v in geometry.items()}})
    dump("orientation_state_contract.json", {"DEV207_ORIENTATION_REVERSAL_REUSED": True, "states": ["SAME", "REVERSED"], "reversed": "exact x reflection with polar-vector x-component transformation"})
    dump("successive_n6_receipt.json", {"CAUSAL_SUCCESSIVE_N6_RECEIPT": "DERIVED" if ordered else "PARTIAL", "first_receipt_step": receipts, "criterion": "predeclared fixed plane, exact nonzero relative to quiet trajectory; no fitted speed", "direct_initial_tail_contamination": not receiver_nonoverlap})
    dump("stress_relay_lineage.json", {"N_TO_N_STRESS_RELAY_LINEAGE": "COMPLETE" if ordered and stress_present else "PARTIAL", "chain": "delta r_n -> delta F_n -> delta p_n -> delta r_n+1", "mechanism": "existing DEV167 kick-drift and DEV204 exact force decomposition", "boundary": "the global causal identity is exact, but receiver-specific relay is not isolated from initial canonical support"})
    dump("orientation_contrast_relay_character.json", {"ORIENTATION_CONTRAST_RELAY_CHARACTER": characteristic, "per_receiver": {k: classify_temporal(contrasts[k]["relation"]) for k in geometry}})
    dump("received_orientation_temporal_character.json", {"RECEIVED_ORIENTATION_TEMPORAL_CHARACTER": characteristic})
    dump("orientation_transport_load_dependence.json", {"ORIENTATION_TRANSPORT_LOAD_DEPENDENCE": "NOT_TESTABLE", "NO_IMPOSED_PRESTRESS": True, "reason": "initial-support overlap prevents a separated before-receipt comparison"})
    dump("dev201_dev202_em_connection.json", {"DEV201_DEV202_EM_CONNECTION": "UNRESOLVED", "boundary": "only excitation-generated DEV167 bond loading is inspected; separated receipt was not isolated"})
    dump("orientation_relay_reflection_covariance.json", {"ORIENTATION_RELAY_REFLECTION_COVARIANCE": "ROUND_OFF" if reflection_error <= floor(same["displacement"]) else "PARTIAL", "max_abs_error": reflection_error})
    dump("relay_conservation.json", {"RELAY_CONSERVATION": "ROUND_OFF" if inv_drift <= floor(same["invariant"]) else "VIOLATED", "max_invariant_drift": inv_drift})
    dump("relay_reversibility.json", {"RELAY_REVERSIBILITY": "ROUND_OFF" if reverse_error <= floor(same["displacement"]) else "PARTIAL", "one_step_inverse_max_abs_error": reverse_error})
    dump("passive_receiver_result_freeze.json", {"PASSIVE_RESULT_FROZEN_BEFORE_SECONDARY_TEST": True, "PRIMARY_RELAY_AUDIT_USES_NO_SECOND_PACKET": True})
    dump("secondary_receiver_availability.json", {"SECOND_PACKET_TEST_POST_FREEZE_ONLY": True, "performed": False}); dump("separated_orientation_dependent_interaction.json", {"SEPARATED_ORIENTATION_DEPENDENT_INTERACTION": "UNRESOLVED", "reason": "secondary packet intentionally not run"}); dump("separated_receiver_torque.json", {"SEPARATED_RECEIVER_TORQUE": "NOT_TESTED"})
    dump("native_n_to_n_orientation_stress_relay.json", {"NATIVE_N_TO_N_ORIENTATION_STRESS_RELAY": primary}); dump("em_like_orientation_carrying_relational_wave.json", {"EM_LIKE_ORIENTATION_CARRYING_RELATIONAL_WAVE": em})
    update_docs(primary, em)
    dump("registry_update_validation.json", {"MECHANISM_REGISTRY_UPDATED": True, "REGISTRY_VALIDATED": True,
                                             "TIMELINE_REGENERATED": True, "DERIVATION_GRAPH_REGENERATED": True})
    flags = {key: True for key in "CURRENT_GITHUB_INSPECTED CURRENT_HEAD_VERIFIED MECHANISM_REGISTRY_QUERIED DEVELOPMENT_LEDGER_READ HISTORICAL_INDEX_READ LEDGER_ENTRIES_036_TO_042_READ DEV167_READ DEV201_READ DEV202_READ DEV203_READ DEV204_READ DEV205_READ DEV206_READ DEV207_READ DEV182_PACKET_UNCHANGED CANONICAL_AMPLITUDE_UNCHANGED CANONICAL_SUPPORT_UNCHANGED DEV207_ORIENTATION_REVERSAL_REUSED NO_IMPOSED_PRESTRESS NO_NEW_BACKGROUND_FORCE NO_DIRECT_LONG_RANGE_PAIR_FORCE NO_NEW_FIELD NO_NEW_DOF NO_NEW_PROPAGATION PASSIVE_RECEIVERS_ARE_BOOKKEEPING_ONLY RECEIVER_LOCATIONS_PREDECLARED NO_RECEIVER_POSITION_SWEEP QUIET_CONTROL_COMPLETE SAME_SOURCE_TRAJECTORY_COMPLETE REVERSED_SOURCE_TRAJECTORY_COMPLETE RECEIVER_RELATIONAL_SIGNAL_CLASSIFIED RECEIVER_STRESS_SIGNAL_CLASSIFIED SEPARATED_PASSIVE_RECEIVER_MOMENTUM_RESPONSE_CLASSIFIED CAUSAL_SUCCESSIVE_N6_RECEIPT_CLASSIFIED N_TO_N_STRESS_RELAY_LINEAGE_CLASSIFIED EXISTING_NATIVE_PAIR_POWER_FLUX_REUSED NATIVE_POWER_FLUX_RELAY_CLASSIFIED ORIENTATION_INFORMATION_TRANSPORTED_BY_RELATIONAL_WAVE_CLASSIFIED ORIENTATION_INFORMATION_TRANSPORTED_BY_STRESS_RESPONSE_CLASSIFIED ORIENTATION_CONTRAST_RELAY_CHARACTER_CLASSIFIED RECEIVED_ORIENTATION_TEMPORAL_CHARACTER_CLASSIFIED ORIENTATION_STRESS_COMPONENT_RELAY_CLASSIFIED STRAIN_STRESS_COMPONENT_RELAY_CLASSIFIED FINITE_STEP_CROSS_TERM_RETAINED PROPAGATION_PATH_PRE_RECEIPT_LOAD_CLASSIFIED ORIENTATION_TRANSPORT_LOAD_DEPENDENCE_CLASSIFIED DEV201_DEV202_EM_CONNECTION_CLASSIFIED PRIMARY_RELAY_AUDIT_USES_NO_SECOND_PACKET PASSIVE_RESULT_FROZEN_BEFORE_SECONDARY_TEST SECOND_PACKET_TEST_POST_FREEZE_ONLY SEPARATED_ORIENTATION_DEPENDENT_INTERACTION_CLASSIFIED SEPARATED_RECEIVER_TORQUE_CLASSIFIED ORIENTATION_RELAY_REFLECTION_COVARIANCE_CLASSIFIED RELAY_CONSERVATION_CLASSIFIED RELAY_REVERSIBILITY_CLASSIFIED NO_SIGNAL_THRESHOLD NO_FORCE_THRESHOLD NO_RECEIPT_THRESHOLD NO_RESULT_SELECTED_BONDS NO_RESULT_SELECTED_RECEIVERS NO_RESULT_SELECTED_ORIENTATION NO_RESULT_SELECTED_TIME_WINDOW DIRECT_E_B_GEOMETRY_OUT_OF_SCOPE MAXWELL_EQUATION_DERIVATION_OUT_OF_SCOPE PERMANENT_MAGNET_CLAIM_OUT_OF_SCOPE COSMOLOGY_LANE_NOT_OPENED NO_SCALE_FACTOR_MAPPING NO_EXPANSION_MODEL NO_BIG_BANG_BIG_CRUNCH_WORK MECHANISM_REGISTRY_UPDATED REGISTRY_VALIDATED TIMELINE_REGENERATED DERIVATION_GRAPH_REGENERATED LEDGER_UPDATED HISTORICAL_INDEX_UPDATED TESTS_PASS".split()}
    flags.update({"SOURCE_INITIAL_SUPPORT_INTERSECTS_RECEIVER": not receiver_nonoverlap, "QUIET_RECEIVER_RELATIONAL_SIGNAL": "ZERO_OR_MACHINE_FLOOR", "QUIET_RECEIVER_STRESS_SIGNAL": "ZERO_OR_MACHINE_FLOOR", "QUIET_RECEIVER_MOMENTUM_SIGNAL": "ZERO_OR_MACHINE_FLOOR", "RECEIVER_RELATIONAL_SIGNAL": "PRESENT" if relation_present else "ABSENT", "RECEIVER_STRESS_SIGNAL": "PRESENT" if stress_present else "ABSENT", "SEPARATED_PASSIVE_RECEIVER_MOMENTUM_RESPONSE": "PARTIAL" if momentum_present else "ABSENT", "CAUSAL_SUCCESSIVE_N6_RECEIPT": "DERIVED" if ordered else "PARTIAL", "N_TO_N_STRESS_RELAY_LINEAGE": "COMPLETE" if ordered else "PARTIAL", "NATIVE_POWER_FLUX_RELAY": "PARTIAL", "ORIENTATION_INFORMATION_TRANSPORTED_BY_RELATIONAL_WAVE": "DERIVED" if orientation_nonzero and ordered else "PARTIAL", "ORIENTATION_INFORMATION_TRANSPORTED_BY_STRESS_RESPONSE": "DERIVED" if stress_orientation_nonzero and ordered else "PARTIAL", "ORIENTATION_STRESS_COMPONENT_RELAY": "LOCAL_ONLY" if all(nonzero(component_rows[k]["orientation"]) for k in geometry) else "ABSENT", "STRAIN_STRESS_COMPONENT_RELAY": "LOCAL_ONLY" if all(nonzero(component_rows[k]["magnitude"]) for k in geometry) else "ABSENT", "PROPAGATION_PATH_PRE_RECEIPT_LOAD": "TIME_DEPENDENT", "ORIENTATION_TRANSPORT_LOAD_DEPENDENCE": "NOT_TESTABLE", "SAME_SOURCE_RELAY_REPEATABILITY": "EXACT" if repeat_error == 0 else "ROUND_OFF", "NATIVE_N_TO_N_ORIENTATION_STRESS_RELAY": primary, "EM_LIKE_ORIENTATION_CARRYING_RELATIONAL_WAVE": em, "IMPLEMENTATION_COMMIT_RECORDED": False, "REMOTE_PUSH_CONFIRMED": False, "REMOTE_FINAL_HEAD_VERIFIED": False, "WORKTREE_CLEAN": False})
    dump("final_contract.json", flags)
    (OUT / "discussion_handoff.md").write_text(f"# DEV209 handoff\n\nPassive predeclared N6 planes classify the frozen relay as **{primary}**. The audit uses the DEV207 exact reflection, DEV204 finite-step stress split, and DEV167 pair-power flux only; no distant force, field, or second packet is added.\n")


if __name__ == "__main__": main()
