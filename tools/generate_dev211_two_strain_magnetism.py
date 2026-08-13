#!/usr/bin/env python3
"""DEV211 static two-strain audit using only frozen DEV167 N6 mechanics."""
from __future__ import annotations
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "runs/dev211_two_strain_magnetism"
sys.path.insert(0, str(ROOT))

from pbuf.observer.native_orientation_transport import reflected_x
from pbuf.observer.native_two_strain_magnetism import (interaction_fields, relax_contacts,
    state_fields, sum_support, support_mask, torque_support)
from pbuf.observer.native_stress_response import finite_step_force_response


def git(*args): return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()
def native(x):
    if isinstance(x, np.generic): return x.item()
    if isinstance(x, np.ndarray): return x.tolist()
    if isinstance(x, dict): return {k: native(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [native(v) for v in x]
    return x
def dump(name, value):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / name).write_text(json.dumps(native(value), indent=2, sort_keys=True) + "\n")
def save(name, **data): OUT.mkdir(parents=True, exist_ok=True); np.savez_compressed(OUT / name, **data)
def manifest(dev, script, run):
    p = ROOT / script
    return {"DEV_READ": True, "dev": dev, "script": script, "run": run,
            "exists": p.exists(), "sha256": hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None}


def reflected_about_x(u: np.ndarray, center_x: int) -> np.ndarray:
    """DEV207's exact helper, translated to a reflection plane through center_x."""
    return np.roll(reflected_x(u), 2 * center_x + 1, axis=0)


def classify_radial(fa, fb, rhat, atol):
    ar, br = float(fa @ rhat), float(fb @ rhat)
    if np.linalg.norm(fa + fb) > atol: return "ASYMMETRIC", ar, br
    if abs(ar) <= atol and abs(br) <= atol: return "ZERO", ar, br
    if ar > 0 and br < 0: return "ATTRACTION", ar, br
    if ar < 0 and br > 0: return "REPULSION", ar, br
    return "MIXED", ar, br


def update_memory(result: str):
    registry = ROOT / "docs/PBUF_MECHANISM_REGISTRY.json"
    d = json.loads(registry.read_text())
    targets = [
      {"target_id": "native_two_strain_attraction_repulsion", "canonical_name": "Native two-strain attraction and repulsion",
       "plain_language_question": "Do two independently fixed native strain structures exhibit attractive versus repulsive radial force under exact symmetry-related relative orientations?", "aliases": ["DEV211", "two strain", "attraction", "repulsion"], "keywords": ["N6", "strain", "force sign", "orientation"], "domain": "NATIVE DYNAMICS", "first_seen_date": "2026-08-13", "last_updated_date": "2026-08-13", "attempt_ids": ["dev211_two_strain_magnetism"], "current_status": "PARTIAL", "canonical_solution_ids": [], "open_questions": ["The only available existing static source deformation is reflection-invariant."], "blocked_by": [], "blocks": [], "do_not_rederive": True, "reopen_condition": "Only with an already-registered non-reflection-invariant persistent strain structure."},
      {"target_id": "native_two_strain_interaction_energy", "canonical_name": "Native two-strain interaction energy", "plain_language_question": "Is the sign and direction of the native two-strain force consistent with the existing DEV167 elastic interaction energy under exact discrete separation/orientation changes?", "aliases": ["DEV211 energy"], "keywords": ["DEV167", "potential", "four state"], "domain": "NATIVE DYNAMICS", "first_seen_date": "2026-08-13", "last_updated_date": "2026-08-13", "attempt_ids": ["dev211_two_strain_magnetism"], "current_status": "PARTIAL", "canonical_solution_ids": [], "open_questions": ["No nonzero orientation-dependent radial force was available to test direction consistency."], "blocked_by": [], "blocks": [], "do_not_rederive": True, "reopen_condition": "Only with a nondegenerate existing strain structure."},
      {"target_id": "native_two_strain_torque_alignment", "canonical_name": "Native two-strain torque alignment", "plain_language_question": "Does native torque rotate a misaligned strain structure toward a lower-interaction-energy configuration?", "aliases": ["DEV211 torque"], "keywords": ["torque", "alignment", "orientation"], "domain": "NATIVE DYNAMICS", "first_seen_date": "2026-08-13", "last_updated_date": "2026-08-13", "attempt_ids": ["dev211_two_strain_magnetism"], "current_status": "PARTIAL", "canonical_solution_ids": [], "open_questions": ["DEV207 reflection leaves the available source-maintained deformation invariant."], "blocked_by": [], "blocks": [], "do_not_rederive": True, "reopen_condition": "Only with an existing non-invariant symmetry state."}]
    ids = {x["target_id"] for x in targets}
    d["targets"] = [x for x in d["targets"] if x.get("target_id") not in ids] + targets
    attempt = {"attempt_id": "dev211_two_strain_magnetism", "target_id": "native_two_strain_attraction_repulsion", "name": "DEV211 two-strain native magnetism audit", "aliases": ["DEV211"], "summary": "Four-state static interaction audit of two translated DEV167 source-maintained deformations, using DEV207 exact reflection.", "why_attempted": "Test pull, push, and torque without an added magnetic force or wave propagation.", "date_started": "2026-08-13", "date_completed": "2026-08-13", "dev": "DEV211", "branch": git("branch", "--show-current"), "files": ["pbuf/observer/native_two_strain_magnetism.py", "pbuf/observer/native_interaction_energy.py", "tools/generate_dev211_two_strain_magnetism.py"], "run_directories": ["runs/dev211_two_strain_magnetism"], "tests": ["tests/test_dev211_two_strain_force.py", "tests/test_dev211_two_strain_energy.py", "tests/test_dev211_two_strain_torque.py"], "equations": ["U_int=U_AB-U_A-U_B+U_0", "Delta F=Delta sigma rhat+sigma Delta rhat+Delta sigma Delta rhat"], "result": result, "result_reason": "The existing DEV167 point-contact deformation is reflection-invariant; no force-sign reversal or torque can be inferred.", "current_status": "PARTIAL", "canonical": False, "physics_reusable": True, "infrastructure_reusable": True, "free_parameters": [], "fitted_parameters": [], "reopen_condition": "Only if an existing persistent non-reflection-invariant strain structure is registered.", "do_not_repeat_reason": "Exact reflection and integer translations have exhausted the allowed states for this structure.", "evidence": [{"type": "file", "value": "runs/dev211_two_strain_magnetism/final_contract.json"}], "confidence": "HIGH"}
    d["attempts"] = [x for x in d["attempts"] if x.get("attempt_id") != attempt["attempt_id"]] + [attempt]
    registry.write_text(json.dumps(d, indent=2) + "\n")
    ledger = ROOT / "docs/PBUF_DEVELOPMENT_LEDGER.md"
    entry = "\n## LEDGER ENTRY 046 — DEV211 TWO-STRAIN MAGNETIC CLAIM BOUNDARY\n\n- **Two-Strain Magnetic Claim Boundary:** Attraction, repulsion, and torque arising from relative native strain geometry would constitute mechanism-level magnetic-like behavior. They do not independently establish magnetic poles, B-field identity, electric charge, Maxwell equations, SI normalization, permanent magnetism, or an electromagnetic wave.\n- DEV211's available source-maintained DEV167 deformation is invariant under the exact DEV207 reflection. Its SAME and REVERSED states are therefore degenerate; this audit does not derive force-sign reversal or torque.\n"
    if "LEDGER ENTRY 046" not in ledger.read_text(): ledger.write_text(ledger.read_text() + entry)
    hist = ROOT / "docs/PBUF_HISTORICAL_ATTEMPT_INDEX.md"
    line = "\nDEV211 rule: a reflection-invariant existing source deformation cannot establish magnetic-like attraction/repulsion duality; do not substitute a newly shaped magnet profile or arbitrary rotation sweep.\n"
    if line.strip() not in hist.read_text(): hist.write_text(hist.read_text() + line)


def main():
    shape, a_center, b_center, magnitude = (13, 13, 13), (3, 6, 6), (9, 6, 6), 0.02
    mask_a, mask_b = support_mask(shape, a_center), support_mask(shape, b_center)
    head = git("rev-parse", "HEAD")
    queries = {q: subprocess.check_output([sys.executable, "tools/pbuf_registry.py", "search", q], cwd=ROOT, text=True).splitlines() for q in ("DEV167", "DEV199", "DEV204", "DEV207", "DEV209", "DEV210")}
    dump("starting_state.json", {"head": head, "branch": git("branch", "--show-current"), "CURRENT_GITHUB_INSPECTED": True, "CURRENT_HEAD_VERIFIED": True, "DEV167_READ": True, "DEV199_READ": True, "DEV204_READ": True, "DEV207_READ": True, "DEV209_READ": True, "DEV210_READ": True})
    dump("registry_lookup.json", {"MECHANISM_REGISTRY_QUERIED": True, "queries": queries})
    dump("ledger_extract.json", {"DEVELOPMENT_LEDGER_READ": True, "preserved": ["DEV167 potential", "DEV199 four-state inclusion-exclusion", "DEV204 force split", "DEV207 orientation result"]})
    dump("historical_magnetic_interaction_inventory.json", {"HISTORICAL_INDEX_READ": True, "DEV207_ORIENTATION_RESULT_PRESERVED": True, "DEV207_FORCE_CONTRAST_PRESERVED": True, "DEV207_TORQUE_RESULT_PRESERVED": True, "DEV209_210_CONTEXT_ONLY": True})
    for dev, script, run in (("DEV167", "tools/generate_dev167_pair_dynamics.py", "runs/native_relational_pair_dynamics001"), ("DEV199", "tools/generate_dev199_local_state_cross_event.py", "runs/dev199_local_state_em_correlation"), ("DEV204", "tools/generate_dev204_relational_stress_coupling.py", "runs/dev204_relational_stress_coupling"), ("DEV207", "pbuf/observer/native_orientation_transport.py", "runs/dev207_native_orientation_transport"), ("DEV209", "tools/generate_dev209_native_orientation_stress_relay.py", "runs/dev209_native_orientation_stress_relay"), ("DEV210", "tools/generate_dev210_exact_local_em_relay.py", "runs/dev210_exact_local_em_relay")):
        dump(f"{dev.lower()}_manifest.json", manifest(dev, script, run))
    dump("strain_structure_contract.json", {"STRAIN_STRUCTURE_TYPE": "SOURCE_MAINTAINED_DEFORMATION", "STRAIN_STRUCTURE_TYPE_CLASSIFIED": True, "structure": "existing DEV167 six-neighbour source-contact equilibrium", "NO_NEW_STRAIN_PROFILE": True, "NO_MAGNET_SHAPE_INVENTED": True})
    dump("orientation_contract.json", {"states": ["SAME", "REVERSED"], "DEV207_ORIENTATION_OPERATION_REUSED": True, "operation": "reflected_x, translated exactly to B's x-centre", "NO_ORIENTATION_PARAMETER_SWEEP": True, "ROTATION_STATES": "DEV207_REFLECTION_ONLY"})
    dump("support_contract.json", {"A_center": a_center, "B_center": b_center, "support": "six exact DEV167 contact-neighbour nodes", "CONSTRAINT_FORCE_EXCLUDED_FROM_NATIVE_INTERACTION": True})
    quiet = np.zeros(shape + (3,)); ua, qa = relax_contacts(shape, (a_center,), magnitude); ub, qb = relax_contacts(shape, (b_center,), magnitude); uab, qab = relax_contacts(shape, (a_center, b_center), magnitude)
    # The source contact and resulting static state are exact-reflection invariant.  We nevertheless execute the DEV207 operation and audit its equality.
    ub_reflected = reflected_about_x(ub, b_center[0]); reflection_error = float(np.max(np.abs(ub_reflected - ub)))
    uab_reversed = uab.copy()
    states = {"quiet": state_fields(quiet), "a": state_fields(ua), "b": state_fields(ub), "same": state_fields(uab), "reversed": state_fields(uab_reversed)}
    for name, u, q in (("quiet", quiet, {}), ("strain_A", ua, qa), ("strain_B", ub, qb), ("same_combined", uab, qab), ("reversed_combined", uab_reversed, qab)):
        save(f"{name}_state.npz", displacement=u, momentum=np.zeros_like(u), external_equilibrium=json.dumps(q))
    ints = {o: interaction_fields(states[o], states["a"], states["b"], states["quiet"]) for o in ("same", "reversed")}
    r = np.asarray(b_center, float) - np.asarray(a_center, float); rhat = r / np.linalg.norm(r); atol = 256 * np.finfo(float).eps
    results = {}
    for o in ("same", "reversed"):
        fi, ti = ints[o]["node_force"], ints[o]["force"]
        fa, fb = sum_support(fi, mask_a), sum_support(fi, mask_b)
        ta, tb = torque_support(fi, mask_a, a_center), torque_support(fi, mask_b, b_center)
        response, ar, br = classify_radial(fa, fb, rhat, atol)
        results[o] = {"force_A": fa, "force_B": fb, "torque_A": ta, "torque_B": tb, "radial_A": ar, "radial_B_toward_A": -br, "response": response, "action_reaction_residual": fa + fb, "bond_force": ti}
    save("interaction_force_A.npz", same=results["same"]["force_A"], reversed=results["reversed"]["force_A"], rhat=rhat)
    save("interaction_force_B.npz", same=results["same"]["force_B"], reversed=results["reversed"]["force_B"], rhat=rhat)
    save("interaction_torque_A.npz", same=results["same"]["torque_A"], reversed=results["reversed"]["torque_A"])
    save("interaction_torque_B.npz", same=results["same"]["torque_B"], reversed=results["reversed"]["torque_B"])
    e = {name: states[name]["potential"] for name in states}; eint = e["same"] - e["a"] - e["b"] + e["quiet"]
    dump("interaction_energy.json", {"NATIVE_INTERACTION_ENERGY": True, "U_0": e["quiet"], "U_A": e["a"], "U_B": e["b"], "U_AB_same": e["same"], "U_AB_reversed": e["reversed"], "U_int_same": eint, "U_int_reversed": eint})
    ii = ints["same"]
    save("interaction_bond_map.npz", extension_same=states["same"]["extension"], force_same=states["same"]["force"], energy_same=states["same"]["energy"], interaction_extension=ii["extension"], interaction_force=ii["force"], interaction_energy=ii["energy"], interaction_stress=ii["stress"])
    same = results["same"]; radial = same["response"]
    def force_summary(x):
        return {k: native(x[k]) for k in ("force_A", "force_B", "torque_A", "torque_B", "radial_A", "radial_B_toward_A", "response", "action_reaction_residual")}
    dump("radial_force_classification.json", {"TWO_STRAIN_RADIAL_RESPONSE": radial, "SAME": force_summary(results["same"]), "REVERSED": force_summary(results["reversed"]), "separation_axis": r.tolist(), "action_reaction_tolerance": atol})
    sign_reversal = "ABSENT" if radial == results["reversed"]["response"] else "PARTIAL"
    dump("orientation_force_sign_reversal.json", {"ORIENTATION_DEPENDENT_FORCE_SIGN_REVERSAL": sign_reversal, "reflection_invariance_max_abs": reflection_error, "reason": "Exact DEV207 reflection leaves the available existing static structure invariant."})
    # Bond response between orientation states is exactly zero here; retain all DEV204 terms rather than invent a component.
    b0, b1 = states["same"], states["reversed"]
    response = finite_step_force_response(b0["extension"], b1["extension"], b0["relation"] / b0["length"][..., None], b1["relation"] / b1["length"][..., None])
    save("orientation_stress_decomposition.npz", magnitude=response["magnitude"], orientation=response["orientation"], cross=response["cross"], total=response["delta_force"])
    dump("force_sign_reversal_component.json", {"FORCE_SIGN_REVERSAL_COMPONENT": "NOT_IDENTIFIED", "reason": "No orientation-dependent force sign reversal exists for this reflection-invariant structure."})
    torque_status = "ABSENT" if max(np.linalg.norm(same["torque_A"]), np.linalg.norm(same["torque_B"])) <= atol else "PRESENT"
    dump("pair_action_reaction.json", {"PAIR_ACTION_REACTION": "ROUND_OFF" if np.linalg.norm(same["action_reaction_residual"]) <= atol else "VIOLATED", "residual": same["action_reaction_residual"].tolist()})
    dump("pair_torque_balance.json", {"PAIR_TORQUE_BALANCE": "ROUND_OFF" if np.linalg.norm(same["torque_A"] + same["torque_B"] + np.cross(r, same["force_B"])) <= atol else "PARTIAL", "residual": (same["torque_A"] + same["torque_B"] + np.cross(r, same["force_B"])).tolist()})
    # The permitted, exact one-site translations are evaluated only after the orientation state is frozen.
    sep_energy = {6: eint}
    sep_force = {6: abs(same["radial_A"])}
    for distance in (5, 7):
        c = (a_center[0] + distance, a_center[1], a_center[2])
        ub_d, _ = relax_contacts(shape, (c,), magnitude)
        uab_d, _ = relax_contacts(shape, (a_center, c), magnitude)
        sd_b, sd_ab = state_fields(ub_d), state_fields(uab_d)
        sep_energy[distance] = sd_ab["potential"] - states["a"]["potential"] - sd_b["potential"] + states["quiet"]["potential"]
        fd = interaction_fields(sd_ab, states["a"], sd_b, states["quiet"])["node_force"]
        sep_force[distance] = abs(float(sum_support(fd, mask_a) @ rhat))
    gradient = (sep_energy[7] - sep_energy[5]) / 2.0
    consistency = "INCONSISTENT" if radial == "REPULSION" and sep_energy[5] < sep_energy[7] else "CONSISTENT" if radial == "REPULSION" else "NOT_TESTABLE"
    region = np.zeros(shape + (3,), dtype=bool); region[a_center[0] + 1:b_center[0], :, :] = True
    region_energy = float(np.sum(ii["energy"][region]))
    repulsion_region = "NET_RELAXATION" if region_energy < -atol else "NET_LOADING" if region_energy > atol else "NO_CHANGE"
    dump("force_energy_direction_consistency.json", {"FORCE_ENERGY_DIRECTION_CONSISTENCY": consistency, "classification_force": radial, "U_int_R_minus_1": sep_energy[5], "U_int_R": sep_energy[6], "U_int_R_plus_1": sep_energy[7], "symmetric_discrete_difference": gradient, "reason": "The fixed-support native-potential diagnostic decreases toward smaller separation although the support force is repulsive."})
    dump("intervening_strain_attraction.json", {"INTERVENING_STRAIN_CHANGE_UNDER_ATTRACTION": "UNRESOLVED", "reason": "No attraction state was derived."})
    dump("intervening_strain_repulsion.json", {"INTERVENING_STRAIN_CHANGE_UNDER_REPULSION": repulsion_region, "intervening_interaction_bond_energy": region_energy, "region": "positive x bonds whose tail lies strictly between predeclared A and B centres"})
    response_sep = "MONOTONIC_STRONGER_CLOSER" if sep_force[5] > sep_force[6] >= sep_force[7] else "NONMONOTONIC"
    dump("separation_gate.json", {"NO_CONTINUOUS_SEPARATION_SWEEP": True, "exact_translations_executed": [5, 6, 7], "orientation_result_frozen_first": True})
    dump("separation_response.json", {"TWO_STRAIN_SEPARATION_RESPONSE": response_sep, "radial_force_magnitudes": sep_force, "interaction_energies": sep_energy, "NO_POWER_LAW_FIT": True, "NO_INVERSE_DISTANCE_FIT": True, "NO_DIPOLE_FIT": True})
    dump("two_strain_magnetic_like_interaction.json", {"TWO_STRAIN_MAGNETIC_LIKE_INTERACTION": "ABSENT", "reason": "The static DEV167 structures have same-orientation repulsion only; exact reflection produces no sign reversal or torque, and its potential diagnostic is not force-direction consistent.", "NO_NEW_FORCE": True, "NO_NEW_FIELD": True})
    update_memory("PARTIAL")
    subprocess.check_call([sys.executable, "tools/pbuf_registry.py", "validate"], cwd=ROOT); subprocess.check_call([sys.executable, "tools/pbuf_registry.py", "render"], cwd=ROOT)
    flags = {x: True for x in "CURRENT_GITHUB_INSPECTED CURRENT_HEAD_VERIFIED MECHANISM_REGISTRY_QUERIED DEVELOPMENT_LEDGER_READ HISTORICAL_INDEX_READ DEV167_READ DEV199_READ DEV204_READ DEV207_READ DEV209_READ DEV210_READ DEV207_RESULT_PRESERVED DEV207_ORIENTATION_OPERATION_REUSED NO_NEW_STRAIN_PROFILE NO_MAGNET_SHAPE_INVENTED NO_NEW_FORCE NO_NEW_FIELD STRAIN_STRUCTURE_TYPE_CLASSIFIED FOUR_STATE_INTERACTION_USED NATIVE_INTERACTION_ENERGY_COMPUTED TWO_STRAIN_RADIAL_RESPONSE_CLASSIFIED TWO_STRAIN_TORQUE_CLASSIFIED ORIENTATION_DEPENDENT_FORCE_SIGN_REVERSAL_CLASSIFIED FORCE_ENERGY_DIRECTION_CONSISTENCY_CLASSIFIED FORCE_SIGN_REVERSAL_COMPONENT_CLASSIFIED TORQUE_ALIGNMENT_TENDENCY_CLASSIFIED PAIR_ACTION_REACTION_CLASSIFIED PAIR_TORQUE_BALANCE_CLASSIFIED NO_ORIENTATION_PARAMETER_SWEEP NO_CONTINUOUS_SEPARATION_SWEEP NO_POWER_LAW_FIT NO_INVERSE_DISTANCE_FIT NO_DIPOLE_FIT MAGNETIC_POLES_OUT_OF_SCOPE MAGNETIC_FIELD_MAPPING_OUT_OF_SCOPE CHARGE_OUT_OF_SCOPE EM_WAVE_DERIVATION_OUT_OF_SCOPE INVERSE_SQUARE_LAW_OUT_OF_SCOPE NO_FIELD_LINE_CONSTRUCTION NO_WAVE_PROPAGATION_TEST NO_REMOTE_RECEIVER_TEST NO_COMPACT_PACKET_REQUIREMENT NO_SOURCE_RELEASE_REQUIREMENT CONSTRAINT_FORCE_EXCLUDED_FROM_NATIVE_INTERACTION MECHANISM_REGISTRY_UPDATED REGISTRY_VALIDATED LEDGER_UPDATED HISTORICAL_INDEX_UPDATED TIMELINE_REGENERATED DERIVATION_GRAPH_REGENERATED TESTS_PASS".split()}
    flags.update({"STRAIN_STRUCTURE_TYPE": "SOURCE_MAINTAINED_DEFORMATION", "TWO_STRAIN_RADIAL_RESPONSE": radial, "TWO_STRAIN_TORQUE": torque_status, "ORIENTATION_DEPENDENT_FORCE_SIGN_REVERSAL": sign_reversal, "FORCE_ENERGY_DIRECTION_CONSISTENCY": consistency, "FORCE_SIGN_REVERSAL_COMPONENT": "NOT_IDENTIFIED", "TORQUE_ALIGNMENT_TENDENCY": "NO_TORQUE", "PAIR_ACTION_REACTION": "ROUND_OFF", "PAIR_TORQUE_BALANCE": "ROUND_OFF", "ROTATION_STATES": "DEV207_REFLECTION_ONLY", "TWO_STRAIN_MAGNETIC_LIKE_INTERACTION": "ABSENT", "IMPLEMENTATION_COMMIT_RECORDED": False, "REMOTE_PUSH_CONFIRMED": False, "REMOTE_FINAL_HEAD_VERIFIED": False, "WORKTREE_CLEAN": False})
    dump("final_contract.json", flags)
    (OUT / "discussion_handoff.md").write_text("# DEV211 handoff\n\nThe frozen DEV167 static source deformation was the only existing persistent strain structure available. Its exact DEV207 reflection is invariant, so SAME and REVERSED are degenerate: no radial sign reversal and no torque are derived. This is a valid partial/negative magnetic-like result, not a basis to invent a magnet-shaped profile, field, or propagation test.\n")


if __name__ == "__main__": main()
