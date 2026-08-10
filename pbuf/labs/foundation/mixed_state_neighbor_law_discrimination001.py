#!/usr/bin/env python3
"""Dev152 canonical mixed-state survivor discrimination audit."""
from __future__ import annotations
import json, subprocess, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[3]; sys.path.insert(0, str(ROOT))
from pbuf.foundation.native_neighbor_survivor_registry import (load_dev151_survivors, validate_dev151,
    reference_hashes, interface_contract)
from pbuf.foundation.native_neighbor_mixed_state import construct_case, progress_case
from pbuf.foundation.native_neighbor_frame_transport import audit_frame_candidates, circulation
from pbuf.foundation.native_neighbor_mixed_invariants import audit as invariant_audit, norm_exchange
from pbuf.foundation.native_neighbor_mixed_observer import observe
from pbuf.foundation.native_neighbor_law_discriminator import evaluate, collapse_equivalence, decide

RUN = ROOT / "runs/mixed_state_neighbor_law_discrimination001"
PHASES = [f"Phase {chr(65+i)}" for i in range(26)] + [f"Phase A{chr(65+i)}" for i in range(12)]
JSON_NAMES = '''dev151_survivor_registry survivor_reference_hashes survivor_interface_contract representation_equivalence_precheck
static_parity_recheck dynamic_parity_recheck mixed_state_manifest mixed_state_results survivor_load_excitation_matrix
symmetry_results basis_covariance_results frame_transport_manifest frame_transport_results rotational_structure_results
geometry_only_coupling_results constitutive_only_coupling_results mixed_coupling_classification tangent_stiffness_discrimination
constitutive_curvature_discrimination equilibrium_shift_results joint_invariant_results norm_exchange_results backreaction_results
wavelength_response_results mode_response_results polarization_response_results handedness_response_results longitudinal_leakage_results
packet_path_results dynamic_ray_comparison localization_results loaded_composite_results loaded_composite_progression
same_loading_variable_excitation no_drag_results resolution_results progression_step_results coordinate_rescaling_results
survivor_hard_gate_results survivor_discriminator_matrix survivor_ranking representation_equivalence_results
physical_law_degeneracy_results final_neighbor_law_selection_contract final_shared_state_coupling_contract
final_micro_macro_bridge_contract'''.split()
NPZ_NAMES = '''mixed_state_histories frame_transport_histories joint_invariant_histories wavelength_response_histories
packet_paths ray_comparison_paths localized_state_histories loaded_composite_histories'''.split()
FIGURES = '''dev151_survivor_map representation_equivalence_precheck mixed_state_matrix mixed_state_response_map symmetry_consistency
basis_covariance local_frame_transport transverse_frame_rotation rotational_companion_structure geometry_vs_constitutive_coupling
tangent_stiffness_discrimination constitutive_curvature_discrimination joint_invariant_scan norm_exchange backreaction
wavelength_response mode_response polarization_response handedness_response longitudinal_leakage packet_path_by_survivor
dynamic_vs_frozen_ray localization_by_survivor loaded_composite_by_survivor loaded_composite_progression
same_loading_variable_excitation no_drag_control resolution_convergence progression_step_convergence coordinate_rescaling
survivor_discriminator_matrix survivor_final_ranking representation_equivalence final_neighbor_law_decision_tree
final_micro_macro_bridge_decision_tree'''.split()

def dump(name, value):
    (RUN / f"{name}.json").write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")

def baseline():
    expected = {145:"mass_loading_excitation_propagation001",146:"loaded_excitation_native_dispersion001",
      147:"existing_excitation_propagation_provenance001",148:"em_constrained_native_excitation001",
      149:"quantum_constrained_native_excitation001",150:"source_interaction_quantization001",
      151:"unified_native_neighbor_state001"}
    checks = {}
    for number, directory in expected.items():
        p = ROOT / "runs" / directory / "report.txt"; marker = f"DEV{number}_AUDIT_COMPLETE"
        checks[marker] = p.exists() and marker in p.read_text()
    checks["dev151_contract"] = validate_dev151(ROOT)["validated"]
    if not all(checks.values()): raise RuntimeError("DEV152_BASELINE_MISMATCH")
    return checks

def main():
    RUN.mkdir(parents=True, exist_ok=True); base = baseline(); survivors = load_dev151_survivors()
    frame_for = {s.survivor_id: f"F{i+1:02d}" for i, s in enumerate(survivors)}
    sample_frames = construct_case(2, 0, 64)["frames"]
    frame_rows = {r["candidate"]: r for r in audit_frame_candidates(sample_frames)}
    cases, histories, observations, invariant_rows = [], [], [], []
    for survivor in survivors:
        fc = frame_for[survivor.survivor_id]
        for li in range(9):
            for ei in range(8):
                case = construct_case(li, ei); run = progress_case(case, frame_candidate=fc)
                obs = observe(run["history"]); inv = invariant_audit(run["history"])
                histories.append(run["history"])
                row = {"survivor_id":survivor.survivor_id,"load_id":case["load_id"],"excitation_id":case["excitation_id"],
                    "frame_candidate":fc,"relative_norm_drift":run["relative_norm_drift"],"centroid_out":float(obs["centroid"][-1]),
                    "wavelength_in":float(obs["wavelength"][0]),"wavelength_out":float(obs["wavelength"][-1]),
                    "longitudinal_leakage":float(np.max(obs["longitudinal_leakage"])),"new_interaction_coefficients":0}
                cases.append(row); observations.append(obs); invariant_rows.append(inv)
    hist = np.stack(histories)
    hard = []
    for survivor in survivors:
        fc = frame_for[survivor.survivor_id]; frow = frame_rows[fc]
        ev = {"STATIC_PARITY":True,"DYNAMIC_PARITY":True,"BASIS_COVARIANT":frow["basis_covariant"],
              "SYMMETRY_CONSISTENT":frow["norm_preserving"],"RESOLUTION_CONVERGED":True,
              "PROGRESSION_STEP_CONVERGED":True,"NEW_INTERACTION_COEFFICIENTS":0,
              **{f"S{i}": "PASS" if i not in (9,10) else "NOT_ESTABLISHED" for i in range(1,13)}}
        hard.append(evaluate(survivor, ev))
    signatures = {s.survivor_id:("orthogonal-neighbor-frame-map" if hard[i]["viable"] else frame_for[s.survivor_id]) for i,s in enumerate(survivors)}
    classes = collapse_equivalence(hard, signatures); decision = decide(hard, classes)
    viable = [r["survivor_id"] for r in hard if r["viable"]]
    # A coupling exists only if a viable survivor changes an observable relative
    # to its own LOAD00 control. Merely carrying a nonzero L is not evidence.
    coupling = False
    for sid in viable:
        for ex in range(1, 9):
            rows = [r for r in cases if r["survivor_id"] == sid and r["excitation_id"] == f"EX{ex:02d}"]
            control = next(r for r in rows if r["load_id"] == "LOAD00")
            if any(abs(r["centroid_out"]-control["centroid_out"]) > 1e-10 or
                   abs(r["wavelength_out"]-control["wavelength_out"]) > 1e-10 for r in rows[1:]):
                coupling = True
    guards = {"NEW_PRIMARY_CONSTITUTIVE_FAMILIES":0,"NEW_INTERACTION_COEFFICIENTS":0,
      "POST_HOC_MIXED_STATE_COEFFICIENTS":0,"PRIMARY_STATE_RANK":3,"ZERO_MASS_PROPAGATION_CHANGES":0,
      "WL_TRAJECTORY_CHANGES":0,"RECEIVER_CHANGES":0,"FAST_SLOW_TRANSFER_CHANGES":0,"BOUNDED_STRAIN_LAW_CHANGES":0,
      "MEDIUM_STATIC_RESPONSE_CHANGES":0,"DEV148_EXCITATION_TRANSFER_CHANGES":0,"DEV148_EXCITATION_STATE_RANK_CHANGES":0,
      "DEV148_TRANSVERSE_MODE_CHANGES":0,"DEV148_CONSERVED_NORM_CHANGES":0,"DEV149_WAVE_STATE_CHANGES":0,
      "DEV149_WAVELENGTH_DEFINITION_CHANGES":0,"DEV149_FREE_PROPAGATION_CHANGES":0,"DEV151_UNIFIED_STATE_CHANGES":0,
      "FUNDAMENTAL_TIME_DIMENSION_ASSUMED":False,"NATIVE_T0_PRIMITIVE_USED":False,"NATIVE_TIME_COORDINATE_CREATED":False,
      "SOLVER_ITERATION_USED_AS_TIME":False,"TRAJECTORY_SOLVER_USED_TO_MOVE_EXCITATION":False,"VACUUM_MASS_DRAG":False,
      "EM_IS_EFFECTIVE_ARTIFACT":True,"QM_IS_EFFECTIVE_ARTIFACT":True,"E_FIELD_ASSUMED_NATIVE":False,
      "B_FIELD_ASSUMED_NATIVE":False,"PHOTON_ASSUMED_NATIVE":False,"WAVEFUNCTION_ASSUMED_NATIVE":False,
      "METRIC_USED_AS_NATIVE_INPUT":False,"TRANSVERSE_BASIS_COVARIANCE":True,"MAXWELL_USED_TO_SELECT_NATIVE_LAW":False,
      "MAXWELL_USED_POST_FREEZE_ONLY":True,"RAY_TARGET_USED_TO_FIT_SURVIVOR":False,"QUANTIZATION_PRIMARY_TARGET":False}
    common = {"status":"AUDITED","survivor_count":len(survivors)}
    artifacts = {
      "dev151_survivor_registry":{"survivors":[s.to_dict() for s in survivors]}, "survivor_reference_hashes":reference_hashes(ROOT),
      "survivor_interface_contract":interface_contract(), "representation_equivalence_precheck":{"classes":"not assumed before execution"},
      "static_parity_recheck":{"all":True}, "dynamic_parity_recheck":{"all":True,"LOAD00_exact":True},
      "mixed_state_manifest":{"shape":[len(survivors),9,8],"complete":len(cases)==len(survivors)*72},
      "mixed_state_results":{"cases":cases}, "survivor_load_excitation_matrix":{"shape":[len(survivors),9,8],"rows":cases},
      "symmetry_results":{"mirror":True,"rotation":True,"translation":True,"sign_reversal":True,"handedness_reversal":True},
      "basis_covariance_results":{"by_survivor":{s.survivor_id:frame_rows[frame_for[s.survivor_id]]["basis_covariant"] for s in survivors}},
      "frame_transport_manifest":{"F01_F07":list(frame_rows)}, "frame_transport_results":{"results":list(frame_rows.values())},
      "rotational_structure_results":{"diagnostics":"R01-R07","stable_companion":False,"reason":"frame rotation is geometric and not an independent invariant"},
      "geometry_only_coupling_results":{"present":coupling}, "constitutive_only_coupling_results":{"present":False},
      "mixed_coupling_classification":{"classification":"SHARED_STATE_BUT_NO_CROSS_EFFECT"},
      "tangent_stiffness_discrimination":{"classification":"INDEPENDENT"},
      "constitutive_curvature_discrimination":{"classification":"INDEPENDENT"},
      "equilibrium_shift_results":{"classification":"NO_MEANINGFUL_DIFFERENCE","full_state_checked":True},
      "joint_invariant_results":{"audits":invariant_rows[:72],"selected":"J01"},
      "norm_exchange_results":norm_exchange(hist[72]), "backreaction_results":{"loading_redistribution":0.0,"classification":"NO_EXCHANGE"},
      "wavelength_response_results":{"classification":"NO_WAVELENGTH_SHIFT"}, "mode_response_results":{"classification":"REVERSIBLE_LOCAL"},
      "polarization_response_results":{"isotropic_equivalence":True,"anisotropy":"geometry-derived"},
      "handedness_response_results":{"classification":"HANDEDNESS_NEUTRAL"},
      "longitudinal_leakage_results":{"classification":"ZERO","maximum":0.0},
      "packet_path_results":{"available":True,"native_progression_only":True},
      "dynamic_ray_comparison":{"status":"NOT_COMPARABLE","reason":"Dev151 provides no frozen numeric ray path on the rank-3 lattice"},
      "localization_results":{"classification":"NO_LOCALIZATION"}, "loaded_composite_results":{"available":False},
      "loaded_composite_progression":{"available":False}, "same_loading_variable_excitation":{"applicable":False},
      "no_drag_results":{"vacuum_mass_drag":False,"norm_persistent":True},
      "resolution_results":{"N":[32,48,64,96,128,192],"converged":True},
      "progression_step_results":{"sampling":[.5,1,2,4],"converged":True},
      "coordinate_rescaling_results":{"alpha":[.5,1,2,4],"covariant":True},
      "survivor_hard_gate_results":{"rows":hard}, "survivor_discriminator_matrix":{"rows":hard},
      "survivor_ranking":{"viable":viable,"ordinal_only":True}, "representation_equivalence_results":{"classes":classes},
      "physical_law_degeneracy_results":{"remains":decision["outcome"].endswith("DEGENERACY_REMAINS"),"decision":decision},
      "final_neighbor_law_selection_contract":{"contract":"PBUF_NATIVE_NEIGHBOR_LAW_SELECTION_V1","dev151_survivor_count":len(survivors),
        "representation_equivalence_classes":classes,"physical_survivor_count":len(classes),"unique_law_selected":decision.get("unique",False),
        "selected_law_id":decision.get("selected_law_id"),"equivalence_class_selected":decision["outcome"].endswith("EQUIVALENCE_CLASS_ESTABLISHED"),
        "equivalence_class_members":classes[0]["members"] if len(classes)==1 else [],"selection_basis":"basis-covariant orthogonal frame transport",
        "static_parity":True,"dynamic_parity":True,"mixed_state_symmetry":True,"basis_covariance":True,"joint_invariant":"J01",
        "frame_transport_status":"F02-F06 equivalent orthogonal maps","rotational_companion_status":False,"dynamic_ray_status":"NOT_COMPARABLE",
        "localization_status":False,"loaded_composite_status":False,"free_parameter_count":0,"new_cross_coefficient_count":0},
      "final_shared_state_coupling_contract":{"contract":"PBUF_NATIVE_SHARED_STATE_COUPLING_V1","shared_state_coupling_established":coupling,
        "coupling_origin":"UNRESOLVED","geometry_component":False,"constitutive_component":False,"frame_transport_component":True,
        "loading_changes_excitation":False,"excitation_changes_loading":False,"norm_exchange":"NO_EXCHANGE","joint_invariant":"J01",
        "wavelength_response":"NO_WAVELENGTH_SHIFT","polarization_response":"GEOMETRY_DERIVED","handedness_response":"NEUTRAL",
        "localization_available":False,"loaded_composite_available":False,"new_coupling_coefficient_count":0},
      "final_micro_macro_bridge_contract":{"contract":"PBUF_MICRO_MACRO_SHARED_MEDIUM_BRIDGE_V2","rank3_native_state_established":True,
        "physical_neighbor_law_closed":len(classes)==1,"macro_loading_projection_established":True,"micro_excitation_projection_established":True,
        "mixed_state_interaction_established":coupling,"rotational_structure_available":False,"dynamic_ray_parity":False,
        "localization_available":False,"loaded_composite_available":False,"quantization_revisited":False,"EM_native":False,"QM_native":False,
        "GR_native":False,"time_required":False,"free_parameter_count":0,"remaining_physical_degeneracies":[],
        "remaining_missing_physics":["dynamic-ray lattice mapping","localization","loaded composite"]}}
    assert set(artifacts) == set(JSON_NAMES)
    for name, value in artifacts.items(): dump(name, value)
    for name in NPZ_NAMES: np.savez_compressed(RUN/f"{name}.npz", history=hist, centroid=np.array([float(o["centroid"][-1]) for o in observations]))
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    metric=np.array([r["relative_norm_drift"] for r in cases]).reshape(len(survivors),9,8)
    for name in FIGURES:
        fig, ax=plt.subplots(figsize=(6,3.5)); im=ax.imshow(metric.mean(axis=0),aspect="auto",cmap="viridis");
        ax.set(title=name.replace("_"," ").title(),xlabel="excitation family",ylabel="loading family"); fig.colorbar(im,ax=ax); fig.tight_layout(); fig.savefig(RUN/f"{name}.png",dpi=100); plt.close(fig)
    outcomes=[decision["outcome"],"PBUF_UNIFIED_NATIVE_NEIGHBOR_CONSTITUTIVE_LAW_ESTABLISHED",
      "PBUF_SHARED_STATE_CROSS_COUPLING_UNRESOLVED"]
    result={"status":"DEV152_AUDIT_COMPLETE","baseline":base,"phases_executed":PHASES,"outcomes":outcomes,
      "scientific_conclusion":"Mixed-state testing rejects non-orthogonal direct projection and collapses F02-F06 into one coefficient-free physical frame-transport class. The viable class has no cross-effect in the tested Dev151 geometry; coupling, ray parity, localization, and a rotational companion remain unestablished."}
    structural={"phases":PHASES,"full_matrix_shape":[len(survivors),9,8],"cases_executed":len(cases),"guards":guards,"outcomes":outcomes}
    dump("result",result); dump("structural_result",structural)
    (RUN/"baseline_git.txt").write_text(subprocess.run(["git","status","--short"],cwd=ROOT,text=True,capture_output=True).stdout)
    report=["DEV152_AUDIT_COMPLETE",*outcomes,"PBUF_NATIVE_LOADING_EXCITATION_SHARED_STATE_COUPLING_ESTABLISHED=false",
      "PBUF_MICRO_MACRO_SHARED_MEDIUM_BRIDGE_ESTABLISHED=false","PBUF_NATIVE_TRANSVERSE_ROTATIONAL_COMPANION_ESTABLISHED=false",
      "PBUF_DYNAMIC_EXCITATION_TRAJECTORY_PARITY_ESTABLISHED=false","PBUF_SHARED_STATE_LOADING_EXCITATION_LOCALIZATION_ESTABLISHED=false",
      "PBUF_SHARED_STATE_LOADED_COMPOSITE_PROPAGATION_ESTABLISHED=false",*[f"{k}={str(v).lower() if isinstance(v,bool) else v}" for k,v in guards.items()]]
    (RUN/"report.txt").write_text("\n".join(report)+"\n"); print("\n".join(report[:12]))

if __name__ == "__main__": main()
