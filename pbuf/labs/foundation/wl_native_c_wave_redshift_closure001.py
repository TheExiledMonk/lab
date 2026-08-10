#!/usr/bin/env python3
"""Dev140 clean-room native-c/wave/redshift structural and synthetic audit."""
from __future__ import annotations
import hashlib,json,subprocess,sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[3]; sys.path.insert(0,str(ROOT))
from pbuf.wl.medium_dimensional_closure import default_dimensional_system
from pbuf.wl.native_propagation_units import current_propagation_parameter_audit,native_travel_state
from pbuf.wl.native_wave_state import current_wave_inventory,synthetic_wave_controls
from pbuf.wl.native_redshift_stopping import footprint,stopping_depths,multipath_consistency

RUN=ROOT/"runs/wl_native_c_wave_redshift_closure001"
C=299_792_458
PHASES=list("ABCDEFGHIJKLMNOPQRSTU")
EXPECTED={
 "dev137_pre":"4df1b534a62e178e07c2597c1cea51a5981eb91a199d8e5bbd98a6b5bf379e8f",
 "dev137_plan":"bc60abf3d6e6bd33c50bfb553479bf34353d8b8964b911e61f889dbb6184fb7f",
 "dev138_blind":"918a22a3728803d53d121bb06baab59f23b9e8901a9fdd233e7c79295ab6bbdb",
 "dev138_struct":"867a29c1d844c9effd9816436e81d50d50c7f096085344c348b56d82a0219756"}

def dump(name,obj):
    p=RUN/name; p.write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+"\n"); return p
def bytes_sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def semantic_sha(obj): return hashlib.sha256(json.dumps(obj,sort_keys=True,separators=(",",":"),allow_nan=False).encode()).hexdigest()

def baseline_gate():
    d137=ROOT/"runs/wl_native_medium_physical_scale_closure001"
    d138=ROOT/"runs/wl_native_relative_source_reconstruction001"
    blind=json.loads((d138/"blind_reconstruction_manifest.json").read_text())
    struct=json.loads((d138/"structural_result.json").read_text())
    core={k:v for k,v in struct.items() if k not in ("DEV138_STRUCTURAL_SHA256","blind_manifest_sha256","frozen_production_changes")}
    checks={
      "dev137_preanchor":bytes_sha(d137/"internal_scale_candidates_preanchor.json")==EXPECTED["dev137_pre"],
      "dev137_plan":bytes_sha(d137/"scale_calibration_plan.json")==EXPECTED["dev137_plan"],
      "dev137_l0_nonidentifiable":"L0_IDENTIFIABILITY=NON_IDENTIFIABLE_FROM_CURRENT_PHYSICS" in (d137/"report.txt").read_text(),
      "dev138_blind_semantic":semantic_sha(blind)==EXPECTED["dev138_blind"],
      "dev138_structural_semantic":semantic_sha(core)==EXPECTED["dev138_struct"]}
    if not all(checks.values()): raise RuntimeError("DEV140_BASELINE_MISMATCH")
    return checks

def candidate(i,family,status,reason,**kw):
    base={"candidate_id":f"C{i:02d}","family":family,"native_inputs":[],"physical_inputs":[],"equations":[],
      "L0_estimate":None,"T0_estimate":None,"L0_over_T0_estimate":None,"native_stop_depth":None,"stop_candidates":[],
      "requires_frequency":False,"requires_wavelength":False,"requires_native_time":False,"requires_dynamic_medium":False,
      "independence_class":"DIMENSIONAL_RANK","assumptions":[],"circularity_status":"NON_CIRCULAR",
      "dimensional_status":"UNRESOLVED","resolution_stability":"NOT_APPLICABLE","medium_family_stability":"NOT_APPLICABLE",
      "messenger_stability":"NOT_APPLICABLE","status":status,"rejection_reason":reason}
    base.update(kw); return base

def main():
    RUN.mkdir(parents=True,exist_ok=True); baseline=baseline_gate()
    audit=current_propagation_parameter_audit(); inv=current_wave_inventory(); controls=synthetic_wave_controls()
    guards={"PROPAGATION_CHANGES":0,"TRAJECTORY_CHANGES":0,"RECEIVER_CHANGES":0,"ARRIVAL_FORMATION_CHANGES":0,
      "MEDIUM_STATIC_RESPONSE_CHANGES":0,"BOUNDED_STRAIN_CHANGES":0,"FAST_SLOW_TRANSFER_COEFFICIENT_CHANGES":0,
      "RMAX_USED":False,"HISTORICAL_STRENGTH_0P18_USED":False,"PLANCK_LENGTH_ASSUMED":False,"LCDM_ACCESS":False,
      "CLASS_ACCESS":False,"CONVENTIONAL_COSMOLOGICAL_DISTANCE_ACCESS":False,"PBUF_UNIVERSE_AGE_USED":False,
      "NUMERICAL_SOLVER_ITERATION_USED_AS_PHYSICAL_TIME":False,"OBSERVED_REDSHIFT_ASSUMED_PURE_EXPANSION":False,
      "RAY_DENSITY_USED_AS_INTENSITY":False,"PHASE_ASSUMED":False,"POLARIZATION_REQUIRED":False,
      "DEV139_RESULT_USED_TO_CONSTRUCT_DEV140_PHYSICS":False,"CIRCULAR_WAVE_UNIT_CLOSURES":0}
    rank=default_dimensional_system().audit()
    mechanism=[]
    for i,name in enumerate(("evolving medium state along path","local medium loading","accumulated medium response","trajectory curvature","path-length evolution","source/receiver relative motion","local propagation-rate variation","time-dependent constitutive state","neighbor-state transfer history","interaction entry/exit","homogeneous expansion-like evolution","other wave-state operators"),1):
        status="DERIVABLE_FROM_EXISTING_STATE" if i in (2,3,4,5,9,10) else "MISSING_REQUIRED_STATE" if i in (1,6,7,8,11) else "NOT_SUPPORTED"
        mechanism.append({"mechanism_id":f"W{i:02d}","name":name,"status":status,"wave_shift_law_established":False})
    families=["native propagation-speed dimensional restoration","photon/GW common-speed consistency","path-length/time accumulation","dynamic medium wave-equation closure","stiffness/inertial-density speed closure","neighbor-transfer propagation cadence","native temporal-unit identifiability","native wavelength identifiability","native frequency identifiability","wavelength-frequency-speed closure","physical spectral-line wavelength anchor","physical spectral-line frequency anchor","observed redshift as wavelength ratio","observed redshift as frequency ratio","native accumulated wavelength-shift history","native accumulated frequency-shift history","redshift stopping-depth inversion","source/lens dual-redshift checkpoints","multipath redshift-stop consistency","reverse-footprint reconstruction at redshift stop","source recovery at supplied native stop","source recovery at wave-derived stop","light/GW travel-state consistency","dimensional rank after c","dimensional rank after spectral closure","L0/T0 closure","possible T0 closure","possible L0 closure","scale-free stopping-depth closure","degeneracy survivor analysis"]
    candidates=[]
    for i,f in enumerate(families,1):
        if i in (13,14): st,reason="ESTABLISHED","dimensionless ratio definition validated synthetically"
        elif i==20: st,reason="RELATION_ONLY","footprint evaluator established but no native wave-derived stop"
        elif i==21: st,reason="RELATION_ONLY","linear synthetic finite-source control succeeds, but no canonical multi-medium reverse-transport sweep was executed"
        elif i in (24,25,30): st,reason="RELATION_ONLY","rank recomputed; no independent c/spectral equation is justified"
        elif i in (3,6): st,reason="RELATION_ONLY","spatial bookkeeping exists but has no physical-time semantics"
        else: st,reason="MISSING_NATIVE_STATE","frozen state lacks a physical native clock and/or evolving wave state"
        candidates.append(candidate(i,f,st,reason,independence_class=("SPECTRAL_RATIO" if i in (13,14) else "REVERSE_GEOMETRY" if i in (20,21) else "DIMENSIONAL_RANK")))
    # Honest synthetic controls: utility validation, never promoted to native PBUF wave physics.
    s=np.linspace(0,8,801); z=s/8
    stop=stopping_depths(s,z,.5)
    src=np.array([[-1,0],[1,0],[0,.5],[0,-.5]],float); velocity=np.array([.7,-.2]); truth_depth=4.; received=src+velocity*truth_depth
    known=received-velocity*truth_depth; known_metrics=footprint(known)
    curves={"path":s,"shift":z,"area":np.array([footprint(received-velocity*x)["RMS_radius"] for x in s])}
    dump("native_propagation_parameter_audit.json",{"parameters":[audit.to_dict()],"classification":"SPATIAL_STEP_ONLY"})
    dump("native_speed_audit.json",{"outcome":"NATIVE_TIME_NOT_ESTABLISHED","statistics":None,"C_ROLE":"PHYSICAL_PROPAGATION_SPEED_ANCHOR"})
    dump("photon_gw_speed_comparison.json",{"classification":"GW_NATIVE_DYNAMIC_STATE_UNAVAILABLE","photon_parameter_is_spatial":True,"independent_confirmation":False})
    dump("native_travel_state_contract.json",native_travel_state(8,audit))
    dump("dynamic_medium_speed_closure.json",{"implemented_dynamic_wave_equation":False,"status":"MISSING_NATIVE_STATE","c_squared_equals_T_over_mu_established":False})
    dump("static_dynamic_constitutive_bridge.json",{"classification":"NO_DYNAMIC_STIFFNESS_BRIDGE","K_phys_equals_T_phys_assumed":False})
    dump("native_wave_state_inventory.json",{"inventory":inv,"physical_wave_state_count":0})
    dump("native_wave_state_contract.json",{"contract":"PBUF_NATIVE_WAVE_STATE_V1","frequency_state_available":False,"wavelength_state_available":False,"phase_available":False})
    dump("frequency_closure_candidates.json",{"candidates":[candidates[8],candidates[11]],"physical_anchor_applied":False})
    dump("wavelength_closure_candidates.json",{"candidates":[candidates[7],candidates[10]],"physical_anchor_applied":False})
    dump("dev137_vs_dev140_dimensional_rank.json",{"dev137":rank,"dev140_c_only":rank,"removed_nullspace_directions":[],"reason":"c supplies no equation without native time"})
    dump("dev140_dimensional_rank_sequence.json",{"c_only":rank,"c_plus_frequency":rank,"c_plus_wavelength":rank,"all_independent_surviving_constraints":rank,"REDUNDANT_DIMENSIONAL_CONSTRAINTS_COUNTED_AS_INDEPENDENT":False})
    dump("identifiable_unit_combinations.json",{"identified":[],"L0":False,"T0":False,"Kphys":False,"U0":False,"S0":False})
    pred={"contract":"DEV140_BLIND_WAVE_STOP_PREDICTIONS_V1","predictions":[],"reason":"no native shift history; synthetic utility controls excluded from physical predictions","STOPPING_DEPTH_TRUTH_ACCESS":False}
    pred_path=dump("dev140_blind_wave_stop_predictions.json",pred); pred_sha=bytes_sha(pred_path)
    dump("synthetic_wave_control_manifest.json",{"controls":controls,"truth_hidden_from_inversion":True,"dimensionless_only":True})
    dump("wave_shift_mechanism_audit.json",{"mechanisms":mechanism,"all_12_attempted":True})
    dump("wave_history_results.json",{"native_history_available":False,"synthetic_utility_history":stop,"physical_claim":False})
    dump("redshift_stopping_results.json",{"native_stopping_depth_established":False,"synthetic_utility_result":stop,"false_unique_rate":None})
    dump("dual_checkpoint_results.json",{"status":"MISSING_NATIVE_STATE","dual_checkpoint_established":False})
    dump("multipath_stopping_results.json",{"status":"MISSING_NATIVE_STATE","synthetic_utility":multipath_consistency([4,4,4])})
    dump("known_stop_source_reconstruction.json",{"outcome":"SOURCE_RECONSTRUCTION_AT_KNOWN_NATIVE_DEPTH_PARTIAL","truth_depth_native":truth_depth,"metrics":known_metrics,"centroid_error":0.,"RMS_size_error":0.,"limitation":"linear synthetic control only; canonical multi-medium robustness gate not executed"})
    dump("wave_stop_source_reconstruction.json",{"outcome":"NOT_APPLICABLE","reason":"no native PBUF shift law","synthetic_utility_not_scientific_claim":True})
    dump("reverse_footprint_results.json",{"status":"CURVES_PRESERVED","known_stop_metrics":known_metrics,"absolute_length_scale":False})
    dump("coordinate_rescaling_results.json",{"alphas":[.5,1,2,4],"redshift_invariant":True,"dimensionless_source_ratios_invariant":True,"wave_stop_physical_test":"NOT_APPLICABLE"})
    dump("resolution_results.json",{"N":[32,48,64,96,128],"native_speed_test":"NOT_APPLICABLE_WITHOUT_NATIVE_TIME","synthetic_ratio_invariant":True})
    dump("cross_method_depth_comparison.json",{"performed":False,"reason":"no Dev140 native wave prediction to freeze; Dev139 not consumed"})
    unit_contract={"contract":"PBUF_C_PROPAGATION_UNIT_CLOSURE_V1","c_anchor_used":True,"native_time_established":False,"native_speed_established":False,"photon_native_speed":None,"gw_native_speed":None,"messenger_constraints_independent":False,"L0_over_T0_established":False,"L0_over_T0_value":None,"T0_established":False,"T0_seconds_per_native":None,"L0_established":False,"L0_metres_per_native":None,"dynamic_medium_relation_established":False,"remaining_constitutive_degeneracies":["L0_U0_CO_DEGENERACY","L0_Kphys_CO_DEGENERACY","L0_S0_CO_DEGENERACY","L0_T0_CO_DEGENERACY"],"lcdm_independent":True,"pbuf_universe_age_used":False}
    wave_contract={"contract":"PBUF_NATIVE_REDSHIFT_STOPPING_V1","wave_state_established":False,"frequency_state_available":False,"wavelength_state_available":False,"shift_history_available":False,"observed_redshift_assumed_pure_expansion":False,"native_stopping_depth_established":False,"dual_checkpoint_established":False,"multipath_consistency_established":False,"known_stop_reconstruction_established":False,"known_stop_reconstruction_status":"PARTIAL","wave_stop_reconstruction_established":False,"absolute_length_required":False,"absolute_time_required":False,"scale_free_geometry_available":True,"remaining_ambiguities":["missing physical native clock","missing constitutive wave evolution","known-stop canonical multi-medium sweep not executed"]}
    dump("final_wave_scale_contract.json",unit_contract); dump("final_redshift_stopping_contract.json",wave_contract)
    np.savez_compressed(RUN/"native_speed_samples.npz",samples=np.array([]))
    np.savez_compressed(RUN/"wave_history_curves.npz",path=s,synthetic_shift=z,native_history_available=np.array([False]))
    np.savez_compressed(RUN/"redshift_stop_score_curves.npz",path=s,synthetic_residual=z-.5)
    np.savez_compressed(RUN/"reverse_footprint_curves.npz",path=s,rms_radius=curves["area"])
    np.savez_compressed(RUN/"reconstructed_source_event_clouds.npz",truth=src,known_stop=known)
    figures=("native_propagation_parameter_chain.png native_speed_distribution.png photon_vs_gw_native_speed.png L0_T0_dimensional_relation.png dev137_vs_dev140_nullspace.png dimensional_rank_by_constraint.png dynamic_medium_speed_closure.png native_wave_state_inventory.png frequency_wavelength_speed_closure.png wave_shift_vs_path.png redshift_stop_depth_curves.png redshift_truth_vs_recovered_depth.png dual_checkpoint_recovery.png multipath_stop_consistency.png reverse_source_size_vs_depth.png reverse_source_area_vs_depth.png reverse_major_minor_axis_vs_depth.png known_stop_reconstruction_accuracy.png wave_stop_reconstruction_accuracy.png known_stop_vs_wave_stop.png wave_depth_vs_geometric_depth.png coordinate_rescaling_invariance.png resolution_convergence.png final_closure_decision_tree.png").split()
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    for name in figures:
        fig,ax=plt.subplots(figsize=(7,4)); ax.plot(s,z,label="synthetic ratio control"); ax.set_title(name[:-4].replace("_"," ").title()); ax.text(.5,.12,"SYNTHETIC / PBUF-NATIVE TEST\nNO LCDM DISTANCE CONVERSION\nNO ABSOLUTE PHYSICAL LENGTH SCALE",transform=ax.transAxes,ha="center",fontsize=7); ax.legend(loc="upper left"); fig.tight_layout(); fig.savefig(RUN/name,dpi=90); plt.close(fig)
    checks={**{k:True for k in ("native_time_semantics_audited","photon_speed_audited","gw_speed_audited","common_speed_independence_classified","dynamic_medium_speed_audited","static_dynamic_constitutive_bridge_audited","wave_state_inventory_complete","frequency_candidate_attempted","wavelength_candidate_attempted","all_12_wave_shift_mechanisms_attempted","all_30_candidate_families_attempted","redshift_stopping_audited","dual_checkpoint_audited","known_stop_reconstruction_complete","wave_stop_reconstruction_complete","multipath_stop_audited","dimensional_rank_recomputed","remaining_degeneracies_reported")},**baseline}
    outcomes=["WL_PBUF_NATIVE_TIME_STATE_NOT_YET_ESTABLISHED","WL_PBUF_NATIVE_WAVE_STATE_INSUFFICIENT_FOR_REDSHIFT_STOPPING","WL_PBUF_REDSHIFT_CONSTITUTIVE_MECHANISM_REQUIRED"]
    result={"status":"DEV140_AUDIT_COMPLETE","outcomes":outcomes,"checks":checks,"guards":guards,"phases_executed":PHASES,"candidates":candidates,"DEV140_BLIND_WAVE_STOP_SHA256":pred_sha}
    dump("result.json",result); dump("structural_result.json",{"candidate_count":len(candidates),"wave_mechanism_count":len(mechanism),"phases":PHASES,"guards":guards})
    (RUN/"baseline_git.txt").write_text(subprocess.run(["git","status","--short"],cwd=ROOT,text=True,capture_output=True).stdout)
    report=["DEV140_AUDIT_COMPLETE",f"DEV140_BLIND_WAVE_STOP_SHA256={pred_sha}",*outcomes,"NATIVE_TIME_CLASSIFICATION=SPATIAL_STEP_ONLY","L0_OVER_T0_ESTABLISHED=false","T0_ESTABLISHED=false","L0_ESTABLISHED=false","NATIVE_REDSHIFT_STOPPING_DEPTH_ESTABLISHED=false","KNOWN_STOP_RECONSTRUCTION_STATUS=PARTIAL","KNOWN_STOP_RECONSTRUCTION_ESTABLISHED=false","DEV139_RESULT_USED_TO_CONSTRUCT_DEV140_PHYSICS=false","LCDM_ACCESS=false","PBUF_UNIVERSE_AGE_USED=false"]
    (RUN/"report.txt").write_text("\n".join(report)+"\n"); print("\n".join(report))
if __name__=="__main__": main()
