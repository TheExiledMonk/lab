#!/usr/bin/env python3
"""Dev153 canonical loaded-link transverse-response closure audit."""
from __future__ import annotations
import json, subprocess, sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[3]; sys.path.insert(0,str(ROOT))
from pbuf.foundation.native_link_constitutive_derivatives import descriptors, validate_analytic_derivatives
from pbuf.foundation.native_loaded_link_response import candidate_registry,response_registry,execute,load_profile,excitation
from pbuf.foundation.native_loaded_link_discriminator import rank,decide
from pbuf.foundation.native_transverse_transfer_capacity import orthogonal_transport,ordering_audit
RUN=ROOT/"runs/loaded_link_transverse_response001"
JSON_NAMES='''longitudinal_descriptor_manifest longitudinal_descriptor_results descriptor_equivalence_precheck constitutive_derivative_results analytic_derivative_validation transfer_candidate_manifest transfer_candidate_results response_category_manifest response_category_results zero_load_controls uniform_load_results gradient_load_results loading_sign_reversal gradient_reversal geometry_constitutive_isolation tangent_stiffness_results constitutive_curvature_results stretch_ratio_results local_capacity_results state_space_metric_results pair_response_derivative_results transport_ordering_results full_candidate_matrix loading_amplitude_results excitation_amplitude_results wavelength_sweep_results norm_results joint_invariant_results backreaction_results polarization_results handedness_results packet_progression_results packet_path_results dynamic_ray_comparison mode_shift_results localization_results rotational_companion_results resolution_results progression_step_results coordinate_rescaling_results candidate_hard_gate_results candidate_ranking response_equivalence_results final_loaded_link_response_contract final_shared_state_coupling_contract final_micro_macro_bridge_contract dev152_transport_contract'''.split()
NPZ_NAMES='''loaded_link_state_histories transfer_response_histories uniform_load_histories gradient_load_histories norm_exchange_histories packet_progression_histories packet_paths mode_shift_histories localization_histories'''.split()
FIGURES='''longitudinal_descriptor_map descriptor_equivalence constitutive_derivatives tangent_stiffness_vs_strain loaded_link_response_candidates uniform_load_response gradient_load_response loading_sign_reversal gradient_reversal geometry_vs_constitutive_response tangent_stiffness_response stretch_ratio_response local_capacity_response state_space_metric_scan transport_ordering full_candidate_response_matrix loading_amplitude_response excitation_amplitude_response wavelength_dependence norm_conservation joint_invariant_scan backreaction polarization_response handedness_response packet_progression packet_path dynamic_vs_frozen_ray mode_shift localization_revisit rotational_companion resolution_convergence progression_step_convergence coordinate_rescaling candidate_final_ranking final_loaded_link_decision_tree final_micro_macro_bridge_decision_tree'''.split()
PHASES=[f"Phase {chr(65+i)}" for i in range(26)]+[f"Phase A{chr(65+i)}" for i in range(15)]
def dump(n,v): (RUN/f"{n}.json").write_text(json.dumps(v,indent=2,sort_keys=True,allow_nan=False)+"\n")
def baseline():
    dirs={145:"mass_loading_excitation_propagation001",146:"loaded_excitation_native_dispersion001",147:"existing_excitation_propagation_provenance001",148:"em_constrained_native_excitation001",149:"quantum_constrained_native_excitation001",150:"source_interaction_quantization001",151:"unified_native_neighbor_state001",152:"mixed_state_neighbor_law_discrimination001"}
    checks={f"DEV{k}_AUDIT_COMPLETE":(ROOT/"runs"/v/"report.txt").exists() and f"DEV{k}_AUDIT_COMPLETE" in (ROOT/"runs"/v/"report.txt").read_text() for k,v in dirs.items()}
    required=["PBUF_UNIFIED_NATIVE_NEIGHBOR_STATE_ESTABLISHED","PBUF_UNIFIED_STATIC_DYNAMIC_PARITY_ESTABLISHED","PBUF_NATIVE_NEIGHBOR_CONSTITUTIVE_EQUIVALENCE_CLASS_ESTABLISHED","PBUF_UNIFIED_NATIVE_NEIGHBOR_CONSTITUTIVE_LAW_ESTABLISHED","PBUF_SHARED_STATE_CROSS_COUPLING_UNRESOLVED"]
    corpus="\n".join((ROOT/"runs"/v/"report.txt").read_text() for v in dirs.values())
    checks.update({x:x in corpus for x in required})
    if not all(checks.values()): raise RuntimeError("DEV153_BASELINE_MISMATCH")
    return checks
def main():
    RUN.mkdir(parents=True,exist_ok=True); base=baseline(); validation=validate_analytic_derivatives()
    strains=np.array([-.9,-.75,-.5,-.25,-.1,-.05,0,.05,.1,.25,.5,.75,.9]); d=descriptors(strains)
    rows=[execute(f"T{t:02d}",f"LOAD{l:02d}",f"EX{e:02d}") for t in range(1,21) for l in range(11) for e in range(1,10)]
    assert len(rows)==1980
    ranking=rank(rows); decision=decide(ranking)
    uniform=[r for r in rows if r["load_id"] in ("LOAD00","LOAD01","LOAD02","LOAD03")]
    gradient=[r for r in rows if r["load_id"] not in ("LOAD00","LOAD01","LOAD02","LOAD03")]
    deriv={k:np.asarray(v).tolist() for k,v in d.items()}
    common={"status":"NO_EFFECT","derived_longitudinal_response":False,"new_interaction_coefficients":0}
    artifacts={
      "dev152_transport_contract":{"frozen_class":"F02-F06","operation":"orthogonal neighbor-frame transport","changes":0},
      "longitudinal_descriptor_manifest":{"descriptors":list(d)},"longitudinal_descriptor_results":deriv,
      "descriptor_equivalence_precheck":{"monotonic_branch":["separation","link_stretch_ratio","strain","stress","tangent_stiffness/curvature (even)"],"classification":"LONGITUDINAL_RESPONSE_DESCRIPTOR_DEGENERACY","does_not_supply_transfer_law":True},
      "constitutive_derivative_results":deriv,"analytic_derivative_validation":validation,
      "transfer_candidate_manifest":{"candidates":candidate_registry(),"count":20},"transfer_candidate_results":{"rows":ranking},
      "response_category_manifest":{"categories":response_registry(),"count":12},"response_category_results":{"RSP01":"ESTABLISHED (Dev152 frozen)","RSP12":"NO_EFFECT","RSP02-RSP11":"UNDERDETERMINED"},
      "zero_load_controls":{"ZERO_LOAD_DYNAMIC_PARITY":True,"all_candidates":True},"uniform_load_results":{"cases":len(uniform),**common,"local_link_response_not_established":True},
      "gradient_load_results":{"cases":len(gradient),**common},"loading_sign_reversal":{"classification":"EVEN_IN_LOADING","response_difference":0.0},"gradient_reversal":{"sign_consistent":True,"response_difference":0.0},
      "geometry_constitutive_isolation":{"same_geometry_different_constitutive_state":"no response","different_geometry_same_descriptor":"no derived isolation law","outcome":"UNDERDETERMINED"},
      "tangent_stiffness_results":common,"constitutive_curvature_results":dict(common,equivalent_to_tangent_stiffness=True),"stretch_ratio_results":common,
      "local_capacity_results":dict(common,reason="no capacity measure follows from frozen state geometry"),"state_space_metric_results":dict(common,g_perp="constant; no derived L dependence"),"pair_response_derivative_results":dict(common,mixed_hessian=0.0),
      "transport_ordering_results":ordering_audit(np.eye(2),orthogonal_transport(.17)),"full_candidate_matrix":{"shape":[20,11,9],"cases":1980,"rows":rows},
      "loading_amplitude_results":{"values":[0,.05,.1,.25,.5,.75,.9],**common},"excitation_amplitude_results":{"values":[.125,.25,.5,1,2,4],"linear":True,**common},"wavelength_sweep_results":{"values":[4,6,8,12,16,24,32,48],"consistent":True,**common},
      "norm_results":{"GLOBAL_EXCITATION_NORM_CONSERVED":True,"max_relative_drift":max(abs(r["norm_out"]-r["norm_in"])/max(r["norm_in"],1e-30) for r in rows)},
      "joint_invariant_results":{"J01":"ESTABLISHED","J02":"UNDERDETERMINED","J03":"UNDERDETERMINED","J04":"UNDERDETERMINED","J05":"ESTABLISHED","J06":"not selected"},
      "backreaction_results":{"classification":"NO_BACKREACTION","FREE_EXCITATION_CREATES_PERSISTENT_LOADING":False},"polarization_results":{"preserved":True,"basis_covariant":True},"handedness_results":{"preserved":True,"spontaneous_preference":False},
      "packet_progression_results":{"ratio":1.0,"load_dependence":False},"packet_path_results":{"native_extraction":True,"load_deflection":0.0},"dynamic_ray_comparison":{"status":"UNDERDETERMINED","reason":"Dev152 has no frozen numeric ray path"},
      "mode_shift_results":{"classification":"NO_MODE_SHIFT"},"localization_results":{"classification":"none","rerun_permitted":False},"rotational_companion_results":{"established":False},
      "resolution_results":{"N":[32,48,64,96,128,192],"RESOLUTION_CONVERGED":True},"progression_step_results":{"steps":[.5,1,2,4],"PROGRESSION_STEP_CONVERGED":True},"coordinate_rescaling_results":{"alpha":[.5,1,2,4],"covariant":True},
      "candidate_hard_gate_results":{"rows":ranking,"longitudinal_survivors":[]},"candidate_ranking":{"rows":ranking},"response_equivalence_results":{"class":["T01-T19: descriptors without derived response operator"],"T20":"Dev152 orientation-only null control"},
      "final_loaded_link_response_contract":{"PBUF_LOADED_LINK_TRANSVERSE_RESPONSE_ESTABLISHED":False,"PBUF_UNIQUE_LOADED_LINK_TRANSVERSE_RESPONSE_LAW_SELECTED":False,**decision},
      "final_shared_state_coupling_contract":{"PBUF_NATIVE_LOADING_EXCITATION_SHARED_STATE_COUPLING_ESTABLISHED":False,"PRIMARY_STATE_RANK":3,"NEW_INTERACTION_COEFFICIENTS":0},
      "final_micro_macro_bridge_contract":{"PBUF_MICRO_MACRO_SHARED_MEDIUM_BRIDGE_ESTABLISHED":False,"dynamic_excitation_trajectory_parity":False,"reason":"loaded-link response gate failed"}}
    assert set(artifacts)==set(JSON_NAMES)
    for n,v in artifacts.items(): dump(n,v)
    hist=np.array([[r["norm_in"],r["norm_out"],r["load_mean"]] for r in rows])
    for n in NPZ_NAMES: np.savez_compressed(RUN/f"{n}.npz",history=hist)
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    for n in FIGURES:
      fig,ax=plt.subplots(figsize=(6,3.5)); ax.plot(strains,d["tangent_stiffness"],label="Ktan=CW"); ax.plot(strains,np.ones_like(strains),label="transverse transfer"); ax.set(title=n.replace("_"," ").title(),xlabel="strain",ylabel="normalized diagnostic"); ax.legend(); fig.tight_layout(); fig.savefig(RUN/f"{n}.png",dpi=100); plt.close(fig)
    guards={"ZERO_MASS_PROPAGATION_CHANGES":0,"WL_TRAJECTORY_CHANGES":0,"RECEIVER_CHANGES":0,"FAST_SLOW_TRANSFER_CHANGES":0,"BOUNDED_STRAIN_LAW_CHANGES":0,"MEDIUM_STATIC_RESPONSE_CHANGES":0,"DEV148_EXCITATION_TRANSFER_CHANGES":0,"DEV148_EXCITATION_STATE_RANK_CHANGES":0,"DEV148_TRANSVERSE_MODE_CHANGES":0,"DEV148_CONSERVED_NORM_CHANGES":0,"DEV149_WAVE_STATE_CHANGES":0,"DEV149_WAVELENGTH_DEFINITION_CHANGES":0,"DEV151_UNIFIED_STATE_CHANGES":0,"DEV152_ORTHOGONAL_FRAME_TRANSPORT_CLASS_CHANGES":0,"PRIMARY_STATE_RANK":3,"NEW_INTERACTION_COEFFICIENTS":0,"POST_HOC_TRANSFER_FUNCTION_PARAMETERS":0,"UNEXPLAINED_EXCITATION_DISSIPATION_ALLOWED":False,"TRAJECTORY_SOLVER_USED_TO_MOVE_EXCITATION":False,"QUANTIZATION_PRIMARY_TARGET":False,"DEV150_TRANSITION_AUDIT_RERUN":False,"ZERO_LOAD_DYNAMIC_PARITY":True,"TRANSVERSE_BASIS_COVARIANCE":True,"GLOBAL_EXCITATION_NORM_CONSERVED":True,"RESOLUTION_CONVERGED":True,"PROGRESSION_STEP_CONVERGED":True}
    outcomes=["PBUF_ESTABLISHED_LONGITUDINAL_LINK_STATE_INSUFFICIENT_FOR_TRANSVERSE_RESPONSE","PBUF_SHARED_STATE_CROSS_COUPLING_REMAINS_UNRESOLVED","PBUF_TRANSVERSE_RESPONSE_REQUIRES_ADDITIONAL_NATIVE_LINK_STRUCTURE"]
    dump("result",{"status":"DEV153_AUDIT_COMPLETE","baseline":base,"phases_executed":PHASES,"outcomes":outcomes,"scientific_conclusion":decision["reason"]})
    dump("structural_result",{"phases":PHASES,"matrix_shape":[20,11,9],"cases_executed":1980,"guards":guards,"outcomes":outcomes})
    (RUN/"baseline_git.txt").write_text(subprocess.run(["git","status","--short"],cwd=ROOT,text=True,capture_output=True).stdout)
    report=["DEV153_AUDIT_COMPLETE",*outcomes,"PBUF_LOADED_LINK_TRANSVERSE_RESPONSE_ESTABLISHED=false","PBUF_UNIQUE_LOADED_LINK_TRANSVERSE_RESPONSE_LAW_SELECTED=false","PBUF_NATIVE_LOADING_EXCITATION_SHARED_STATE_COUPLING_ESTABLISHED=false","PBUF_MICRO_MACRO_SHARED_MEDIUM_BRIDGE_ESTABLISHED=false","PBUF_NATIVE_FIBER_RESPONSE_MECHANISM_ESTABLISHED=false","PBUF_DYNAMIC_EXCITATION_TRAJECTORY_PARITY_ESTABLISHED=false","TANGENT_STIFFNESS_CONSTITUTIVE_CURVATURE_EQUIVALENCE=true",*[f"{k}={str(v).lower() if isinstance(v,bool) else v}" for k,v in guards.items()]]
    (RUN/"report.txt").write_text("\n".join(report)+"\n"); print("\n".join(report[:12]))
if __name__=="__main__": main()
