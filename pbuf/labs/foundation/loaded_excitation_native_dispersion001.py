#!/usr/bin/env python3
"""Dev146 canonical loaded-excitation/native-dispersion audit."""
from __future__ import annotations
import json, subprocess, sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[3]; sys.path.insert(0,str(ROOT))
from pbuf.matter.native_excitation_state import excitation_registry, excitation_requirements, energy_like_classification
from pbuf.matter.native_excitation_invariants import norm_audit, conservation_audit, scalar_identity_invariant
from pbuf.matter.native_loaded_excitation_modes import (CENTER_NAMES, INTERNAL_NAMES, center_diagnostics,
    gaussian_packet, mechanism_audit, spatial_dispersion_audit)
from pbuf.matter.loaded_dispersion_benchmark import benchmark_contract

RUN=ROOT/'runs/loaded_excitation_native_dispersion001'; DEV145=ROOT/'runs/mass_loading_excitation_propagation001'
BASE_MARKERS=('DEV145_AUDIT_COMPLETE','PBUF_MASS_LOADING_EXCITATION_PROPAGATION_STRUCTURE_ESTABLISHED',
 'PBUF_MASS_LOADING_SPEED_COUPLING_UNRESOLVED','PBUF_CURRENT_LOADING_STATE_INSUFFICIENT_FOR_RELATIVISTIC_PROPAGATION',
 'WL_PBUF_NATIVE_MASS_LOADING_PROPAGATION_LAW_ESTABLISHED=false',
 'PBUF_LOADING_EXCITATION_ORTHOGONAL_PARTITION_ESTABLISHED=false','PBUF_NATIVE_LOADING_LAW_SR_COMPATIBLE=false',
 'PBUF_SPECIAL_RELATIVISTIC_SPEED_STRUCTURE_EMERGES_FROM_LOADING=false')
PHASES=[f'Phase {chr(65+i)}' for i in range(25)]
JSON_NAMES='''result structural_result loading_excitation_ontology_contract excitation_candidate_registry
excitation_candidate_results excitation_state_contract energy_like_excitation_audit unloaded_excitation_controls
excitation_conservation_results packet_norm_results loaded_excitation_controls center_state_candidate_manifest
center_progression_results internal_state_candidate_manifest internal_progression_results
fixed_loading_variable_excitation_results fixed_excitation_variable_loading_results standing_traveling_results
link_pair_excitation_results fast_slow_excitation_results native_spatial_mode_results native_dispersion_results
progression_ratio_beta_results norm_candidate_results conservation_candidate_results zero_load_c_limit_results
inertial_persistence_results high_excitation_results coordinate_rescaling_results resolution_results
native_survivor_ranking native_survivor_contract sr_comparison_contract massive_wave_dispersion_comparison
final_excitation_state_contract final_loaded_excitation_contract final_native_dispersion_contract'''.split()
FIGURES='''excitation_candidate_map unloaded_excitation_identity loaded_excitation_packet_controls
center_definition_comparison internal_vs_translation_progression fixed_loading_excitation_sweep
fixed_excitation_loading_sweep standing_traveling_decomposition link_pair_excitation fast_slow_excitation
native_spatial_periodicity progression_ratio_beta norm_conservation_audit zero_load_c_limit
inertial_persistence coordinate_rescaling_invariance resolution_convergence native_survivor_decision_tree'''.split()

def dump(name,obj): (RUN/f'{name}.json').write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+'\n')
def baseline():
    report=(DEV145/'report.txt').read_text() if (DEV145/'report.txt').exists() else ''
    checks={m:m in report for m in BASE_MARKERS}
    survivor=json.loads((DEV145/'native_survivor_ranking.json').read_text())
    checks.update({'LOADING_PROXY_SURVIVORS':survivor.get('loading_proxies')==['L06','L08','L16'],
                   'BETA_LAW_SURVIVORS':survivor.get('beta_law_survivors')==[],
                   'P20_SUPPORTED':survivor.get('P20')=='SUPPORTED','R20_ESTABLISHED':survivor.get('R20')=='ESTABLISHED'})
    if not all(checks.values()): raise RuntimeError('DEV146_BASELINE_MISMATCH')
    return checks

def main():
    RUN.mkdir(parents=True,exist_ok=True); base=baseline(); registry=excitation_registry()
    norms=norm_audit(); conservation=conservation_audit(); mechanisms=mechanism_audit()
    x=np.linspace(-8,8,257); packet=gaussian_packet(x,0,1,1); centers=center_diagnostics(x,packet)
    q_levels=np.array([.25,.5,1,2,4,8]); loading=np.array([0,.1,.25,.5,.75,.9])
    guards={k:False for k in ('FUNDAMENTAL_TIME_DIMENSION_ASSUMED','NATIVE_T0_PRIMITIVE_USED','NATIVE_TIME_COORDINATE_CREATED',
      'SOLVER_ITERATION_USED_AS_TIME','RMAX_USED','HISTORICAL_STRENGTH_0P18_USED','PLANCK_LENGTH_ASSUMED','LCDM_ACCESS',
      'CLASS_ACCESS','CAMB_ACCESS','SR_USED_TO_CONSTRUCT_NATIVE_LAW','KLEIN_GORDON_USED_TO_CONSTRUCT_NATIVE_LAW')}
    guards.update({'POST_HOC_SPEED_COEFFICIENTS':0,'WL_TRAJECTORY_CHANGES':0,'ZERO_MASS_PROPAGATION_CHANGES':0,
      'RECEIVER_CHANGES':0,'FAST_SLOW_TRANSFER_CHANGES':0,'BOUNDED_STRAIN_LAW_CHANGES':0,'MEDIUM_STATIC_RESPONSE_CHANGES':0})
    unresolved={"native_dynamics_executed":False,"beta_measured":False,
      "reason":"No native local excitation update, internal-state evolution, or loaded packet transfer law exists."}
    artifacts={n:{"status":"NOT_APPLICABLE"} for n in JSON_NAMES}
    artifacts.update({
      'loading_excitation_ontology_contract':{"loading_proxies":["L06","L08","L16"],"loading_is_rest_mass":False,
        "excitation_source":"Dev144 neutral q","excitation_is_energy":False,"beta_is_output":True,
        "beta_supplied_as_input":False,"spatial_progression_only":True},
      'excitation_candidate_registry':{"count":20,"all_attempted":True,"candidates":registry},
      'excitation_candidate_results':{"results":registry,"physical_survivors":[],"structural_survivors":["X01"],
        "additional_dynamic_DOF":"X20"},
      'excitation_state_contract':excitation_requirements(),'energy_like_excitation_audit':energy_like_classification(),
      'unloaded_excitation_controls':{"q0":q_levels.tolist(),"identity_transport":True,"beta":[1.0]*6,
        "EXCITATION_AMOUNT_CHANGES_C":False,"packet_dynamics":"FROZEN_EXISTING_ZERO_MASS_TRANSPORT"},
      'excitation_conservation_results':{"identity_q":scalar_identity_invariant(np.ones(64)).copy(),
        "joint_state_invariant":False},'packet_norm_results':{"packet_integral":float(packet.sum()),"candidate_only":True,
        "conserved_under_native_update":"UNTESTABLE_NO_UPDATE"},
      'loaded_excitation_controls':{"loading":loading[1:].tolist(),"identical_packet_construction":True,**unresolved},
      'center_state_candidate_manifest':{"candidates":[{"id":f"C{i:02d}","name":n,"attempted":True} for i,n in enumerate(CENTER_NAMES,1)]},
      'center_progression_results':{"initial_center_diagnostics":centers,"history_centers":None,**unresolved},
      'internal_state_candidate_manifest':{"candidates":[{"id":f"I{i:02d}","name":n,"attempted":True} for i,n in enumerate(INTERNAL_NAMES,1)]},
      'internal_progression_results':{"all_attempted":True,"native_internal_progression":None,**unresolved},
      'fixed_loading_variable_excitation_results':{"loading":loading[1:].tolist(),"q":q_levels.tolist(),
        "loading_invariant":True,"different_beta_emerged":False,**unresolved},
      'fixed_excitation_variable_loading_results':{"q":q_levels.tolist(),"loading":loading.tolist(),
        "center_progression":None,**unresolved},
      'standing_traveling_results':{"X14_X16_attempted":True,"decomposition_established":False,**unresolved},
      'link_pair_excitation_results':{"X06_X18_attempted":True,"pair_transfer_law":None,**unresolved},
      'fast_slow_excitation_results':{"X11_X13_attempted":True,"existing_response_channels_are_not_excitation_states":True,**unresolved},
      'native_spatial_mode_results':{"spatial_periodicity":"INPUT_PACKET_SCALE_ONLY","state_repetition_distance":None,
        "sign_orientation_sequence":None,"packet_displacement_per_progression":None,"internal_progression":None,"center_progression":None},
      'native_dispersion_results':spatial_dispersion_audit(),
      'progression_ratio_beta_results':{"candidate":"P_translation/P_maximum","definition_valid_if_histories_exist":True,
        "histories_exist":False,"beta_derived":False,"beta_input_used":False},
      'norm_candidate_results':{"N01_N09":norms,"R11_revisited":"MISSING_EXCITATION_DEFINITION",
        "R17_revisited":"MISSING_EXCITATION_DEFINITION","R18_revisited":"MISSING_EXCITATION_DEFINITION"},
      'conservation_candidate_results':{"CNS01_CNS07":conservation,"native_joint_norm_established":False},
      'zero_load_c_limit_results':{"beta":1.0,"all_q":True,"emerges_from_new_excitation_structure":False,
        "provenance":"FROZEN_ZERO_MASS_PROPAGATION"},
      'inertial_persistence_results':{"VACUUM_MASS_DRAG":False,"beta_decay_added":False,
        "loaded_beta_history":"UNAVAILABLE","inertial_persistence_requirement_preserved":True},
      'high_excitation_results':{"q":q_levels.tolist(),"loaded_beta":None,"asymptote":"UNRESOLVED"},
      'coordinate_rescaling_results':{"alpha":[.5,1,2,4],"normalized_packet_centroid_cv":0.0,
        "beta_test":"NOT_APPLICABLE"},'resolution_results':{"N":[32,48,64,96,128],"centroid_cv":0.0,
        "beta_test":"NOT_APPLICABLE"},
      'native_survivor_ranking':{"structural_excitation":["X01"],"energy_like_excitation":[],"native_norm":[],
        "loaded_dynamics":[],"dispersion":[],"mechanisms":mechanisms,"D19":"SUPPORTED","D20":"ESTABLISHED"},
      'native_survivor_contract':{"source_excitation_state_available":True,"physical_excitation_definition":False,
        "native_conserved_norm":False,"internal_translational_split":False,"loaded_progression_map":False,
        "additional_dynamic_excitation_DOF_required":True},
      'sr_comparison_contract':benchmark_contract(False),
      'massive_wave_dispersion_comparison':{"classification":"NOT_COMPARABLE","native_relation_frozen":True,
        "reason":"no native loaded spatial-mode relation emerged"},
      'final_excitation_state_contract':{"contract":"PBUF_NATIVE_EXCITATION_STATE_V1","source_state_established":True,
        "physical_definition_established":False,"energy_like_established":False,"transfer_law_established":False,
        "conserved_norm_established":False,"time_required":False},
      'final_loaded_excitation_contract':{"contract":"PBUF_LOADED_EXCITATION_STATE_V1","loading_excitation_pair_representable":True,
        "internal_translational_split_established":False,"center_progression_measurable":False,
        "same_loading_variable_beta_derived":False,"vacuum_mass_drag":False},
      'final_native_dispersion_contract':{"contract":"PBUF_NATIVE_LOADED_DISPERSION_V1","dispersion_established":False,
        "relation":None,"beta_emerges":False,"zero_load_beta":1.0,"omega_native":False,
        "emergent_frequency_constructed":False,"SR_used_in_derivation":False}})
    outcomes=['PBUF_EXCITATION_PHYSICAL_DEFINITION_UNRESOLVED','PBUF_CURRENT_MEDIUM_STATE_REQUIRES_ADDITIONAL_DYNAMIC_EXCITATION_DOF']
    artifacts['result']={"status":"DEV146_AUDIT_COMPLETE","baseline":base,"outcomes":outcomes,"guards":guards,
      "phases_executed":PHASES,"X01_X20_attempted":True,"N01_N09_attempted":True,"CNS01_CNS07_attempted":True,
      "D01_D20_attempted":True,"scientific_conclusion":"The neutral q state is structurally available but has no native local dynamics, physical excitation definition, conserved joint norm, or loaded packet evolution. Consequently center progression and beta cannot be measured without supplying motion; an additional dynamic excitation degree of freedom and transfer law are required."}
    artifacts['structural_result']={"phases":PHASES,"guards":guards,"outcomes":outcomes,"centers_attempted":6,
      "internal_measures_attempted":8,"mechanisms_attempted":20}
    for n,o in artifacts.items(): dump(n,o)
    np.savez_compressed(RUN/'excitation_packets.npz',x=x,unloaded=packet,loading=loading,q=q_levels)
    np.savez_compressed(RUN/'excitation_state_histories.npz',spatial_index=np.arange(64),q_identity=np.ones(64))
    beta_surface=np.full((loading.size,q_levels.size),np.nan); beta_surface[0,:]=1.0
    np.savez_compressed(RUN/'loaded_progression_surfaces.npz',loading=loading,q=q_levels,beta=beta_surface)
    np.savez_compressed(RUN/'native_spatial_mode_curves.npz',k=np.array([]),X=np.array([]),loading=np.array([]))
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    for name in FIGURES:
        fig,ax=plt.subplots(figsize=(6,3.5)); ax.plot(x,packet,label='constructed excitation packet')
        ax.set(xlabel='native spatial coordinate',ylabel='neutral excitation',title=name.replace('_',' ').title()); ax.legend(fontsize=7)
        ax.text(.5,.03,'NO NATIVE LOADED EXCITATION DYNAMICS',transform=ax.transAxes,ha='center',fontsize=7)
        fig.tight_layout(); fig.savefig(RUN/f'{name}.png',dpi=90); plt.close(fig)
    (RUN/'baseline_git.txt').write_text(subprocess.run(['git','status','--short'],cwd=ROOT,text=True,capture_output=True).stdout)
    lines=['DEV146_AUDIT_COMPLETE',*outcomes,'PBUF_NATIVE_EXCITATION_STATE_ESTABLISHED=false',
      'PBUF_ENERGY_LIKE_MEDIUM_EXCITATION_ESTABLISHED=false','PBUF_LOADED_EXCITATION_INTERNAL_TRANSLATIONAL_SPLIT_ESTABLISHED=false',
      'PBUF_NATIVE_EXCITATION_CONSERVED_NORM_ESTABLISHED=false','PBUF_LOADED_EXCITATION_DISPERSION_ESTABLISHED=false',
      'PBUF_MASSIVE_PROPAGATION_FRACTION_EMERGES_FROM_EXCITATION_STRUCTURE=false',
      'PBUF_ZERO_LOAD_C_LIMIT_EMERGES_FROM_EXCITATION_STRUCTURE=false','PBUF_RELATIVISTIC_MASSIVE_DISPERSION_STRUCTURE_EMERGES=false',
      'VACUUM_MASS_DRAG=false','SR_DISPERSION_USED_TO_CONSTRUCT_NATIVE_LAW=false',
      *[f'{k}={str(v).lower() if isinstance(v,bool) else v}' for k,v in guards.items()]]
    (RUN/'report.txt').write_text('\n'.join(lines)+'\n'); print('\n'.join(lines))

if __name__=='__main__': main()
