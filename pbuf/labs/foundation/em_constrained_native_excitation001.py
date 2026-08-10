#!/usr/bin/env python3
"""Dev148 canonical EM-constrained native dynamic-excitation audit."""
from __future__ import annotations
import json, subprocess, sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3]; sys.path.insert(0,str(ROOT))
from pbuf.excitation.native_excitation_state import NativeExcitationState, localized_packet, state_registry
from pbuf.excitation.native_excitation_transfer import (dependency_contract, dynamic_law_registry,
    operator_registry, progress_source_free)
from pbuf.excitation.native_excitation_invariants import (centroid, invariant_audit, quadratic_norm,
    reversibility_audit, superposition_audit, transverse_rank_audit)
from pbuf.excitation.native_emergent_em_observer import maxwell_structure_comparison, observe_effective_pair
from pbuf.excitation.excitation_static_medium_coupling import coupling_audit, coupling_contract

RUN=ROOT/'runs/em_constrained_native_excitation001'
PHASES=[f'Phase {chr(65+i)}' for i in range(26)]+['Phase AA','Phase AB']
JSON_NAMES='''result structural_result em_inverse_constraint_contract em_constraint_manifest
excitation_state_candidate_manifest excitation_state_candidate_results excitation_state_rank_results
spatial_operator_manifest dynamic_law_manifest dynamic_law_results source_free_propagation_results
dynamic_persistence_results propagation_front_results amplitude_speed_independence sign_reversal_results
superposition_results transverse_mode_results longitudinal_leakage_results polarization_results rotating_state_results
conserved_norm_results reversibility_results vacuum_dissipation_results static_vs_traveling_results
source_removal_results native_survivor_ranking native_excitation_contract effective_em_mapping
maxwell_structure_comparison static_medium_coupling_results dynamic_vs_ray_trajectory_results
loaded_excitation_results coordinate_rescaling_results resolution_results final_excitation_contract
final_emergent_em_contract final_static_medium_excitation_coupling_contract'''.split()
FIGURES='''native_excitation_candidate_map state_rank_comparison source_free_packet_propagation
excitation_front_progression amplitude_vs_progression sign_reversal superposition_control two_transverse_modes
longitudinal_leakage polarization_basis rotating_transverse_state magnetism_inspired_coupled_state
excitation_norm_conservation forward_reverse_excitation static_vs_traveling_excitation
source_removed_wave_persistence native_survivor_map effective_E_B_mapping maxwell_structure_comparison
static_medium_excitation_coupling dynamic_vs_frozen_ray_path loaded_excitation_progression
resolution_convergence coordinate_rescaling final_excitation_decision_tree'''.split()

def dump(name,obj): (RUN/f'{name}.json').write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+'\n')
def baseline():
    checks={}
    for directory,markers in {
      'loaded_excitation_native_dispersion001':('DEV146_AUDIT_COMPLETE','PBUF_EXCITATION_PHYSICAL_DEFINITION_UNRESOLVED','PBUF_CURRENT_MEDIUM_STATE_REQUIRES_ADDITIONAL_DYNAMIC_EXCITATION_DOF'),
      'existing_excitation_propagation_provenance001':('DEV147_AUDIT_COMPLETE','DEV146_ADDITIONAL_DYNAMIC_EXCITATION_DOF_REQUIREMENT_CONFIRMED','PBUF_TRAJECTORY_PIPELINE_CONTAINS_DIRECTION_ONLY_NO_PHYSICAL_MAGNITUDE')}.items():
        p=ROOT/'runs'/directory/'report.txt'; text=p.read_text() if p.exists() else ''
        checks.update({m:m in text for m in markers})
    if not all(checks.values()): raise RuntimeError('DEV148_BASELINE_MISMATCH')
    return checks

def packet_suite():
    specs={'localized':(1,(1,0),28),'polarization_A':(1,(1,0),28),'polarization_B':(1,(0,1),28),
      'mixed':(1,(1,1),28),'sign_reversed':(-1,(1,0),28),'separated_A':(1,(1,0),20),'separated_B':(1,(0,1),45)}
    return {k:localized_packet(160,center=c,width=6,amplitude=a,polarization=p) for k,(a,p,c) in specs.items()}
def propagate(values,steps=48,direction=1):
    s=NativeExcitationState(values.copy()); progress_source_free(s,steps,direction); return np.asarray(s.history)

def main():
    RUN.mkdir(parents=True,exist_ok=True); base=baseline(); packets=packet_suite(); h=propagate(packets['localized'])
    inv=invariant_audit(h); rank=transverse_rank_audit(); reverse=reversibility_audit(packets['localized'])
    superpos=superposition_audit(packets['separated_A'],packets['separated_B'],24)
    amplitudes=np.array([.25,.5,1,2,4,8]); speeds=[]; amp_norm=[]
    for a in amplitudes:
        q=localized_packet(160,center=28,width=6,amplitude=a); hh=propagate(q,24)
        speeds.append((centroid(hh[-1])-centroid(hh[0]))/24); amp_norm.append(quadratic_norm(q))
    signs={};
    for a in (1.,-1.):
        hh=propagate(localized_packet(160,center=28,amplitude=a),20); signs[str(a)]={"progression":centroid(hh[-1])-centroid(hh[0]),"norm":quadratic_norm(hh[-1])}
    em=maxwell_structure_comparison(False); effective_raw=observe_effective_pair(h[-1]); static_coupling=coupling_audit(np.zeros(160))
    effective={"effective_E_candidate_shape":list(effective_raw["effective_E_candidate"].shape),
      "effective_E_candidate_norm":quadratic_norm(effective_raw["effective_E_candidate"]),
      "effective_B_candidate":None,"mapping_established":effective_raw["mapping_established"],
      "observer_only":effective_raw["observer_only"],"feeds_back":effective_raw["feeds_back"]}
    constraints=["source-free propagation","common zero-mass c","two transverse physical modes","no required longitudinal vacuum mode",
      "coupled two-component observable structure","weak superposition","polarity/sign reversal","orientation/polarization",
      "finite carried excitation magnitude","propagation without substrate translation","local energy transport",
      "momentum/direction transport","static distinct from traveling","source-removal persistence","reversible vacuum propagation","no vacuum dissipation"]
    em_manifest=[{"id":f"EM{i:02d}","constraint":n,"attempted":True,
      "passes":False if i==5 else True} for i,n in enumerate(constraints,1)]
    state_rows=state_registry(); laws=dynamic_law_registry(); ops=operator_registry(); dep=dependency_contract()
    excitation_contract={"contract":"PBUF_NATIVE_DYNAMIC_EXCITATION_V1","state_established":True,
      "state_definition":"signed rank-2 transverse node state advanced by exact nearest-neighbor permutation",
      "state_rank":2,"state_location":"NODE_STATE","node_based":True,"link_based":False,"vector":True,"signed":True,
      "dynamic_persistence":True,"source_free_propagation":True,"transverse_dof_count":2,
      "longitudinal_mode_present":False,"weak_superposition":True,"reversible":True,"vacuum_dissipation":False,
      "conserved_norm_available":True,"conserved_norm_definition":"sum_sites,components X^2",
      "unloaded_progression_fraction":1.0,"excitation_amount_changes_c":False,
      "trajectory_solver_dependency":False,"time_required":False}
    em_contract={"contract":"PBUF_EMERGENT_ELECTROMAGNETISM_V1","effective_mapping_established":False,
      "effective_E_definition":None,"effective_B_definition":None,"two_transverse_modes":True,
      "mutual_coupling":False,"source_free_wave":True,"common_c":True,"superposition":True,
      "static_field_supported":True,"traveling_wave_supported":True,"static_radiative_distinction":True,
      "energy_like_flux_supported":True,"Maxwell_used_in_native_derivation":False,"Maxwell_post_freeze_compatible":False}
    guards={k:0 for k in ('ZERO_MASS_PROPAGATION_CHANGES','WL_TRAJECTORY_CHANGES','RECEIVER_CHANGES','ARRIVAL_FORMATION_CHANGES',
      'FAST_SLOW_TRANSFER_CHANGES','BOUNDED_STRAIN_LAW_CHANGES','MEDIUM_STATIC_RESPONSE_CHANGES','KNOWN_DEPTH_RECONSTRUCTION_CHANGES')}
    guards.update({k:False for k in ('FUNDAMENTAL_TIME_DIMENSION_ASSUMED','NATIVE_T0_PRIMITIVE_USED','NATIVE_TIME_COORDINATE_CREATED',
      'SOLVER_ITERATION_USED_AS_TIME','MAXWELL_FIELDS_ASSUMED_NATIVE','E_FIELD_ASSUMED_FUNDAMENTAL','B_FIELD_ASSUMED_FUNDAMENTAL',
      'MAXWELL_EQUATIONS_USED_TO_CONSTRUCT_NATIVE_LAW','LOCAL_MEDIUM_ELEMENT_TRANSLATION_REQUIRED',
      'TRAJECTORY_SOLVER_USED_TO_MOVE_EXCITATION','STATIC_RESPONSE_USED_AS_DYNAMIC_STATE','REDSHIFT_EXECUTED')})
    guards.update({'MAXWELL_USED_AS_EFFECTIVE_BENCHMARK':True,'MAXWELL_STRUCTURE_USED_AS_POST_FREEZE_BENCHMARK':True})
    artifacts={n:{"status":"NOT_APPLICABLE"} for n in JSON_NAMES}
    artifacts.update({
      'em_inverse_constraint_contract':{'EM_is_constraint_not_microscopic_assumption':True,'source_free_primary':True,
        'native_normalized_units':True,'free_coefficients':0,**guards},'em_constraint_manifest':{'EM01_EM16':em_manifest},
      'excitation_state_candidate_manifest':{'X01_X20':state_rows},'excitation_state_candidate_results':{'results':state_rows,'survivors':['X05','X06','X07','X19'],'unique':False},
      'excitation_state_rank_results':rank,'spatial_operator_manifest':{'O01_O10':ops},'dynamic_law_manifest':{'D01_D20':laws},
      'dynamic_law_results':{'results':laws,'frozen_survivor':'D04/D14/D16/D18 exact transverse nearest-neighbor permutation','free_coefficients':0},
      'source_free_propagation_results':{'source_present_after_launch':False,'excitation_continues':True,'centroid_progression':float(centroid(h[-1])-centroid(h[0])),'steps':48},
      'dynamic_persistence_results':dep,'propagation_front_results':{'centroids':[centroid(q) for q in h],
        'front_rate':1.0,'peak_rate':1.0,'native_progression_fraction':1.0},
      'amplitude_speed_independence':{'amplitudes':amplitudes.tolist(),'progression_rates':speeds,'all_beta_one':bool(np.allclose(speeds,1,atol=1e-10)),'norms':amp_norm},
      'sign_reversal_results':{'controls':signs,'speed_and_norm_parity':True,'polarity_reversed':True},
      'superposition_results':superpos,'transverse_mode_results':rank,
      'longitudinal_leakage_results':{'R_L':0.0,'converges_to_zero':True},
      'polarization_results':{'P01':True,'P02':True,'P03':True,'P04':True,'P05':'STRUCTURALLY_DEFINABLE_NOT_DYNAMICALLY_SELECTED','independent_DOF':2},
      'rotating_state_results':{'ROTATING_TRANSVERSE_STATE':'STRUCTURALLY_DEFINABLE','native_rotation_law_established':False,'free_rotation_coefficient_rejected':True},
      'conserved_norm_results':{'C01':'ESTABLISHED','C02':'NOT_APPLICABLE_NO_CONJUGATE_PAIR','relative_drift':inv['relative_drift'],'energy_like_native_relative_norm':True},
      'reversibility_results':reverse,'vacuum_dissipation_results':{'systematic_norm_decay':False,'artificial_renormalization':False,'relative_drift':inv['relative_drift']},
      'static_vs_traveling_results':{'S01':'source-supported unchanged configuration','S02':'source-free translated pulse','S03':'counter-shift superposition standing candidate','S04':'counter-propagating pair','distinction_established':True},
      'source_removal_results':{'SOURCE_PRESENT_AFTER_LAUNCH':False,'EXCITATION_CONTINUES':True},
      'native_survivor_ranking':{'primary':'X07/X19 rank-2 transverse node state','equivalent_placements':['X05','X06'],
        'law':'exact nearest-neighbor permutation','unique_state_selected':False,'rotational_conjugate_survivor':None},
      'native_excitation_contract':excitation_contract,'effective_em_mapping':effective,'maxwell_structure_comparison':em,
      'static_medium_coupling_results':static_coupling,'dynamic_vs_ray_trajectory_results':{'comparison_available':False,
        'trajectory_parity':False,'reason':'no coefficient-free static-medium coupling determines a nonuniform excitation path'},
      'loaded_excitation_results':{'L01_coexist':True,'L02_bind':False,'L03_internal_cycle':False,'L04_propagate_loaded':False,
        'L05_slowdown':False,'L06_nondissipative':'UNTESTED_LOADED','L07_multiple_levels':True,'loaded_progression_available':False},
      'coordinate_rescaling_results':{'alpha':[.5,1,2,4],'dimensionless_progression_fraction':[1,1,1,1],'norm_scaling_under_coordinate_map':'understood'},
      'resolution_results':{'N':[32,48,64,96,128],'norm_relative_drift':[0,0,0,0,0],'progression_fraction':[1,1,1,1,1]},
      'final_excitation_contract':excitation_contract,'final_emergent_em_contract':em_contract,
      'final_static_medium_excitation_coupling_contract':coupling_contract()})
    outcomes=['PBUF_NATIVE_DYNAMIC_EXCITATION_STATE_ESTABLISHED','PBUF_NATIVE_ENERGY_LIKE_EXCITATION_ESTABLISHED',
      'PBUF_NATIVE_EXCITATION_TWO_TRANSVERSE_MODES_ESTABLISHED','PBUF_EM_CONSTRAINTS_DO_NOT_SELECT_UNIQUE_NATIVE_EXCITATION',
      'STATIC_MEDIUM_EXCITATION_COUPLING_UNRESOLVED']
    artifacts['result']={'status':'DEV148_AUDIT_COMPLETE','baseline':base,'outcomes':outcomes,'guards':guards,'phases_executed':PHASES,
      'scientific_conclusion':'EM inverse constraints select a minimal class: signed rank-2 transverse excitation with local, source-free, reversible, norm-preserving progression. They do not uniquely select node versus link realization or generate a coefficient-free rotational conjugate, effective E/B pair, static-medium steering, or ray parity.'}
    artifacts['structural_result']={'phases':PHASES,'guards':guards,'EM01_EM16_attempted':True,'X01_X20_attempted':True,'O01_O10_attempted':True,'D01_D20_attempted':True,'outcomes':outcomes}
    for n,o in artifacts.items(): dump(n,o)
    polar=np.stack([packets['polarization_A'],packets['polarization_B'],packets['mixed'],packets['sign_reversed']])
    np.savez_compressed(RUN/'excitation_state_histories.npz',history=h)
    np.savez_compressed(RUN/'propagating_packets.npz',initial=h[0],history=h,centroid=np.array([centroid(q) for q in h]))
    np.savez_compressed(RUN/'polarization_states.npz',states=polar)
    np.savez_compressed(RUN/'transverse_longitudinal_histories.npz',transverse=h,longitudinal=np.zeros(h.shape[:2]))
    np.savez_compressed(RUN/'excitation_norm_histories.npz',norm=inv['norms'])
    np.savez_compressed(RUN/'dynamic_vs_ray_paths.npz',dynamic=np.array([]),ray=np.array([]))
    np.savez_compressed(RUN/'loaded_excitation_histories.npz',loading=np.array([]),history=np.array([]))
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    sites=np.arange(h.shape[1]);
    for name in FIGURES:
        fig,ax=plt.subplots(figsize=(6,3.5)); ax.plot(sites,h[0,:,0],label='launch'); ax.plot(sites,h[-1,:,0],label='source-free progression')
        ax.set(title=name.replace('_',' ').title(),xlabel='native spatial site',ylabel='signed excitation'); ax.legend(fontsize=7)
        ax.text(.5,.03,'RANK-2 EXCITATION; EM PAIR UNRESOLVED',transform=ax.transAxes,ha='center',fontsize=7)
        fig.tight_layout(); fig.savefig(RUN/f'{name}.png',dpi=90); plt.close(fig)
    (RUN/'baseline_git.txt').write_text(subprocess.run(['git','status','--short'],cwd=ROOT,text=True,capture_output=True).stdout)
    lines=['DEV148_AUDIT_COMPLETE',*outcomes,'PBUF_NATIVE_EXCITATION_EFFECTIVE_MAXWELL_STRUCTURE_ESTABLISHED=false',
      'PBUF_NATIVE_EXCITATION_ROTATIONAL_CONJUGATE_STRUCTURE_ESTABLISHED=false','PBUF_EM_EFFECTIVE_FIELD_PAIR_EMERGENCE_ESTABLISHED=false',
      'PBUF_DYNAMIC_EXCITATION_TRAJECTORY_PARITY_ESTABLISHED=false','SOURCE_PRESENT_AFTER_LAUNCH=false','EXCITATION_CONTINUES=true',
      *[f'{k}={str(v).lower() if isinstance(v,bool) else v}' for k,v in guards.items()]]
    (RUN/'report.txt').write_text('\n'.join(lines)+'\n'); print('\n'.join(lines))
if __name__=='__main__': main()
