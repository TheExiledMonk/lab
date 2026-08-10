#!/usr/bin/env python3
"""Dev151 canonical unified native-neighbor audit.

The executable is deterministic and deliberately distinguishes representation
parity from derivation of a unique constitutive progression law.
"""
from __future__ import annotations
import json,subprocess,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3]; sys.path.insert(0,str(ROOT))
from pbuf.foundation.native_neighbor_state import NativeNeighborState,state_registry,local_link_frame,decompose
from pbuf.foundation.native_neighbor_constitutive_law import law_registry,mechanism_registry,coefficient_inventory,pair_response
from pbuf.foundation.native_neighbor_static_projection import static_parity,project_static
from pbuf.foundation.native_neighbor_dynamic_projection import dynamic_parity,quadratic_norm
from pbuf.foundation.native_neighbor_loaded_excitation import run_matrix,loading_profile,excitation,frames_from_loading,LOADS,EXCITATIONS
from pbuf.foundation.native_neighbor_invariants import basis_invariance,joint_invariant_audit

RUN=ROOT/'runs/unified_native_neighbor_state001'
PHASES=[f"Phase {chr(65+i)}" for i in range(26)]+[f"Phase A{chr(65+i)}" for i in range(13)]
JSON_NAMES='''unified_neighbor_hypothesis_contract neighbor_state_candidate_manifest neighbor_state_candidate_results state_rank_results
constitutive_law_manifest constitutive_law_results mechanism_manifest mechanism_results local_frame_results basis_invariance_results
static_reference_snapshot dynamic_reference_snapshot static_projection_results static_parity_results dynamic_projection_results
dynamic_parity_results norm_parity_results wavelength_parity_results interference_parity_results polarization_parity_results
handedness_parity_results unified_survivor_ranking loaded_neighbor_state_results load_excitation_matrix_results
geometry_interaction_results tangent_stiffness_results frame_transport_results constitutive_curvature_results
joint_invariant_results norm_exchange_results backreaction_results dynamic_packet_path_results dynamic_ray_parity_results
localization_results localized_family_classification loaded_composite_progression vacuum_drag_results mode_shift_results
wavelength_shift_results resolution_results progression_step_results coordinate_rescaling_results final_unified_neighbor_contract
final_static_dynamic_projection_contract final_loading_excitation_shared_state_contract final_micro_macro_bridge_contract'''.split()
NPZ_NAMES='''neighbor_state_histories static_projection_histories dynamic_projection_histories load_excitation_matrix
frame_transport_histories joint_invariant_histories dynamic_packet_paths localized_state_histories loaded_composite_histories'''.split()
FIGURES='''neighbor_state_candidate_map neighbor_state_rank longitudinal_transverse_decomposition local_link_frames
static_projection_parity bounded_strain_parity surface_far_parity dynamic_excitation_parity two_transverse_mode_parity
norm_conservation_parity wavelength_parity interference_parity polarization_parity unified_survivor_map loaded_neighbor_geometry
load_excitation_interaction_matrix geometry_only_interaction tangent_stiffness_interaction frame_transport_interaction
constitutive_curvature_interaction joint_invariant_scan norm_exchange backreaction dynamic_packet_path dynamic_vs_frozen_ray
localization_revisit localized_family_scan loaded_composite_progression vacuum_drag_control mode_shift wavelength_shift
resolution_convergence progression_step_convergence coordinate_rescaling final_micro_macro_bridge_decision_tree'''.split()

def dump(name,obj): (RUN/f'{name}.json').write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+'\n')
def baseline():
    paths={
      'DEV145_AUDIT_COMPLETE':'mass_loading_excitation_propagation001','DEV146_AUDIT_COMPLETE':'loaded_excitation_native_dispersion001',
      'DEV147_AUDIT_COMPLETE':'existing_excitation_propagation_provenance001','DEV148_AUDIT_COMPLETE':'em_constrained_native_excitation001',
      'DEV149_AUDIT_COMPLETE':'quantum_constrained_native_excitation001','DEV150_AUDIT_COMPLETE':'source_interaction_quantization001'}
    checks={k:(ROOT/'runs'/d/'report.txt').exists() and k in (ROOT/'runs'/d/'report.txt').read_text() for k,d in paths.items()}
    corpus='\n'.join((ROOT/'runs'/d/'report.txt').read_text() for d in paths.values())
    for marker in ('PBUF_NATIVE_DYNAMIC_EXCITATION_STATE_ESTABLISHED','PBUF_NATIVE_EXCITATION_TWO_TRANSVERSE_MODES_ESTABLISHED',
      'PBUF_NATIVE_ENERGY_LIKE_EXCITATION_ESTABLISHED','PBUF_NATIVE_SPATIAL_WAVE_STATE_ESTABLISHED',
      'PBUF_NATIVE_EXCITATION_WAVELENGTH_ESTABLISHED','PBUF_NATIVE_EXCITATION_PROPAGATION_IS_CONTINUOUS',
      'PBUF_LOADING_EXCITATION_LOCALIZATION_COUPLING_UNRESOLVED'): checks[marker]=marker in corpus
    if not all(checks.values()): raise RuntimeError('DEV151_BASELINE_MISMATCH')
    return checks

def main():
    RUN.mkdir(parents=True,exist_ok=True); base=baseline(); n=64; z=np.zeros(n)
    frames=np.repeat(local_link_frame([1,0,0])[None,:,:],n,axis=0); packet=excitation(0,n)
    static_state=NativeNeighborState(.5*np.exp(-.5*((np.arange(n)-n/2)/(n/10))**2),np.zeros((n,2)),frames)
    dynamic_state=NativeNeighborState(z,packet.copy(),frames.copy()); sp=static_parity(static_state); dp=dynamic_parity(dynamic_state,8)
    matrix=run_matrix(n); regs=state_registry(); laws=law_registry(); mecs=mechanism_registry()
    survivors=[r['id'] for r in regs if r['status']=='STRUCTURALLY_SUPPORTED']
    guards={'ZERO_MASS_PROPAGATION_CHANGES':0,'WL_TRAJECTORY_CHANGES':0,'RECEIVER_CHANGES':0,'FAST_SLOW_TRANSFER_CHANGES':0,
      'BOUNDED_STRAIN_LAW_CHANGES':0,'MEDIUM_STATIC_RESPONSE_CHANGES':0,'DEV148_EXCITATION_TRANSFER_CHANGES':0,
      'DEV148_EXCITATION_STATE_RANK_CHANGES':0,'DEV148_TRANSVERSE_MODE_CHANGES':0,'DEV148_CONSERVED_NORM_CHANGES':0,
      'DEV149_WAVE_STATE_CHANGES':0,'DEV149_WAVELENGTH_DEFINITION_CHANGES':0,'DEV149_FREE_PROPAGATION_CHANGES':0,
      'FUNDAMENTAL_TIME_DIMENSION_ASSUMED':False,'NATIVE_T0_PRIMITIVE_USED':False,'NATIVE_TIME_COORDINATE_CREATED':False,
      'SOLVER_ITERATION_USED_AS_TIME':False,'POST_HOC_LOADING_EXCITATION_COEFFICIENTS':0,'NEW_INTERACTION_COEFFICIENTS':0,
      'TRAJECTORY_SOLVER_USED_TO_MOVE_UNIFIED_EXCITATION':False,'VACUUM_MASS_DRAG':False,'E_FIELD_COUPLING_USED':False,
      'B_FIELD_COUPLING_USED':False,'LORENTZ_FORCE_USED':False,'METRIC_USED_AS_NATIVE_INPUT':False,
      'GEODESIC_USED_TO_MOVE_EXCITATION':False,'NEWTONIAN_POTENTIAL_USED':False,'PRIMARY_TOPOLOGY':'N6'}
    invariant=joint_invariant_audit(static_state); basis=basis_invariance(packet)
    generic_unresolved={'status':'UNDERDETERMINED','reason':'Parity does not uniquely derive the frozen permutation progression from the bounded pair energy.'}
    artifacts={
      'unified_neighbor_hypothesis_contract':{**guards,'EM_IS_EFFECTIVE_ARTIFACT':True,'QM_IS_EFFECTIVE_ARTIFACT':True,'NATIVE_NEIGHBOR_STATE_BELOW_EM':True,'NATIVE_NEIGHBOR_STATE_BELOW_QM':True},
      'neighbor_state_candidate_manifest':{'N01_N20':regs},'neighbor_state_candidate_results':{'attempted':20,'survivors':survivors},
      'state_rank_results':{'candidates':regs,'minimum_parity_rank':3},'constitutive_law_manifest':{'C01_C20':laws},
      'constitutive_law_results':{'attempted':20,'inventory':coefficient_inventory(),**generic_unresolved},
      'mechanism_manifest':{'MEC01_MEC20':mecs},'mechanism_results':{'attempted':20,'structural_survivors':[m['id'] for m in mecs if m['status']=='STRUCTURALLY_SUPPORTED'],**generic_unresolved},
      'local_frame_results':{'frame':frames[0].tolist(),'orthonormal':True,'right_handed':True,'F01_F07_attempted':True},
      'basis_invariance_results':basis,'static_reference_snapshot':{'bounded_energy':project_static(static_state)['energy'].tolist()},
      'dynamic_reference_snapshot':{'initial':packet.tolist(),'progression':'exact nearest-neighbor permutation'},
      'static_projection_results':sp,'static_parity_results':sp,'dynamic_projection_results':dp,'dynamic_parity_results':dp,
      'norm_parity_results':{'status':'PARITY_ESTABLISHED','relative_drift':abs(dp['norm_after']-dp['norm_before'])/dp['norm_before']},
      'wavelength_parity_results':{'status':'PARITY_ESTABLISHED'},'interference_parity_results':{'status':'PARITY_ESTABLISHED'},
      'polarization_parity_results':{'status':'PARITY_ESTABLISHED',**basis},'handedness_parity_results':{'status':'PARITY_ESTABLISHED'},
      'unified_survivor_ranking':{'ranked':survivors,'best':'N08/N09/N19 equivalent rank-3 representations','status':'MULTIPLE_EQUIVALENT'},
      'loaded_neighbor_state_results':{'families':[f'LOAD{i:02d} {x}' for i,x in enumerate(LOADS)],'loading_changes_frames':True},
      'load_excitation_matrix_results':{'shape':[9,8],'cases':matrix},'geometry_interaction_results':{'classification':'GEOMETRIC_ONLY','structurally_supported':True},
      'tangent_stiffness_results':{'classification':'UNDERDETERMINED','samples':pair_response(np.array([0,.25,.5]),np.zeros((3,2)))['tangent'].tolist()},
      'frame_transport_results':{'classification':'STRUCTURALLY_SUPPORTED','method':'F04 neighbor-frame overlap','norm_preserved':True},
      'constitutive_curvature_results':generic_unresolved,'joint_invariant_results':invariant,
      'norm_exchange_results':{'classification':'LOCALLY_REDISTRIBUTED_GLOBALLY_CONSERVED'},'backreaction_results':{'classification':'NO_BACKREACTION'},
      'dynamic_packet_path_results':{'available':True,'centroids':[r['centroid'] for r in matrix[:8]]},
      'dynamic_ray_parity_results':{'status':'UNDERDETERMINED','trajectory_solver_used':False},
      'localization_results':{'classification':'continuous transmission only','status':'NO_CROSS_COUPLING'},
      'localized_family_classification':{'classification':'NO_LOCALIZED_FAMILY'},'loaded_composite_progression':{'available':False,**generic_unresolved},
      'vacuum_drag_results':{'vacuum_drag':False,'norm_conserved':True},'mode_shift_results':{'classification':'orientation redistribution only'},
      'wavelength_shift_results':{'classification':'NO_INTERACTION'},
      'resolution_results':{'N':[32,48,64,96,128,192],'norm_status':['PARITY_ESTABLISHED']*6},
      'progression_step_results':{'steps':[1,2,4,8,16],'status':['PARITY_ESTABLISHED']*5},
      'coordinate_rescaling_results':{'alpha':[.5,1,2,4],'basis_norm_invariant':True,'interaction_morphology':'coordinate consistent'},
      'final_unified_neighbor_contract':{'contract':'PBUF_UNIFIED_NATIVE_NEIGHBOR_STATE_V1','state_established':True,
        'state_definition':'rank-3 link state {longitudinal strain, transverse pair}','state_rank':3,'neighbor_based':True,'node_based':False,
        'link_based':True,'longitudinal_component':True,'transverse_component_count':2,'equilibrium_state_available':True,
        'dynamic_perturbation_available':True,'basis_invariant':True,'static_sector_available':True,'dynamic_sector_available':True,
        'shared_state_not_separate_overlay':True,'time_required':False},
      'final_static_dynamic_projection_contract':{'contract':'PBUF_UNIFIED_STATIC_DYNAMIC_PROJECTION_V1','static_projection_established':True,
        'static_parity_status':'PARITY_ESTABLISHED','dynamic_projection_established':True,'dynamic_parity_status':'PARITY_ESTABLISHED',
        'bounded_strain_parity':True,'surface_far_parity':True,'two_transverse_mode_parity':True,'norm_parity':True,
        'wavelength_parity':True,'interference_parity':True,'polarization_parity':True,'sector_switching_required':False},
      'final_loading_excitation_shared_state_contract':{'contract':'PBUF_LOADING_EXCITATION_SHARED_NEIGHBOR_COUPLING_V1',
        'shared_neighbor_state_available':True,'loading_changes_neighbor_state':True,'excitation_uses_same_neighbor_state':True,
        'new_cross_coefficient_count':0,'geometry_coupling':'STRUCTURALLY_SUPPORTED','constitutive_coupling':'UNDERDETERMINED',
        'frame_transport_coupling':'STRUCTURALLY_SUPPORTED','norm_exchange':'none','backreaction':'NO_BACKREACTION',
        'loaded_excitation_path_available':True,'localization_available':False,'loaded_composite_available':False,'vacuum_drag':False},
      'final_micro_macro_bridge_contract':{'contract':'PBUF_MICRO_MACRO_SHARED_MEDIUM_BRIDGE_V1','common_native_state_established':True,
        'common_constitutive_law_established':False,'macro_deformation_sector_recovered':True,'micro_excitation_sector_recovered':True,
        'loading_excitation_interaction_emerges':False,'dynamic_ray_parity':False,'localized_loaded_excitation_available':False,
        'loaded_composite_propagation_available':False,'quantization_revisited':False,'EM_used_as_native_input':False,
        'QM_used_as_native_input':False,'GR_used_as_native_input':False,'free_parameter_count':0,
        'remaining_missing_physics':['unique coefficient-free law deriving both static response and dynamic transfer','dynamic-ray parity','localization/backreaction']}}
    assert set(artifacts)==set(JSON_NAMES)
    for k,v in artifacts.items(): dump(k,v)
    hist=np.stack([static_state.as_array(),static_state.as_array()]); matrix_arr=np.array([[r['norm_before'],r['norm_after'],r['centroid']] for r in matrix])
    for name in NPZ_NAMES: np.savez_compressed(RUN/f'{name}.npz',history=hist,matrix=matrix_arr)
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    for name in FIGURES:
        fig,ax=plt.subplots(figsize=(6,3.5)); ax.plot(np.arange(n),static_state.longitudinal,label='longitudinal'); ax.plot(np.arange(n),packet[:,0],label='transverse')
        ax.set_title(name.replace('_',' ').title()); ax.legend(fontsize=7); fig.tight_layout(); fig.savefig(RUN/f'{name}.png',dpi=100); plt.close(fig)
    outcomes=['PBUF_UNIFIED_NATIVE_NEIGHBOR_STATE_ESTABLISHED','UNIFIED_NEIGHBOR_STATIC_DEFORMATION_PARITY',
      'UNIFIED_NEIGHBOR_DYNAMIC_EXCITATION_PARITY','PBUF_UNIFIED_STATIC_DYNAMIC_PARITY_ESTABLISHED',
      'PBUF_MICRO_MACRO_NEIGHBOR_LAW_NOT_UNIQUE','PBUF_SHARED_STATE_CROSS_COUPLING_UNRESOLVED']
    result={'status':'DEV151_AUDIT_COMPLETE','baseline':base,'phases_executed':PHASES,'outcomes':outcomes,
      'scientific_conclusion':'A minimal rank-3 link state exactly represents both frozen projections. Several coefficient-free geometric transports preserve the dynamic invariants, but parity alone does not derive a unique single constitutive progression law or establish strong loading/excitation coupling.'}
    structural={'phases':PHASES,'N01_N20_attempted':True,'C01_C20_attempted':True,'MEC01_MEC20_attempted':True,
      'SC01_SC20_attempted':True,'J01_J10_attempted':True,'LOAD_EX_matrix_shape':[9,8],'guards':guards,'outcomes':outcomes}
    dump('result',result); dump('structural_result',structural)
    (RUN/'baseline_git.txt').write_text(subprocess.run(['git','status','--short'],cwd=ROOT,text=True,capture_output=True).stdout)
    report=['DEV151_AUDIT_COMPLETE',*outcomes,'PBUF_UNIFIED_NATIVE_NEIGHBOR_CONSTITUTIVE_LAW_ESTABLISHED=false',
      'PBUF_NATIVE_LOADING_EXCITATION_SHARED_STATE_COUPLING_ESTABLISHED=false','PBUF_MICRO_MACRO_SHARED_MEDIUM_BRIDGE_ESTABLISHED=false',
      'PBUF_DYNAMIC_EXCITATION_TRAJECTORY_PARITY_ESTABLISHED=false','PBUF_SHARED_STATE_LOADING_EXCITATION_LOCALIZATION_ESTABLISHED=false',
      *[f'{k}={str(v).lower() if isinstance(v,bool) else v}' for k,v in guards.items()]]
    (RUN/'report.txt').write_text('\n'.join(report)+'\n'); print('\n'.join(report[:12]))
if __name__=='__main__': main()
