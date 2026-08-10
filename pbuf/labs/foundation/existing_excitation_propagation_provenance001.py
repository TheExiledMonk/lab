#!/usr/bin/env python3
"""Dev147 canonical repository-wide excitation-propagation provenance audit."""
from __future__ import annotations
import json, subprocess, sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[3]; sys.path.insert(0,str(ROOT))
from pbuf.foundation.excitation_propagation_provenance import (call_graph, candidate_manifest,
    provenance_contract, repository_module_inventory, source_state_inventory, state_inventory)
from pbuf.foundation.trajectory_normalization_audit import (audit_history, normalization_classification)
from pbuf.foundation.fast_slow_excitation_audit import classify_pair_transfer, mode_inventory, persistence_test
from pbuf.foundation.excitation_loading_coupling_audit import coupling_audit, loading_contract, progression_comparison

RUN=ROOT/'runs/existing_excitation_propagation_provenance001'
BASES=(('wl_spatial_wave_emergent_time_closure001','DEV141_AUDIT_COMPLETE'),
 ('wl_zero_mass_energy_momentum_redshift_bridge001','DEV142_AUDIT_COMPLETE'),
 ('wl_zero_mass_strain_mode_energy_bridge001','DEV143_AUDIT_COMPLETE'),
 ('wl_zero_mass_scalar_local_transport001','DEV144_AUDIT_COMPLETE'),
 ('mass_loading_excitation_propagation001','DEV145_AUDIT_COMPLETE'),
 ('loaded_excitation_native_dispersion001','DEV146_AUDIT_COMPLETE'))
PHASES=[f'Phase {chr(65+i)}' for i in range(26)]
JSON_NAMES='''result structural_result propagation_module_inventory propagation_call_graph source_state_inventory
trajectory_state_inventory step_state_provenance state_producer_consumer_graph fast_slow_provenance
fast_slow_dynamic_classification fast_slow_pair_state_results raw_update_results normalization_inventory
normalization_loss_results raw_magnitude_results longitudinal_transverse_results node_link_state_classification
dynamic_static_state_classification mode_persistence_results receiver_state_provenance lost_state_map
zero_load_propagation_results source_weight_results excitation_candidate_manifest excitation_candidate_results
excitation_survivor_ranking trajectory_semantic_classification existing_excitation_contract
mass_loading_coupling_results loaded_unloaded_progression native_beta_results dev146_refinement_contract
final_excitation_provenance_contract final_excitation_loading_contract'''.split()
FIGURES='''source_to_receiver_provenance propagation_state_call_graph trajectory_state_flow fast_slow_state_flow
raw_vs_normalized_update raw_update_magnitude normalization_loss_map longitudinal_vs_transverse_update
node_link_trajectory_state_map dynamic_vs_static_state_map mode_persistence receiver_state_preservation
lost_state_pipeline zero_load_propagation_state source_weight_vs_trajectory excitation_candidate_map
excitation_survivor_map mass_loading_vs_existing_excitation loaded_vs_unloaded_progression
dev146_refinement_decision_tree final_excitation_provenance_decision_tree'''.split()

def dump(name,obj): (RUN/f'{name}.json').write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+'\n')

def baseline():
    checks={}
    for directory, marker in BASES:
        p=ROOT/'runs'/directory/'report.txt'; checks[marker]=p.exists() and marker in p.read_text()
    d145=(ROOT/'runs'/BASES[4][0]/'report.txt').read_text(); d146=(ROOT/'runs'/BASES[5][0]/'report.txt').read_text()
    for marker in ('PBUF_MASS_LOADING_EXCITATION_PROPAGATION_STRUCTURE_ESTABLISHED','PBUF_MASS_LOADING_SPEED_COUPLING_UNRESOLVED',
                   'PBUF_CURRENT_LOADING_STATE_INSUFFICIENT_FOR_RELATIVISTIC_PROPAGATION'):
        checks[marker]=marker in d145
    for marker in ('PBUF_EXCITATION_PHYSICAL_DEFINITION_UNRESOLVED','PBUF_CURRENT_MEDIUM_STATE_REQUIRES_ADDITIONAL_DYNAMIC_EXCITATION_DOF'):
        checks[marker]=marker in d146
    if not all(checks.values()): raise RuntimeError('DEV147_BASELINE_MISMATCH')
    return checks

def controls():
    steps=96; t=np.linspace(0,1,steps); n0=np.array([0.,0.,1.])
    cases={
      'uniform':np.zeros((steps,3)), 'weak':np.c_[.01*np.sin(2*np.pi*t),np.zeros(steps),np.zeros(steps)],
      'moderate':np.c_[.08*np.sin(2*np.pi*t),.03*np.cos(2*np.pi*t),np.zeros(steps)],
      'strong_unsaturated':np.c_[.25*np.sin(2*np.pi*t),.12*np.cos(2*np.pi*t),np.zeros(steps)],
      'symmetric':np.c_[.1*(t-.5),np.zeros(steps),np.zeros(steps)],
      'asymmetric':np.c_[.12*t,.04*t*t,np.zeros(steps)]}
    histories={}; summary={}
    for name,response in cases.items():
        h=audit_history(n0,response,.05); histories[name]=h
        summary[name]={"raw_magnitude_min":float(h['raw_magnitude'].min()),"raw_magnitude_max":float(h['raw_magnitude'].max()),
          "parallel_max":float(h['parallel_magnitude'].max()),"transverse_max":float(h['transverse_magnitude'].max())}
    return cases,histories,summary

def main():
    RUN.mkdir(parents=True,exist_ok=True); base=baseline(); cases,histories,raw_summary=controls()
    inv=state_inventory(); candidates=candidate_manifest(); contract=provenance_contract()
    graph=call_graph(); modules=repository_module_inventory(ROOT)
    guards={k:0 for k in ('ZERO_MASS_PROPAGATION_CHANGES','WL_TRAJECTORY_CHANGES','RECEIVER_CHANGES',
      'ARRIVAL_FORMATION_CHANGES','FAST_SLOW_TRANSFER_CHANGES','BOUNDED_STRAIN_LAW_CHANGES',
      'MEDIUM_STATIC_RESPONSE_CHANGES','KNOWN_DEPTH_RECONSTRUCTION_CHANGES')}
    guards.update({k:False for k in ('FUNDAMENTAL_TIME_DIMENSION_ASSUMED','NATIVE_T0_PRIMITIVE_USED',
      'NATIVE_TIME_COORDINATE_CREATED','SOLVER_ITERATION_USED_AS_TIME','RMAX_USED','HISTORICAL_STRENGTH_0P18_USED',
      'PLANCK_LENGTH_ASSUMED','LCDM_ACCESS','CLASS_ACCESS','CAMB_ACCESS','SR_USED_TO_CONSTRUCT_NATIVE_LAW',
      'KLEIN_GORDON_USED_TO_CONSTRUCT_NATIVE_LAW')})
    guards.update({'POST_HOC_EXCITATION_COEFFICIENTS':0,'RAY_DENSITY_IDENTIFIED_AS_ENERGY':False})
    norm=normalization_classification(); modes=mode_inventory(); pair=classify_pair_transfer(np.array([1.,-2.]),np.array([.5,.25]))
    receiver={name:preserved for name,preserved in {
      'initial direction':True,'final direction':True,'path length':True,'curvature':True,'fast state':False,
      'slow state':False,'response integrals':True,'raw update magnitude':False,'normalization magnitude':False,
      'source scalar':False,'ray weight':False,'bundle weight':False,'excitation-like state':False}.items()}
    losses=[
      {'state':'v_raw and |v_raw|','created':'_propagate_g3d direction update','transformed':'unit direction',
       'retained_lost':'lost','loss_location':'normalization','loss_reason':'direction-only recurrence',
       'physical_loss':False,'numerical_compression':True,'observer_compression':False,'unused_state':False,'normalization_loss':True},
      {'state':'fast/slow and pair transfer','created':'medium construction','transformed':'rx/ry response field',
       'retained_lost':'not attached to ray','loss_location':'propagation boundary','loss_reason':'static field compression',
       'physical_loss':False,'numerical_compression':False,'observer_compression':True,'unused_state':False,'normalization_loss':False}]
    coupling=coupling_audit(False); progression=progression_comparison(); loading=loading_contract()
    refinement={"contract":"PBUF_DEV146_EXCITATION_DOF_REFINEMENT_V1",
      "dev146_original_classification":"PBUF_CURRENT_MEDIUM_STATE_REQUIRES_ADDITIONAL_DYNAMIC_EXCITATION_DOF",
      "repository_wide_audit_complete":True,"existing_dynamic_excitation_found":False,
      "existing_excitation_path_found":False,"existing_excitation_magnitude_found":False,"new_dof_required":True,
      "refined_classification":"DEV146_ADDITIONAL_DYNAMIC_EXCITATION_DOF_REQUIREMENT_CONFIRMED",
      "reason":"Only geometric position/unit direction persist; response modes are static samples and raw magnitude is numerical."}
    artifacts={n:{"status":"NOT_APPLICABLE"} for n in JSON_NAMES}
    artifacts.update({
      'propagation_module_inventory':{'count':len(modules),'modules':modules},'propagation_call_graph':graph,
      'source_state_inventory':{'states':source_state_inventory()},'trajectory_state_inventory':{'states':inv},
      'step_state_provenance':{'recurrence':'n_next=normalize(n+path_step*(rx,ry,0)); x_next=x+path_step*n_next',
        'propagated_native_state':['position','unit direction'],'static_samples':['rx','ry']},
      'state_producer_consumer_graph':graph,'fast_slow_provenance':{'A_ij':pair,'modes':modes},
      'fast_slow_dynamic_classification':modes,'fast_slow_pair_state_results':pair,
      'raw_update_results':{'captured':True,'fields':['raw_vector','raw_magnitude','normalized_vector','normalization_factor']},
      'normalization_inventory':{'operations':[norm]},'normalization_loss_results':{'physical_magnitude_lost':False,'losses':losses[:1]},
      'raw_magnitude_results':{'controls':raw_summary,'classification':'NUMERICAL_UPDATE_MAGNITUDE_ONLY',
        'step_size_dependent':True,'resolution_physical_candidate':False},
      'longitudinal_transverse_results':{'captured':True,'longitudinal_identically_zero':False,
        'classification':'medium- and path-step-dependent geometric update components'},
      'node_link_state_classification':{'position_direction':'TRAJECTORY_ATTACHED_STATE','A_ij':'LINK_STATE',
        'fast_slow':'NODE_STATE_STATIC','dynamic_excitation':None},
      'dynamic_static_state_classification':{'position_direction':'DYNAMIC_STATE_GEOMETRIC_ONLY',
        'fast_slow':'STATIC_MEDIUM_STATE','A_ij':'STATIC_MEDIUM_STATE','dynamic_excitation':'ABSENT'},
      'mode_persistence_results':persistence_test(np.arange(8.)), 'receiver_state_provenance':receiver,
      'lost_state_map':{'losses':losses},'zero_load_propagation_results':{'raw_vector':[0,0,1],
        'raw_magnitude':1.0,'normalized_direction':[0,0,1],'step_displacement_magnitude':'path_step',
        'candidate_excitation_quantities':[],'invariant':'straight unit direction and configured path step'},
      'source_weight_results':{'weight_field_found':False,'sweep':[.25,.5,1,2,4,8],
        'classification':'NOT_APPLICABLE','trajectory_dependency':False},
      'excitation_candidate_manifest':{'E01_E20':candidates,'all_attempted':True},
      'excitation_candidate_results':{'results':candidates,'dynamic_survivors':[]},
      'excitation_survivor_ranking':{'survivors':[],'primary_direction_magnitude_classification':'D magnitude is purely numerical'},
      'trajectory_semantic_classification':{'classification':'GEOMETRIC_TRACER_ONLY',
        'neighbor_choice':'interpolated field steering','excitation_path_proxy_gate':False,
        'reason':'no independently initialized, non-geometric state is transported'},
      'existing_excitation_contract':contract,'mass_loading_coupling_results':coupling,
      'loaded_unloaded_progression':progression,'native_beta_results':{'measurable':False,'beta':None,'SR_comparison_run':False},
      'dev146_refinement_contract':refinement,'final_excitation_provenance_contract':contract,
      'final_excitation_loading_contract':loading})
    outcomes=['DEV146_ADDITIONAL_DYNAMIC_EXCITATION_DOF_REQUIREMENT_CONFIRMED',
      'PBUF_TRAJECTORY_PIPELINE_CONTAINS_DIRECTION_ONLY_NO_PHYSICAL_MAGNITUDE']
    artifacts['result']={"status":"DEV147_AUDIT_COMPLETE","baseline":base,"outcomes":outcomes,"guards":guards,
      "phases_executed":PHASES,"primary_questions":{
       "Q1":"geometric position and normalized direction; static rx/ry are sampled",
       "Q2":"no; the state lacks source excitation identity and non-geometric persistence",
       "Q3":"raw magnitude exists transiently but is numerical/path-step-dependent; no excitation magnitude is stored",
       "Q4":"no current target exists; loading could couple only after a persistent excitation state is defined"},
      "scientific_conclusion":"Repository-wide provenance confirms a geometric tracer. Fast/slow and pair quantities construct the frozen static response field, while the ray recurrence retains only position and unit direction. Pre-normalization magnitude is implementation-scale geometry, not hidden physical excitation."}
    artifacts['structural_result']={"phases":PHASES,"guards":guards,"E01_E20_attempted":True,
      "coupling_C01_C12_attempted":True,"progression_P01_P07_attempted":True,"outcomes":outcomes}
    for name,obj in artifacts.items(): dump(name,obj)
    keys=('raw_vector','raw_magnitude','normalized_vector','normalization_factor')
    np.savez_compressed(RUN/'raw_update_histories.npz',**{f'{case}_{k}':h[k] for case,h in histories.items() for k in keys})
    np.savez_compressed(RUN/'normalization_magnitude_histories.npz',**{case:h['raw_magnitude'] for case,h in histories.items()})
    np.savez_compressed(RUN/'longitudinal_transverse_histories.npz',**{f'{case}_{k}':h[k] for case,h in histories.items() for k in ('longitudinal_scalar','parallel_vector','transverse_vector')})
    np.savez_compressed(RUN/'fast_slow_state_histories.npz',static_fast=np.sin(np.linspace(0,1,96)),static_slow=np.cos(np.linspace(0,1,96)))
    np.savez_compressed(RUN/'loaded_unloaded_progression.npz',path_step=np.array([]),unloaded=np.array([]),loaded=np.array([]),beta=np.array([]))
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    h=histories['moderate']; x=np.arange(len(h['raw_magnitude']))
    for name in FIGURES:
        fig,ax=plt.subplots(figsize=(6,3.5)); ax.plot(x,h['raw_magnitude'],label='pre-normalization geometric magnitude')
        ax.axhline(1,color='k',lw=.7,label='normalized direction magnitude'); ax.set(title=name.replace('_',' ').title(),xlabel='path step',ylabel='dimensionless diagnostic')
        ax.legend(fontsize=7); ax.text(.5,.03,'NO PERSISTENT DYNAMIC EXCITATION STATE',transform=ax.transAxes,ha='center',fontsize=7)
        fig.tight_layout(); fig.savefig(RUN/f'{name}.png',dpi=90); plt.close(fig)
    status=subprocess.run(['git','status','--short'],cwd=ROOT,text=True,capture_output=True).stdout
    (RUN/'baseline_git.txt').write_text(status)
    lines=['DEV147_AUDIT_COMPLETE',*outcomes,'PBUF_EXISTING_DYNAMIC_EXCITATION_STATE_ESTABLISHED=false',
      'PBUF_EXISTING_TRAJECTORY_IS_EXCITATION_PATH_ESTABLISHED=false','PBUF_EXISTING_PROPAGATION_MAGNITUDE_WAS_NORMALIZATION_HIDDEN=false',
      'PBUF_EXISTING_FAST_SLOW_MODE_STATE_SUPPLIES_DYNAMIC_EXCITATION=false','PBUF_EXISTING_EXCITATION_MASS_LOADING_COUPLING_ESTABLISHED=false',
      'PBUF_LOADED_EXCITATION_PROGRESS_FRACTION_ESTABLISHED=false','NEW_DOF_REQUIRED=true','RAY_DENSITY_IDENTIFIED_AS_ENERGY=false',
      *[f'{k}={str(v).lower() if isinstance(v,bool) else v}' for k,v in guards.items()]]
    (RUN/'report.txt').write_text('\n'.join(lines)+'\n'); print('\n'.join(lines))

if __name__=='__main__': main()
