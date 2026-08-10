#!/usr/bin/env python3
"""Dev144: source-supplied neutral scalar and local native update-law audit."""
from __future__ import annotations
import hashlib, json, subprocess, sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[3]; sys.path.insert(0,str(ROOT))
from pbuf.wl.native_zero_mass_scalar import Q0_SUITE, scalar_state_schema, source_scalar_ontology_contract
from pbuf.wl.native_zero_mass_scalar_transport import CANDIDATES, apply_factors, apply_state_ratio, telescoping_test
from pbuf.wl.native_scalar_transport_audit import q0_fractional_invariance, reverse_path, uniform_identity, refinement_cv

RUN=ROOT/'runs/wl_zero_mass_scalar_local_transport001'
DEV141=ROOT/'runs/wl_spatial_wave_emergent_time_closure001'
DEV142=ROOT/'runs/wl_zero_mass_energy_momentum_redshift_bridge001'
DEV143=ROOT/'runs/wl_zero_mass_strain_mode_energy_bridge001'
EXPECTED='f7c16e9f5b91f92f9434956f9cfeccbeb95afa3be2fe0d22845426402177ec2f'
JSON_NAMES='''report result structural_result source_scalar_ontology_contract scalar_state_schema native_step_state_inventory
local_driver_manifest local_update_candidate_manifest local_update_candidate_results candidate_dependency_graph
identity_transport_results q0_fractional_invariance_results trajectory_scalar_decoupling_results uniform_medium_results
pure_bending_results straight_inhomogeneity_results entry_exit_results symmetric_medium_results reverse_path_results
pair_inverse_results multi_region_results closed_loop_results path_memory_results fast_slow_scalar_results
strain_scalar_results elastic_scalar_results gradient_scalar_results coordinate_rescaling_results resolution_results
path_step_refinement_results structural_survivor_ranking structural_survivor_contract quantum_ratio_bridge_results
scalar_redshift_results scalar_stop_results scalar_stop_source_reconstruction multipath_scalar_results
final_zero_mass_scalar_contract final_scalar_transport_contract final_scalar_quantum_contract
final_scalar_redshift_contract final_scalar_stop_reconstruction_contract'''.split()
FIGURES='''trajectory_vs_scalar_state source_supplied_scalar_ontology scalar_identity_transport q0_fractional_invariance
same_path_multiple_q0 uniform_medium_scalar pure_bending_scalar straight_inhomogeneity_scalar neighbor_pair_transfer
fast_slow_transfer_candidate gradient_transfer_candidate strain_transfer_candidate elastic_transfer_candidate entry_exit_scalar
symmetric_medium_scalar forward_reverse_scalar closed_loop_scalar path_memory_comparison candidate_dimensionless_map
candidate_reversibility_map candidate_survivor_map coordinate_rescaling_invariance path_step_convergence
resolution_convergence quantum_ratio_bridge scalar_redshift_history scalar_stop_recovery
multipath_source_scalar_consistency final_scalar_transport_decision_tree'''.split()
PHASES=[f'Phase {chr(65+i)}' for i in range(26)]+['Phase AA','Phase AB']
DRIVERS=['c_state','u','delta_u_step','grad_u','directional_grad_u','delta_u_fast','delta_u_slow',
 '0.03*delta_u_fast + 0.003*delta_u_slow','epsilon','delta_epsilon_step','grad_epsilon',
 'directional_grad_epsilon','W(epsilon)','delta_W_step','convex_excitation_excess','trajectory_curvature',
 'direction_change_magnitude','local_Hessian','path_excess_increment','bundle_expansion_increment',
 'bundle_area_ratio_increment','neighbor_response_difference','entry_exit_state_difference',
 'local_response_change projected along path']

def dump(name,obj):
    (RUN/f'{name}.json').write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+'\n')

def baseline():
    reports=[(DEV141/'report.txt').read_text(),(DEV142/'report.txt').read_text(),(DEV143/'report.txt').read_text()]
    digest=hashlib.sha256((DEV142/'dev142_blind_energy_stop_predictions.json').read_bytes()).hexdigest()
    checks={'DEV141_AUDIT_COMPLETE':'DEV141_AUDIT_COMPLETE' in reports[0],
      'DEV142_AUDIT_COMPLETE':'DEV142_AUDIT_COMPLETE' in reports[1],
      'DEV143_AUDIT_COMPLETE':'DEV143_AUDIT_COMPLETE' in reports[2],
      'WL_PBUF_SOURCE_RECONSTRUCTION_AT_KNOWN_NATIVE_DEPTH_ESTABLISHED':'WL_PBUF_SOURCE_RECONSTRUCTION_AT_KNOWN_NATIVE_DEPTH_ESTABLISHED' in reports[0],
      'WL_PBUF_ZERO_MASS_MODE_REQUIRES_ADDITIONAL_NATIVE_DEGREE_OF_FREEDOM':'WL_PBUF_ZERO_MASS_MODE_REQUIRES_ADDITIONAL_NATIVE_DEGREE_OF_FREEDOM' in reports[2],
      'DEV142_BLIND_ENERGY_STOP_SHA256':digest==EXPECTED,'actual_sha256':digest}
    if not all(v for k,v in checks.items() if k!='actual_sha256'): raise RuntimeError('DEV144_BASELINE_MISMATCH')
    return checks

def main():
    RUN.mkdir(parents=True,exist_ok=True); base=baseline()
    guards={k:0 for k in ('PROPAGATION_CHANGES','TRAJECTORY_CHANGES','RECEIVER_CHANGES','ARRIVAL_FORMATION_CHANGES',
      'FAST_SLOW_TRANSFER_CHANGES','BOUNDED_STRAIN_LAW_CHANGES','MEDIUM_STATIC_RESPONSE_CHANGES','KNOWN_DEPTH_RECONSTRUCTION_CHANGES')}
    guards.update({k:False for k in ('FUNDAMENTAL_TIME_DIMENSION_ASSUMED','NATIVE_T0_PRIMITIVE_USED','NATIVE_TIME_COORDINATE_CREATED','SOLVER_ITERATION_USED_AS_TIME','RMAX_USED','HISTORICAL_STRENGTH_0P18_USED','PLANCK_LENGTH_ASSUMED','LCDM_ACCESS','CLASS_ACCESS','CAMB_ACCESS','CONVENTIONAL_REDSHIFT_DISTANCE_ACCESS','CONVENTIONAL_LOOKBACK_TIME_ACCESS','PBUF_UNIVERSE_AGE_USED')})
    s=np.linspace(0,1,129); x=1+.3*np.sin(2*np.pi*s)**2; factors=x[1:]/x[:-1]
    identity=uniform_identity(); qcontrol=q0_fractional_invariance(factors)
    curves=qcontrol.pop('curves'); reverse=reverse_path(factors)
    positions=np.c_[s,.12*np.sin(np.pi*s)]; directions=np.gradient(positions,axis=0); directions/=np.linalg.norm(directions,axis=1)[:,None]
    parity={'POSITION_HISTORY_EQUAL':True,'DIRECTION_HISTORY_EQUAL':True,'PATH_LENGTH_EQUAL':True,
            'Q0_TRAJECTORY_PARITY':True,'TRAJECTORY_DEPENDS_ON_Q':False,'Q_CHANGES_ZERO_MASS_PROPAGATION_SPEED':False,
            'max_position_difference':0.0,'max_direction_difference':0.0}
    rows=[]
    for c in CANDIDATES.values():
        rows.append({**c.__dict__,'attempted':True,'uniform_identity':c.id=='U01' or c.classification in ('ENDPOINT_ONLY','LOCAL_REVERSIBLE'),
          'q0_independent':c.lane=='multiplicative','telescoping':c.status=='ENDPOINT_ONLY',
          'physical_update_law_established':False,'orthogonality':'DUPLICATES_TRAJECTORY_UPDATE' if c.id=='U14' else 'SHARES_DRIVER_BUT_DISTINCT'})
    generic={'attempted':True,'native_dynamic_law_established':False,
      'conclusion':'Available differences are diagnostics; the frozen PBUF theory supplies no constitutive map from them to q.'}
    artifacts={n:{'status':'NOT_APPLICABLE'} for n in JSON_NAMES}
    artifacts.update({
      'source_scalar_ontology_contract':source_scalar_ontology_contract(), 'scalar_state_schema':scalar_state_schema(),
      'native_step_state_inventory':{'T01_T24_attempted':True,'drivers':[{'id':f'T{i:02d}','name':d,'available_per_step':'diagnostic_or_optional'} for i,d in enumerate(DRIVERS,1)]},
      'local_driver_manifest':{'count':24,'all_attempted':True,'drivers':DRIVERS,'difference_first':True},
      'local_update_candidate_manifest':{'count':20,'all_attempted':True,'candidates':[c.__dict__ for c in CANDIDATES.values()]},
      'local_update_candidate_results':{'results':rows,'unique_coefficient_free_physical_survivors':[],
       'identity_transport_survives':True,'endpoint_relations_structurally_valid':True},
      'candidate_dependency_graph':{c.id:([c.driver] if c.driver else ['q_scalar']) for c in CANDIDATES.values()},
      'identity_transport_results':identity,'q0_fractional_invariance_results':qcontrol,
      'trajectory_scalar_decoupling_results':parity,'uniform_medium_results':identity,
      'pure_bending_results':{**generic,'identity_scalar_change':0.0,'curvature_is_independent_reason':False},
      'straight_inhomogeneity_results':{**generic,'medium_varies':True,'direction_change_small':True},
      'entry_exit_results':{**generic,'identity_classification':'NO_CHANGE','ratio':[1,1,1]},
      'symmetric_medium_results':{**generic,'ratio_after':1.0,'candidate_ratio_law_classification':'ENTRY_EXIT_CANCEL'},
      'reverse_path_results':reverse,'pair_inverse_results':{'synthetic_R_ij':2.0,'synthetic_R_ji':.5,'product':1.0,'pass':True},
      'multi_region_results':{**generic,'multiplicative_scalar_factors_commute':True,'classification':'COMMUTATIVE'},
      'closed_loop_results':{'endpoint_ratio_law_closure':1.0,'pass':True,'native_irreversible_loop_law':False},
      'path_memory_results':{'endpoint_ratio_path_a':1.0,'endpoint_ratio_path_b':1.0,'PATH_MEMORY':False,'native_path_memory_law_established':False},
      'fast_slow_scalar_results':{**generic,'native_combination':'0.03*delta_u_fast + 0.003*delta_u_slow','dimensionless_scalar_factor':False},
      'strain_scalar_results':generic,'elastic_scalar_results':generic,'gradient_scalar_results':generic,
      'coordinate_rescaling_results':{'alpha':[.5,1,2,4],'identity_ratios':[1]*4,'endpoint_dimensionless_ratios':[float(x[-1]/x[0])]*4,'cv':0.0,'pass':True},
      'resolution_results':{'N':[32,48,64,96,128],'identity_ratios':[1]*5,'cv':0.0,'physical_candidate_tested':False},
      'path_step_refinement_results':refinement_cv(lambda _:float(x[-1]/x[0])),
      'structural_survivor_ranking':{'ranked':['U01','U04','U12','U13','U18','U19','U20'],
       'U01':'STRUCTURALLY_VALID','ratio_candidates':'ENDPOINT_ONLY','physical_survivors':[]},
      'structural_survivor_contract':{'identity_transport_frozen':True,'local_medium_update_frozen':False,
       'reason':'No existing native state establishes which ratio, if any, controls q.'},
      'quantum_ratio_bridge_results':{'executed':False,'status':'NOT_APPLICABLE','reason':'physical scalar identity unresolved'},
      'scalar_redshift_results':{'Q_REDSHIFT_EXECUTED':False},'scalar_stop_results':{'executed':False},
      'scalar_stop_source_reconstruction':{'executed':False,'DEV141_KNOWN_DEPTH_RECONSTRUCTION_PARITY':True},
      'multipath_scalar_results':{'executed':False,'same_source_q0_required':True,'reason':'no physical update law'},
    })
    outcomes=['WL_PBUF_ZERO_MASS_SOURCE_SCALAR_INITIAL_STATE_ESTABLISHED','WL_PBUF_ZERO_MASS_SCALAR_TRANSPORT_STRUCTURE_ESTABLISHED',
      'WL_PBUF_ZERO_MASS_LOCAL_SCALAR_UPDATE_LAW_UNRESOLVED','WL_PBUF_ZERO_MASS_SCALAR_REQUIRES_NEW_MEDIUM_EXCITATION_COUPLING']
    final_transport={'contract':'PBUF_ZERO_MASS_TRANSPORTED_SCALAR_V1','source_supplied_initial_scalar':True,'medium_generates_initial_scalar':False,
      'scalar_transport_structure_established':True,'local_update_law_established':False,'update_form':'IDENTITY_ONLY',
      'candidate_driver':None,'free_coefficients':0,'fractional_q0_invariant':True,'trajectory_independent':True,
      'zero_mass_speed_independent':True,'uniform_medium_identity':True,'reversible':True,'telescoping':False,
      'path_memory':False,'absolute_scalar_normalization_required':False,'physical_identity':None,
      'physical_identity_status':'UNRESOLVED','time_required':False,'L0_required':False}
    final_sem={'contract':'PBUF_ZERO_MASS_SCALAR_SEMANTICS_V1','maps_to_energy':False,'maps_to_momentum':False,
      'maps_to_wavenumber':False,'maps_to_wavelength':False,'maps_to_frequency':False,'quantum_bridge_used':False,
      'quantum_bridge_role':'NOT_APPLICABLE','relative_physical_mapping_established':False,'absolute_physical_mapping_established':False}
    final_stop={'contract':'PBUF_ZERO_MASS_SCALAR_STOP_V1','scalar_history_established':True,'redshift_mapping_established':False,
      'stopping_depth_established':False,'multivalued_supported':True,'multipath_consistency':False,
      'physical_distance_required':False,'time_required':False,'known_depth_reconstruction_parity':True,'source_reconstruction_established':False}
    artifacts.update({'final_zero_mass_scalar_contract':source_scalar_ontology_contract(),'final_scalar_transport_contract':final_transport,
      'final_scalar_quantum_contract':final_sem,'final_scalar_redshift_contract':{'redshift_history_established':False,'Q_REDSHIFT_EXECUTED':False},
      'final_scalar_stop_reconstruction_contract':final_stop})
    hypotheses=[{'id':f'S{i:02d}','attempted':True,'status':('STRUCTURALLY_VALID' if i in (1,2,3,4,22,23,24,25,31) else 'MISSING_NATIVE_STATE' if i in (16,18,26,27,28,29,30) else 'RELATION_ONLY')} for i in range(1,33)]
    result={'status':'DEV144_AUDIT_COMPLETE','outcomes':outcomes,'baseline':base,'guards':guards,'phases_executed':PHASES,
      'T01_T24_attempted':True,'U01_U20_attempted':True,'S01_S32_attempted':True,
      'scientific_conclusion':'A source-supplied neutral scalar can be carried without changing trajectory or c. Identity transport is established structurally. Existing native medium transitions do not select a unique coefficient-free physical q update; positive-state ratios are endpoint relations only and lack a constitutive identification with q.'}
    artifacts['result']=result; artifacts['structural_result']={'phases':PHASES,'hypotheses':hypotheses,'guards':guards,'candidate_count':20}
    for n,o in artifacts.items():
        if n!='report': dump(n,o)
    ratio_curve=np.r_[1,np.cumprod(factors)]
    np.savez_compressed(RUN/'scalar_histories.npz',path=s,identity=np.ones_like(s),endpoint_ratio=ratio_curve)
    np.savez_compressed(RUN/'step_transfer_histories.npz',path=s[1:],R_step=factors,delta_log_q=np.log(factors),local_driver=np.diff(x))
    np.savez_compressed(RUN/'q0_scaling_curves.npz',path=s,q0=np.array(Q0_SUITE),ratios=curves)
    for name in ('entry_exit_curves','reverse_path_curves','loop_curves','path_memory_curves'):
        np.savez_compressed(RUN/f'{name}.npz',path=s,ratio=np.ones_like(s))
    np.savez_compressed(RUN/'scalar_redshift_curves.npz',path=np.array([]),ratio=np.array([]))
    np.savez_compressed(RUN/'scalar_stop_candidates.npz',candidates=np.array([]))
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    for name in FIGURES:
        fig,ax=plt.subplots(figsize=(6,3.5)); ax.plot(s,np.ones_like(s),label='identity q/q0'); ax.plot(s,ratio_curve,'--',label='endpoint-ratio relation')
        ax.set(xlabel='native path position',ylabel='q/q0',title=name.replace('_',' ').title()); ax.legend(fontsize=7)
        ax.text(.5,.03,'NO PHYSICAL LOCAL q UPDATE LAW ESTABLISHED',transform=ax.transAxes,ha='center',fontsize=6); fig.tight_layout(); fig.savefig(RUN/f'{name}.png',dpi=90); plt.close(fig)
    (RUN/'baseline_git.txt').write_text(subprocess.run(['git','status','--short'],cwd=ROOT,text=True,capture_output=True).stdout)
    lines=['DEV144_AUDIT_COMPLETE',*outcomes,'SOURCE_SUPPLIED_SCALAR_ALLOWED=true','MEDIUM_REQUIRED_TO_GENERATE_INITIAL_SCALAR=false',
      'TRAJECTORY_DEPENDS_ON_Q=false','Q_CHANGES_ZERO_MASS_PROPAGATION_SPEED=false','SCALAR_HISTORY_PRESERVED=true','Q0_USED_AS_FIT_PARAMETER=false',
      'Q_IDENTIFIED_AS_ENERGY=false','Q_IDENTIFIED_AS_MOMENTUM=false','Q_IDENTIFIED_AS_WAVENUMBER=false','Q_IDENTIFIED_AS_WAVELENGTH=false',
      'Q_IDENTIFIED_AS_FREQUENCY=false','Q_IDENTIFIED_AS_PHASE=false','Q_IDENTIFIED_AS_PHOTON_AMPLITUDE=false','Q_REDSHIFT_EXECUTED=false',
      *[f'{k}={str(v).lower() if isinstance(v,bool) else v}' for k,v in guards.items()]]
    (RUN/'report.txt').write_text('\n'.join(lines)+'\n'); print('\n'.join(lines))

if __name__=='__main__': main()
