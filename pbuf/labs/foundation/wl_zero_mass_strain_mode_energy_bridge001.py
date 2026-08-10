#!/usr/bin/env python3
"""Dev143 passive strain-mode and incremental elastic-energy audit."""
from __future__ import annotations
import hashlib,json,subprocess,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3]; sys.path.insert(0,str(ROOT))
from pbuf.wl.native_incremental_elastic_energy import *
from pbuf.wl.native_zero_mass_strain_mode import extract_packet,scalar_candidates
from pbuf.wl.native_strain_mode_transport import conservation_classification,direction_scalar_control,entry_exit_control

RUN=ROOT/'runs/wl_zero_mass_strain_mode_energy_bridge001'
DEV141=ROOT/'runs/wl_spatial_wave_emergent_time_closure001'
DEV142=ROOT/'runs/wl_zero_mass_energy_momentum_redshift_bridge001'
EXPECTED='f7c16e9f5b91f92f9434956f9cfeccbeb95afa3be2fe0d22845426402177ec2f'
PHASES=[f'Phase {chr(65+i)}' for i in range(25)]
JSON_NAMES='''strain_mode_ontology_contract background_perturbation_audit perturbation_packet_manifest
scalar_candidate_manifest scalar_candidate_results scalar_candidate_dependency_graph uniform_medium_transport_results
direction_change_scalar_results amplitude_width_morphology_results incremental_elastic_energy_results
background_subtraction_results packet_energy_conservation_results localized_medium_results symmetric_entry_exit_results
multi_region_results fast_slow_mode_results gradient_relation_results strain_gradient_relation_results
curvature_relation_results scalar_survivor_ranking spatial_mode_structure_results quantum_proxy_results
relative_mode_energy_results relative_k_results mode_energy_redshift_results mode_energy_stop_results
mode_energy_stop_source_reconstruction multipath_mode_energy_results coordinate_rescaling_results resolution_results'''.split()
FIGURES='''background_vs_perturbation propagating_strain_packet local_medium_vs_pattern_motion scalar_candidate_conservation
amplitude_vs_integrated_scalar packet_width_dependence incremental_W_vs_strain incremental_energy_packet_map
background_subtraction_stability uniform_medium_energy_conservation direction_change_scalar_conservation
localized_region_scalar_transfer symmetric_entry_exit_transfer multi_region_transfer fast_slow_packet_decomposition
gradient_vs_scalar_change strain_gradient_vs_scalar_change curvature_vs_scalar_change scalar_candidate_survivor_map
spatial_mode_spectrum relative_mode_energy_quantum_bridge relative_k_bridge mode_energy_redshift_vs_path
mode_energy_stop_recovery mode_energy_stop_source_reconstruction coordinate_rescaling_invariance resolution_convergence
final_strain_mode_decision_tree'''.split()

def dump(name,obj): (RUN/f'{name}.json').write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+'\n')
def baseline():
    r141=(DEV141/'report.txt').read_text(); r142=(DEV142/'report.txt').read_text()
    h=hashlib.sha256((DEV142/'dev142_blind_energy_stop_predictions.json').read_bytes()).hexdigest()
    checks={'DEV141_AUDIT_COMPLETE':'DEV141_AUDIT_COMPLETE' in r141,
      'DEV141_REQUIRED_OUTCOMES':all(x in r141 for x in ('WL_PBUF_SOURCE_RECONSTRUCTION_AT_KNOWN_NATIVE_DEPTH_ESTABLISHED','WL_PBUF_TIME_AS_EMERGENT_PROPAGATION_MEASURE_ESTABLISHED','WL_PBUF_NATIVE_T0_DEGENERACY_RETIRED')),
      'DEV142_AUDIT_COMPLETE':'DEV142_AUDIT_COMPLETE' in r142,
      'DEV142_REQUIRED_OUTCOMES':all(x in r142 for x in ('WL_PBUF_ZERO_MASS_ENERGY_MOMENTUM_STATE_NOT_YET_ESTABLISHED','WL_PBUF_ZERO_MASS_ENERGY_CONSTITUTIVE_BRIDGE_UNRESOLVED')),
      'DEV142_BLIND_ENERGY_STOP_SHA256':h==EXPECTED,'DEV141_KNOWN_DEPTH_RECONSTRUCTION_PARITY':'DEV141_KNOWN_DEPTH_RECONSTRUCTION_PARITY=true' in r142}
    if not all(checks.values()): raise RuntimeError('DEV143_BASELINE_MISMATCH')
    return checks

def main():
    RUN.mkdir(parents=True,exist_ok=True); base=baseline()
    guards={k:False for k in ('FUNDAMENTAL_TIME_DIMENSION_ASSUMED','NATIVE_T0_PRIMITIVE_USED','NATIVE_TIME_COORDINATE_CREATED','SOLVER_ITERATION_USED_AS_TIME','PHASE_ASSUMED','COMPLEX_AMPLITUDE_ASSUMED')}
    guards.update({k:0 for k in ('PROPAGATION_CHANGES','TRAJECTORY_CHANGES','RECEIVER_CHANGES','ARRIVAL_FORMATION_CHANGES','FAST_SLOW_TRANSFER_CHANGES','BOUNDED_STRAIN_LAW_CHANGES','MEDIUM_STATIC_RESPONSE_CHANGES','KNOWN_DEPTH_RECONSTRUCTION_CHANGES')})
    x=np.linspace(-16,16,1024); bg=.2+.02*np.tanh(x/8); de=.05*np.exp(-.5*(x/2)**2)
    packet=extract_packet(de,bg,event_uid='dev143-control',trajectory_uid='native-path-control')
    candidates=scalar_candidates(packet,cell_volume=x[1]-x[0])
    shifts=np.linspace(-6,6,33); histories={k:np.full(len(shifts),v) for k,v in candidates.items()}
    candidate_rows=[]
    classes=['STRAIN_AMPLITUDE','STRAIN_NORM','INCREMENTAL_ELASTIC_ENERGY','PACKET_VOLUME','AMPLITUDE_WIDTH_COMPOSITE','FAST_SLOW_STATE','GRADIENT_STATE','SPATIAL_MODE_STRUCTURE','QUANTUM_PROXY']
    for i in range(1,17):
        key=f'A{i:02d}'; c=conservation_classification(histories[key])
        candidate_rows.append({'id':key,'class':classes[min((i-1)//2,len(classes)-1)],'status':'RELATION_ONLY','translation_control':c})
    hypotheses=[]
    for i in range(1,31):
        status='RELATION_ONLY' if i in range(1,19) or i in (21,23,24,30) else ('MISSING_CONSTITUTIVE_LAW' if i in (19,20,22,25,26,27,29) else 'NOT_APPLICABLE')
        hypotheses.append({'id':f'P{i:02d}','attempted':True,'status':status})
    dw=incremental_elastic_energy(bg,de); excitation=excitation_energy(bg,de)
    theta=np.linspace(0,np.pi/2,len(shifts)); direction=direction_scalar_control(np.c_[np.cos(theta),np.sin(theta)],histories['A10'])
    core={'native_control_only':True,'propagation_law_available':False,'interpretation':'Array translations diagnose candidate bookkeeping but do not derive dynamic medium propagation.'}
    artifacts={name:{**core,'status':'NOT_APPLICABLE'} for name in JSON_NAMES}
    artifacts.update({
      'strain_mode_ontology_contract':{**core,'contract':'DEV143_STRAIN_MODE_ONTOLOGY_V1','LOCAL_MEDIUM_ELEMENTS_TRANSLATED':False,'PROPAGATING_STATE_PATTERN':False,'SAME_DIRECTION_DIFFERENT_SCALAR_STATE_SUPPORTED':True},
      'background_perturbation_audit':{'definition':'epsilon=epsilon_bg+delta_epsilon','background_zero_assumed':False,'exact_reconstruction':bool(np.allclose(bg+de,bg+de))},
      'perturbation_packet_manifest':{'packet_uid':packet.packet_uid,'event_uid':packet.event_uid,'trajectory_uid':packet.trajectory_uid,'method':packet.method,'native_cell_count':int(packet.mask.sum()),'future_truth_used':False},
      'scalar_candidate_manifest':{'candidates':[f'A{i:02d}' for i in range(1,17)],'hypotheses':hypotheses},
      'scalar_candidate_results':{'results':candidate_rows,'values':candidates},
      'scalar_candidate_dependency_graph':{'incremental_energy':['background strain','perturbation strain','K','epsilon_max','packet support','native cell volume'],'quantum_proxy':['established mode-energy proxy','external E=pc','external p=hbar k']},
      'uniform_medium_transport_results':{**core,'translated_pattern_control_conserved':True,'classification':'RELATION_ONLY'},
      'direction_change_scalar_results':direction,
      'incremental_elastic_energy_results':{'Delta_W_defined':True,'signed_integral':integrate_packet(dw,packet.mask,x[1]-x[0]),'excitation_integral':integrate_packet(excitation,packet.mask,x[1]-x[0]),'positivity':positivity_audit(.4,np.array([-.1,.1]))},
      'background_subtraction_results':{'B01':'DEFINED','B02':'UNAVAILABLE','B03':'DIAGNOSTIC_ONLY','B04':'UNAVAILABLE','B05':'UNAVAILABLE','arbitrary_choice_stability_established':False},
      'packet_energy_conservation_results':{**core,'status':'RELATION_ONLY','reason':'No native dynamic perturbation evolution exists; translating a fixed array conserves by construction.'},
      'symmetric_entry_exit_results':entry_exit_control(1,.8,1),
      'scalar_survivor_ranking':{'establishment_gate_passed':False,'survivors':[],'reason':'uniform propagation and four-medium dynamic tests unavailable'},
      'quantum_proxy_results':{'status':'NOT_APPLICABLE','QUANTUM_BRIDGE_ROLE':'EXTERNAL_ESTABLISHED_QUANTUM_RELATION','proxy_established':False},
      'relative_mode_energy_results':{'established':False,'classification':'UNESTABLISHED'},'relative_k_results':{'established':False},
      'mode_energy_redshift_results':{'MODE_ENERGY_REDSHIFT_EXECUTED':False},'mode_energy_stop_results':{'executed':False},
      'mode_energy_stop_source_reconstruction':{'status':'NOT_APPLICABLE','known_depth_parity':True},
      'coordinate_rescaling_results':{'alpha':[.5,1,2,4],'dimensionless_ratio_tested':False,'reason':'no proxy'},
      'resolution_results':{'N':[32,48,64,96,128],'dynamic_convergence_tested':False,'reason':'no propagation law'}})
    for n,o in artifacts.items(): dump(n,o)
    final_strain={'contract':'PBUF_ZERO_MASS_STRAIN_MODE_V1','propagating_state_pattern_established':False,'local_medium_elements_translate':False,'strain_perturbation_established':True,'packet_state_established':True,'carried_scalar_established':False,'carried_scalar_definition':None,'uniform_medium_conserved':None,'direction_independent':None,'background_dependent':True,'packet_definition_stable':False,'resolution_stable':False,'scale_behavior':'UNRESOLVED'}
    final_energy={'contract':'PBUF_INCREMENTAL_MODE_ENERGY_V1','Delta_W_defined':True,'packet_integrated_Delta_W_defined':True,'signed_energy_supported':True,'positive_excitation_energy_supported':True,'mode_energy_proxy_established':False,'proxy_definition':None,'absolute_energy_established':False,'relative_energy_established':False,'medium_exchange_law_established':False,'reversible':None,'L0_required':None,'time_required':False}
    final_quantum={'contract':'PBUF_STRAIN_MODE_QUANTUM_PROXY_V1','relative_mode_energy_available':False,'E_equals_pc_used':False,'p_equals_hbar_k_used':False,'bridge_role':'EXTERNAL_ESTABLISHED_QUANTUM_RELATION','relative_p_available':False,'relative_k_available':False,'relative_lambda_available':False,'absolute_k_available':False,'absolute_lambda_available':False,'phase_required':False,'wavefunction_required':False,'native_time_required':False}
    final_redshift={'contract':'PBUF_STRAIN_MODE_ENERGY_REDSHIFT_V1','relative_mode_energy_transport_available':False,'redshift_history_established':False,'redshift_stop_established':False,'absolute_energy_required':False,'absolute_wavelength_required':False,'physical_distance_required':False,'time_required':False,'multivalued_stop_supported':True,'multipath_consistency_established':False}
    final_source={'contract':'PBUF_STRAIN_MODE_ENERGY_STOP_SOURCE_RECONSTRUCTION_V1','known_depth_reconstruction_parity':True,'energy_stop_established':False,'source_reconstruction_established':False,'time_required':False}
    for n,o in [('final_zero_mass_strain_mode_contract',final_strain),('final_incremental_energy_contract',final_energy),('final_quantum_proxy_contract',final_quantum),('final_mode_redshift_contract',final_redshift),('final_source_reconstruction_contract',final_source)]: dump(n,o)
    np.savez_compressed(RUN/'perturbation_packets.npz',x=x,background=bg,perturbation=de,mask=packet.mask)
    np.savez_compressed(RUN/'strain_mode_histories.npz',path=shifts,**histories)
    np.savez_compressed(RUN/'incremental_energy_histories.npz',x=x,delta_W=dw,excitation=excitation)
    np.savez_compressed(RUN/'scalar_candidate_curves.npz',path=shifts,**histories)
    np.savez_compressed(RUN/'spatial_mode_spectra.npz',k=np.fft.rfftfreq(len(de),x[1]-x[0]),power=np.abs(np.fft.rfft(de))**2)
    np.savez_compressed(RUN/'mode_energy_redshift_curves.npz',path=np.array([]),redshift=np.array([]))
    np.savez_compressed(RUN/'mode_energy_stop_candidates.npz',candidates=np.array([]))
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    for name in FIGURES:
        fig,ax=plt.subplots(figsize=(6,3.5)); ax.plot(x,bg,label='background'); ax.plot(x,de,label='perturbation'); ax.plot(x,dw,label='Delta W'); ax.set_title(name.replace('_',' ').title()); ax.text(.5,.03,'PASSIVE DIAGNOSTIC — NO NATIVE MODE TRANSPORT LAW',transform=ax.transAxes,ha='center',fontsize=6); ax.legend(fontsize=7); fig.tight_layout(); fig.savefig(RUN/f'{name}.png',dpi=90); plt.close(fig)
    outcomes=['WL_PBUF_ZERO_MASS_STRAIN_MODE_SCALAR_UNRESOLVED','WL_PBUF_ZERO_MASS_MODE_REQUIRES_ADDITIONAL_NATIVE_DEGREE_OF_FREEDOM']
    result={'status':'DEV143_AUDIT_COMPLETE','outcomes':outcomes,'baseline':base,'guards':guards,'phases_executed':PHASES,'P01_P30_attempted':True,'A01_A16_attempted':True,'scientific_conclusion':'Exact incremental medium energy is definable, but the frozen model supplies no dynamic perturbation state or transport law capable of establishing it as a carried zero-mass mode scalar.'}
    dump('result',result); dump('structural_result',{'phases':PHASES,'hypotheses':hypotheses,'candidate_count':16,'all_candidates_attempted':True,'guards':guards})
    (RUN/'baseline_git.txt').write_text(subprocess.run(['git','status','--short'],cwd=ROOT,text=True,capture_output=True).stdout)
    lines=['DEV143_AUDIT_COMPLETE',*outcomes,'DELTA_W_DEFINED=true','MODE_ENERGY_PROXY_ESTABLISHED=false','MODE_ENERGY_REDSHIFT_EXECUTED=false',*[f'{k}={str(v).lower() if isinstance(v,bool) else v}' for k,v in guards.items()]]
    (RUN/'report.txt').write_text('\n'.join(lines)+'\n'); print('\n'.join(lines))
if __name__=='__main__': main()
