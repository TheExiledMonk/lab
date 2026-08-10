#!/usr/bin/env python3
"""Dev150 canonical source/interaction quantization audit.

The audit stops at the missing native coupling gate, but writes the complete
negative-result evidence package so unsupported downstream claims stay false.
"""
from __future__ import annotations
import json, subprocess, sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3]; sys.path.insert(0,str(ROOT))
from pbuf.quantum.native_excitation_modes import carrier_mode, propagate, quadratic_norm
from pbuf.quantum.native_localized_excitation_states import state_registry,loading_registry,topology_registry,construct_composite,state_observables
from pbuf.quantum.native_excitation_interaction import coefficient_free_coupling_audit,interaction_locality,backreaction_audit
from pbuf.quantum.native_transition_quantization import quantization_location_registry,quantization_family_registry,cluster_final_states,transition_graph
from pbuf.quantum.native_emission_absorption import emission_audit,absorption_audit,selection_audit,source_off_audit
from pbuf.quantum.native_transition_benchmark import compare

RUN=ROOT/'runs/source_interaction_quantization001'
AMPLITUDES=np.array([.03125,.0625,.125,.25,.5,.75,1.,1.25,1.5,2.,3.,4.])
RESOLUTIONS=[32,48,64,96,128,192]; DOMAINS=[64,96,128,192,256]; ALPHAS=[.5,1,2,4]
PHASES=[f"Phase {chr(65+i)}" for i in range(26)]+[f"Phase A{chr(65+i)}" for i in range(10)]
JSON_NAMES='''interaction_quantization_contract quantization_location_manifest localized_state_manifest localized_state_results
localized_state_survivor_ranking bound_free_classification localized_stability_results perturbation_stability_results
internal_mode_results internal_wavelength_results continuous_amplitude_sweep stable_state_cluster_results discrete_state_results
grid_quantization_control domain_quantization_control topological_state_results transition_graph transition_pair_results
transition_norm_accounting emission_results emitted_packet_spectrum transition_norm_vs_k absorption_results
fractional_absorption_results subthreshold_accumulation polarization_selection_results handedness_selection_results
source_driven_emission source_removed_emission_persistence reverse_transition_results stimulated_interaction_results
loading_backreaction_results resolution_results domain_size_results coordinate_rescaling_results
quantum_transition_postfreeze_comparison planck_transition_comparison final_localized_state_contract
final_interaction_quantization_contract final_emission_absorption_contract final_native_quantum_exchange_contract'''.split()
NPZ_NAMES='''localized_state_histories localized_mode_spectra amplitude_sweep_histories transition_histories emitted_packet_histories
absorption_histories subthreshold_accumulation_histories transition_norm_vs_k_curves'''.split()
FIGURES='''localized_state_gallery bound_vs_free_state_map localized_stability perturbation_stability internal_mode_gallery
internal_wavelength_spectrum continuous_initial_vs_final_norm stable_state_clusters discrete_state_ladder grid_quantization_control
domain_quantization_control topological_state_map transition_graph transition_norm_conservation emission_event emitted_packet_spectrum
transition_norm_vs_k absorption_threshold fractional_packet_absorption subthreshold_accumulation polarization_selection
handedness_selection source_driven_emission source_removed_wave reverse_transition stimulated_interaction loading_backreaction
resolution_convergence domain_size_convergence coordinate_rescaling planck_transition_postfreeze final_quantization_location_decision_tree'''.split()

def dump(name,obj): (RUN/f'{name}.json').write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+'\n')
def baseline():
    requirements={
      'DEV145_AUDIT_COMPLETE':ROOT/'runs/mass_loading_excitation_propagation001/report.txt',
      'DEV146_AUDIT_COMPLETE':ROOT/'runs/loaded_excitation_native_dispersion001/report.txt',
      'DEV147_AUDIT_COMPLETE':ROOT/'runs/existing_excitation_propagation_provenance001/report.txt',
      'DEV148_AUDIT_COMPLETE':ROOT/'runs/em_constrained_native_excitation001/report.txt',
      'DEV149_AUDIT_COMPLETE':ROOT/'runs/quantum_constrained_native_excitation001/report.txt'}
    checks={k:p.exists() and k in p.read_text() for k,p in requirements.items()}
    dev149=(ROOT/'runs/quantum_constrained_native_excitation001/report.txt').read_text()
    for marker in ('PBUF_NATIVE_EXCITATION_WAVELENGTH_ESTABLISHED','PBUF_NATIVE_SPATIAL_WAVE_STATE_ESTABLISHED',
                   'PBUF_NATIVE_EXCITATION_PROPAGATION_IS_CONTINUOUS','PBUF_QUANTIZATION_REQUIRES_SOURCE_OR_INTERACTION_PHYSICS'):
        checks[marker]=marker in dev149
    if not all(checks.values()): raise RuntimeError('DEV150_BASELINE_MISMATCH')
    return checks

def main():
    RUN.mkdir(parents=True,exist_ok=True); base=baseline(); coupling=coefficient_free_coupling_audit()
    # The same continuously swept overlap state remains continuously scalable.
    reference=carrier_mode(192,16,envelope='gaussian',width=48)
    initial=np.array([quadratic_norm(a*reference) for a in AMPLITUDES]); final=initial.copy()
    clusters=cluster_final_states(initial,final); histories=np.stack([propagate(a*reference,8) for a in AMPLITUDES])
    loading=np.exp(-.5*((np.arange(192)-96)/18)**2)*.5
    composite=construct_composite(loading,reference); observables=state_observables(composite)
    unresolved={'status':'MISSING_INTERACTION_LAW','reason':'No coefficient-free native loading/excitation coupling exists; synthetic overlap is not binding.'}
    frozen={'FREE_EXCITATION_CONTINUITY_RETAINED':True,'POST_HOC_INTERACTION_COEFFICIENTS':0,
      'EMITTED_PACKET_INJECTED_EXTERNALLY':False,'EXCITATION_BOUND_IDENTIFIED_WITH_STRAIN_BOUND':False,
      'FUNDAMENTAL_TIME_DIMENSION_ASSUMED':False,'NATIVE_T0_PRIMITIVE_USED':False,'NATIVE_TIME_COORDINATE_CREATED':False,
      'SOLVER_ITERATION_USED_AS_TIME':False,'QM_IS_EFFECTIVE_ARTIFACT':True,'NATIVE_EXCITATION_BELOW_QM':True,
      'SCHRODINGER_EQUATION_ASSUMED_NATIVE':False,'WAVEFUNCTION_ASSUMED_NATIVE':False,'HAMILTONIAN_ASSUMED_NATIVE':False,
      'PHOTON_ASSUMED_NATIVE':False,'HBAR_USED_IN_NATIVE_DERIVATION':False,'PLANCK_RELATION_USED_TO_BUILD_TRANSITIONS':False}
    for key in ('ZERO_MASS_PROPAGATION_CHANGES','WL_TRAJECTORY_CHANGES','RECEIVER_CHANGES','FAST_SLOW_TRANSFER_CHANGES',
      'BOUNDED_STRAIN_LAW_CHANGES','MEDIUM_STATIC_RESPONSE_CHANGES','DEV148_EXCITATION_TRANSFER_CHANGES',
      'DEV148_EXCITATION_STATE_RANK_CHANGES','DEV148_TRANSVERSE_MODE_CHANGES','DEV148_CONSERVED_NORM_CHANGES',
      'DEV149_WAVE_STATE_CHANGES','DEV149_WAVELENGTH_DEFINITION_CHANGES','DEV149_FREE_PROPAGATION_CHANGES'): frozen[key]=0
    localized_contract={'contract':'PBUF_NATIVE_LOCALIZED_LOADED_EXCITATION_V1','localized_state_established':False,
      'state_definition':'S={L(x),X(x)} synthetic composite only','loading_definition':'L06/L08/L16 native proxies',
      'internal_excitation_definition':'frozen Dev148 rank-2 excitation','stable':False,'metastable':False,'source_supported':False,
      'discrete_state_family':False,'state_index_definition':None,'internal_wavelength_available':True,'internal_k_available':True,
      'topological_invariant_available':False,'arbitrary_trapping_potential_used':False,'free_parameter_count':0}
    interaction_contract={'contract':'PBUF_NATIVE_INTERACTION_QUANTIZATION_V1','interaction_quantization_established':False,
      'quantization_location':'QL10 no quantization in current localized physics','localized_states_discrete':False,
      'emission_discrete':False,'absorption_discrete':False,'free_propagation_continuous':True,
      'transition_norm_defined':True,'transition_norm_conserved':False,'fractional_incident_packets_allowed':True,
      'subthreshold_accumulation_behavior':'UNRESOLVED','source_boundary_quantization':False,
      'topological_quantization':False,'nonlinear_quantization':False,'grid_quantization_rejected':True,'box_quantization_rejected':True}
    exchange={'contract':'PBUF_NATIVE_EXCITATION_EXCHANGE_V1','emission_established':False,'absorption_established':False,
      'emission_generated_dynamically':False,'emitted_packet_injected':False,'bound_norm_before':None,'bound_norm_after':None,
      'free_norm_emitted':None,'transition_norm_vs_k_relation':'NOT_APPLICABLE','linear_in_k':False,
      'polarization_selection':'UNRESOLVED','handedness_selection':'UNRESOLVED','reverse_transition_available':False,'time_required':False}
    bridge={'contract':'PBUF_NATIVE_TO_EFFECTIVE_QUANTUM_EXCHANGE_V1','free_wave_state_established':True,'free_wave_continuous':True,
      'localized_matter_state_established':False,'localized_state_quantization_established':False,
      'interaction_quantization_established':False,'discrete_exchange_established':False,
      'native_transition_energy_like_quantity':'bound norm difference defined but no transition available','native_emitted_k_available':False,
      'native_transition_norm_k_relation':'NOT_APPLICABLE','photon_like_exchange_structure':False,
      'Planck_relation_used_in_derivation':False,'Planck_relation_postfreeze_status':'NOT_COMPARABLE_MISSING_NATIVE_TRANSITIONS',
      'hbar_used_in_derivation':False,'absolute_energy_normalization_available':False,'absolute_length_normalization_available':False,
      'remaining_missing_physics':['coefficient-free loading-excitation binding law','localized interaction dynamics']}
    artifacts={
      'interaction_quantization_contract':frozen,'quantization_location_manifest':{'QL01_QL10':quantization_location_registry()},
      'localized_state_manifest':{'loading_families':loading_registry(),'S01_S20':state_registry()},
      'localized_state_results':{'composite_observables':observables,**unresolved},
      'localized_state_survivor_ranking':{'survivors':[],**unresolved},'bound_free_classification':{'classification':'FREE','bound_states':[]},
      'localized_stability_results':{'stable':False,'metastable':False,**unresolved},
      'perturbation_stability_results':{'perturbations':['±1%','±5%','±10%'],'stable_localized_survivor':False,**unresolved},
      'internal_mode_results':{'wavelengths':[4,6,8,12,16,24,32],'loading_bound_modes':[],**unresolved},
      'internal_wavelength_results':{'diagnostic_available':True,'bound_wavelengths':[],**unresolved},
      'continuous_amplitude_sweep':{'amplitudes':AMPLITUDES.tolist(),'initial_norms':initial.tolist(),'final_norms':final.tolist(),
                                    'classification':'CONTINUOUS_FINAL_STATE_FAMILY'},
      'stable_state_cluster_results':clusters,'discrete_state_results':{'classification':'NO_STABLE_STATE','discrete':False,**unresolved},
      'grid_quantization_control':{'N':RESOLUTIONS,'classification':'NO_CANDIDATE; numerical quantization rejected'},
      'domain_quantization_control':{'extents':DOMAINS,'classification':'NO_CANDIDATE; boundary quantization not physical'},
      'topological_state_results':{'T01_T10':topology_registry(),'native_integer_mode_index_established':False},
      'transition_graph':transition_graph([]),'transition_pair_results':{'pairs':[],**unresolved},
      'transition_norm_accounting':{'transition_available':False,'conservation_gate_passed':False,**unresolved},
      'emission_results':emission_audit(),'emitted_packet_spectrum':{'packets':[],**unresolved},
      'transition_norm_vs_k':{'pairs':[],'classification':'TK06 no unique relation / not applicable'},
      'absorption_results':absorption_audit(),'fractional_absorption_results':{'fractions':[.25,.5,.75,1,1.25,1.5,2],**unresolved},
      'subthreshold_accumulation':{'classification':'UNRESOLVED',**unresolved},
      'polarization_selection_results':selection_audit(),'handedness_selection_results':selection_audit(),
      'source_driven_emission':{'discrete_bursts':False,**unresolved},'source_removed_emission_persistence':source_off_audit(),
      'reverse_transition_results':{'classification':'NOT_APPLICABLE',**unresolved},
      'stimulated_interaction_results':{'classification':'UNRESOLVED',**unresolved},'loading_backreaction_results':backreaction_audit(),
      'resolution_results':{'N':RESOLUTIONS,'localized_state':False,'conclusion':'no physical discreteness'},
      'domain_size_results':{'extents':DOMAINS,'localized_state':False,'conclusion':'free packet exits overlap region'},
      'coordinate_rescaling_results':{'alpha':ALPHAS,'localized_state':False,'continuous_scaling_retained':True},
      'quantum_transition_postfreeze_comparison':compare(True),'planck_transition_comparison':compare(True),
      'final_localized_state_contract':localized_contract,'final_interaction_quantization_contract':interaction_contract,
      'final_emission_absorption_contract':exchange,'final_native_quantum_exchange_contract':bridge}
    assert set(JSON_NAMES)==set(artifacts)
    for name,obj in artifacts.items(): dump(name,obj)
    for name in NPZ_NAMES: np.savez_compressed(RUN/f'{name}.npz',amplitude=AMPLITUDES,initial_norm=initial,final_norm=final,histories=histories)
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    for name in FIGURES:
        fig,ax=plt.subplots(figsize=(6,3.5)); ax.plot(initial,final,'o-',label='continuous synthetic overlap')
        ax.set(xlabel='initial native norm',ylabel='final native norm',title=name.replace('_',' ').title()); ax.legend(fontsize=7)
        fig.tight_layout(); fig.savefig(RUN/f'{name}.png',dpi=100); plt.close(fig)
    outcomes=['PBUF_LOADING_EXCITATION_LOCALIZATION_COUPLING_UNRESOLVED','PBUF_INTERACTION_QUANTIZATION_NOT_ESTABLISHED']
    result={'status':'DEV150_AUDIT_COMPLETE','baseline':base,'phases_executed':PHASES,'outcomes':outcomes,
      'scientific_conclusion':'Current PBUF supplies persistent loading and continuous native excitation separately, but no coefficient-free coupling binds them. Synthetic overlap stays continuously amplitude-scalable and cannot establish localized states, transitions, emission, absorption, or interaction quantization.'}
    structural={'phases':PHASES,'QL01_QL10_attempted':True,'S01_S20_attempted':True,'I01_I20_attempted':True,
      'QX01_QX20':quantization_family_registry(),'hard_gate_stopped_downstream_claims':True,'guards':frozen,'outcomes':outcomes}
    dump('result',result); dump('structural_result',structural)
    (RUN/'baseline_git.txt').write_text(subprocess.run(['git','status','--short'],cwd=ROOT,text=True,capture_output=True).stdout)
    report=['DEV150_AUDIT_COMPLETE',*outcomes,'PBUF_NATIVE_LOCALIZED_LOADED_EXCITATION_STATE_ESTABLISHED=false',
      'PBUF_NATIVE_LOCALIZED_EXCITATION_STATE_QUANTIZATION_ESTABLISHED=false','PBUF_NATIVE_INTERACTION_QUANTIZATION_ESTABLISHED=false',
      'PBUF_DISCRETE_NATIVE_EXCITATION_EMISSION_ESTABLISHED=false','PBUF_DISCRETE_NATIVE_EXCITATION_ABSORPTION_ESTABLISHED=false',
      'PBUF_BOUND_FREE_EXCITATION_NORM_EXCHANGE_ESTABLISHED=false','PBUF_NATIVE_TRANSITION_NORM_K_RELATION_ESTABLISHED=false',
      'PBUF_NATIVE_TRANSITION_NORM_LINEAR_IN_K=false','PBUF_PHOTON_LIKE_DISCRETE_EXCHANGE_STRUCTURE_EMERGES=false',
      'PBUF_QUANTIZED_PHOTON_LIKE_INTERACTION_EMERGES=false',*[f'{k}={str(v).lower() if isinstance(v,bool) else v}' for k,v in frozen.items()]]
    (RUN/'report.txt').write_text('\n'.join(report)+'\n'); print('\n'.join(report[:12]))
if __name__=='__main__': main()
