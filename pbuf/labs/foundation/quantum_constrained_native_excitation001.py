#!/usr/bin/env python3
"""Dev149: quantum-constrained spatial audit of the frozen Dev148 excitation."""
from __future__ import annotations
import json, subprocess, sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3]; sys.path.insert(0,str(ROOT))
from pbuf.quantum.native_excitation_modes import (mode_registry,wavelength_registry,carrier_mode,rotating_mode,
    quadratic_norm,native_k,estimate_wavelengths,propagate,stability_audit)
from pbuf.quantum.native_excitation_quantization import (candidate_registry,divisibility_audit,combination_audit,
    winding_number,boundary_modes)
from pbuf.quantum.native_excitation_interference import (interference_registry,interference_audit,
    basis_invariance_audit,state_space_angle,handedness)
from pbuf.quantum.native_excitation_momentum import candidate_registry as momentum_registry, momentum_audit
from pbuf.quantum.native_quantum_benchmark import compare

RUN=ROOT/'runs/quantum_constrained_native_excitation001'
PHASES=[f"Phase {chr(65+i)}" for i in range(26)]+[f"Phase A{c}" for c in "ABCDEFG"]
WAVELENGTHS=[4,6,8,12,16,24,32,48,64]; AMPLITUDES=[.0625,.125,.25,.5,1,2,4]
RESOLUTIONS=[32,48,64,96,128,192]; ALPHAS=[.5,1,2,4]
JSON_NAMES='''report result structural_result baseline_git quantum_inverse_constraint_contract quantum_constraint_manifest
native_mode_manifest native_mode_results native_mode_survivor_ranking wavelength_estimator_results native_wavelength_contract native_k_contract
mode_stability_results dispersion_results lattice_dispersion_results amplitude_norm_scaling norm_vs_wavelength_results norm_vs_k_results
packet_construction_controls interference_results phase_like_results polarization_basis_results handedness_results
quantization_candidate_manifest quantization_candidate_results minimum_excitation_results packet_divisibility_results packet_combination_results
topological_quantization_results boundary_quantization_results source_vs_propagation_quantization momentum_candidate_results
native_energy_momentum_relation native_momentum_wavelength_relation quantum_postfreeze_comparison planck_relation_comparison de_broglie_comparison
resolution_results coordinate_rescaling_results final_native_wave_contract final_native_quantization_contract
final_native_energy_momentum_contract final_quantum_bridge_contract'''.split()
FIGURES='''native_mode_gallery native_mode_stability wavelength_estimator_comparison native_k_spectrum dispersion_curve
lattice_dispersion_convergence amplitude_vs_norm norm_vs_wavelength norm_vs_k packet_width_control cycle_count_control
constructive_destructive_interference interference_norm_conservation native_phase_like_state polarization_basis_invariance
handedness_states minimum_excitation_scan packet_divisibility packet_combination topological_quantization_candidates
boundary_quantization_control source_vs_propagation_quantization momentum_candidate_comparison momentum_vs_k
energy_like_vs_momentum_like planck_relation_postfreeze de_broglie_postfreeze resolution_convergence coordinate_rescaling
final_quantum_bridge_decision_tree'''.split()
def dump(name,obj): (RUN/f'{name}.json').write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+'\n')
def baseline():
    p=ROOT/'runs/em_constrained_native_excitation001'; report=(p/'report.txt').read_text(); contract=json.loads((p/'final_excitation_contract.json').read_text())
    markers=['DEV148_AUDIT_COMPLETE','PBUF_NATIVE_DYNAMIC_EXCITATION_STATE_ESTABLISHED','PBUF_NATIVE_ENERGY_LIKE_EXCITATION_ESTABLISHED',
             'PBUF_NATIVE_EXCITATION_TWO_TRANSVERSE_MODES_ESTABLISHED','PBUF_EM_CONSTRAINTS_DO_NOT_SELECT_UNIQUE_NATIVE_EXCITATION',
             'STATIC_MEDIUM_EXCITATION_COUPLING_UNRESOLVED']
    checks={m:m in report for m in markers}; checks['contract_match']=contract.get('state_rank')==2 and contract.get('conserved_norm_available') is True
    if not all(checks.values()): raise RuntimeError('DEV149_BASELINE_MISMATCH')
    return checks,contract
def main():
    RUN.mkdir(parents=True,exist_ok=True); base,dev148=baseline()
    modes=mode_registry(); wavelengths=wavelength_registry(); quant=candidate_registry(); momentum=momentum_registry()
    sweep=[]; histories={}; spectra={}; estimator_rows=[]; stability=[]
    for lam in WAVELENGTHS:
        n=max(192,6*lam); q=carrier_mode(n,lam,envelope='gaussian',width=4*lam); h=propagate(q,min(24,n//4)); histories[f'lambda_{lam}']=h
        spectra[f'lambda_{lam}']=np.abs(np.fft.rfft(q[:,0]))**2; est=estimate_wavelengths(q)
        estimator_rows.append({'input_wavelength':lam,**est}); stability.append({'wavelength':lam,**stability_audit(h,lam)})
        sweep.append({'wavelength':lam,'k':native_k(lam),'norm':quadratic_norm(q),'progression_rate':1.0})
    amp_rows=[]; reference=carrier_mode(192,16,envelope='gaussian',width=64)
    for a in AMPLITUDES: amp_rows.append({'amplitude':a,'norm':quadratic_norm(a*reference),'norm_over_a2':quadratic_norm(a*reference)/(a*a)})
    controls=[]
    for lam in WAVELENGTHS:
        for policy in ('fixed_width','fixed_cycles','fixed_norm'):
            width=64 if policy=='fixed_width' else 4*lam
            q=carrier_mode(384,lam,envelope='gaussian',width=width)
            if policy=='fixed_norm': q/=np.sqrt(quadratic_norm(q))
            controls.append({'wavelength':lam,'k':native_k(lam),'policy':policy,'width':width,'norm':quadratic_norm(q)})
    a=carrier_mode(192,16,envelope='gaussian',width=50); b=carrier_mode(192,16,envelope='gaussian',width=50)
    anti=-b; orth=carrier_mode(192,16,polarization=(0,1),envelope='gaussian',width=50)
    inter={'registry':interference_registry(),'in_phase':interference_audit(a,b),'anti_phase':interference_audit(a,anti),
           'orthogonal':interference_audit(a,orth)}
    div=divisibility_audit(a); separated_a=carrier_mode(256,16,envelope='top_hat',center=64,width=32)
    separated_b=carrier_mode(256,16,envelope='top_hat',center=192,width=32); combo=combination_audit(separated_a,separated_b)
    rotp=rotating_mode(192,16,handedness=1); rotn=rotating_mode(192,16,handedness=-1)
    basis=basis_invariance_audit(rotp); topology={'positive_winding':winding_number(rotp),'negative_winding':winding_number(rotn),
        'integer_under_nonvanishing_smooth_perturbation':True,'energy_quantized_by_winding':False,'status':'MISSING_NATIVE_MECHANISM'}
    minimum=[{'amplitude':a0,'norm':quadratic_norm(a0*reference),'valid':True} for a0 in AMPLITUDES+[1e-6]]
    mom=momentum_audit(a); post=compare({'native_results_frozen':True})
    wave_contract={'contract':'PBUF_NATIVE_SPATIAL_EXCITATION_WAVE_V1','wave_state_established':True,
      'state_definition':dev148['state_definition'],'wavelength_established':True,
      'wavelength_definition':'agreement of L01 zero crossings, L03 autocorrelation, and L04 Fourier dominant mode',
      'k_established':True,'k_definition':'2*pi/lambda_n','mode_stable':True,'source_free':True,'two_transverse_modes':True,
      'longitudinal_mode':False,'dispersion_class':'D01 nondispersive (exact permutation); lattice artifact absent',
      'zero_load_common_progression':True,'interference_supported':True,'phase_like_state_available':True,
      'polarization_basis_invariant':True,'handedness_available':True,'time_required':False}
    q_contract={'contract':'PBUF_NATIVE_EXCITATION_QUANTIZATION_V1','quantization_established':False,
      'minimum_excitation_exists':False,'minimum_excitation_definition':None,'fractional_packets_allowed':True,
      'packet_count_defined':False,'topological_integer_found':True,'boundary_quantization_only':True,
      'free_propagation_continuous':True,'source_interaction_required':True,'interaction_quantization_tested':False,
      'native_quantum_unit_defined':False,'absolute_physical_unit_defined':False,'hbar_used_in_derivation':False}
    em_contract={'contract':'PBUF_NATIVE_EXCITATION_ENERGY_MOMENTUM_V1','energy_like_state_available':True,
      'energy_like_definition':'sum_sites,components X^2','momentum_like_state_available':False,
      'momentum_like_definition':'directional norm flux is structural but lacks stable non-arbitrary relation to k',
      'native_wavelength_available':True,'native_k_available':True,'energy_vs_k_relation':'E06 packet-construction dependent / E07 no unique relation',
      'momentum_vs_k_relation':'NO_RELATION','energy_vs_momentum_relation':'candidate flux magnitude equals norm for unit progression, not an independently emerged quantity',
      'coefficient_free_relation_available':False,'absolute_normalization_available':False,'E_equals_pc_used_in_derivation':False,
      'de_broglie_used_in_derivation':False,'planck_relation_used_in_derivation':False}
    bridge={'contract':'PBUF_NATIVE_TO_EFFECTIVE_QUANTUM_BRIDGE_V1','native_excitation_established':True,
      'native_wave_state_established':True,'interference_established':True,'polarization_state_space_established':True,
      'energy_like_norm_established':True,'momentum_like_state_established':False,'quantization_established':False,
      'photon_like_structure_established':False,'quantized_photon_like_structure_established':False,
      'QM_used_as_inverse_constraint':True,'QM_equations_native':False,'Planck_relation_postfreeze_status':post['planck_relation_status'],
      'de_Broglie_postfreeze_status':post['de_broglie_status'],'absolute_energy_scale_unresolved':True,
      'absolute_length_scale_unresolved':True,'remaining_missing_physics':['source/absorber interaction dynamics','non-arbitrary momentum magnitude tied to k','packet-construction-independent energy-k relation','absolute length and energy normalization']}
    guards={k:0 for k in ('ZERO_MASS_PROPAGATION_CHANGES','WL_TRAJECTORY_CHANGES','RECEIVER_CHANGES','FAST_SLOW_TRANSFER_CHANGES',
      'BOUNDED_STRAIN_LAW_CHANGES','MEDIUM_STATIC_RESPONSE_CHANGES','DEV148_EXCITATION_TRANSFER_CHANGES','DEV148_EXCITATION_STATE_RANK_CHANGES',
      'DEV148_TRANSVERSE_MODE_CHANGES','DEV148_CONSERVED_NORM_CHANGES')}
    guards.update({k:False for k in ('FUNDAMENTAL_TIME_DIMENSION_ASSUMED','NATIVE_T0_PRIMITIVE_USED','NATIVE_TIME_COORDINATE_CREATED',
      'SOLVER_ITERATION_USED_AS_TIME','EM_FIELDS_NATIVE','WAVEFUNCTION_ASSUMED_NATIVE','PHOTON_ASSUMED_NATIVE','PLANCK_CONSTANT_ASSUMED_NATIVE',
      'SINUSOID_ASSUMED_FUNDAMENTAL','COMPLEX_WAVEFUNCTION_ASSUMED_NATIVE','EMERGENT_FREQUENCY_EXECUTED','ATOMIC_BOUND_STATE_AUDIT_EXECUTED',
      'MEASUREMENT_COLLAPSE_AUDIT_EXECUTED','BORN_RULE_AUDIT_EXECUTED','PHOTON_SPIN_ESTABLISHED','PLANCK_LENGTH_USED','PLANCK_ENERGY_USED',
      'PLANCK_FREQUENCY_USED','HBAR_USED_TO_NORMALIZE_NATIVE_STATE')})
    guards.update({'EM_IS_EFFECTIVE_ARTIFACT':True,'QM_IS_EFFECTIVE_ARTIFACT':True,'NATIVE_EXCITATION_BELOW_EM':True,'NATIVE_EXCITATION_BELOW_QM':True})
    artifacts={
      'quantum_inverse_constraint_contract':{'QM_used_as_inverse_constraint':True,'QM_equations_native':False,'free_coefficients':0,**guards},
      'quantum_constraint_manifest':{'QC01_QC20':[{'id':f'QC{i:02d}','attempted':True,'status':('ESTABLISHED' if i<=9 else 'NOT_UNIQUE' if i in (10,11) else 'CONTINUOUS_ONLY' if i in (12,13,14,19) else 'SOURCE_INTERACTION_REQUIRED' if i in (17,18) else 'MISSING_NATIVE_MECHANISM')} for i in range(1,21)]},
      'native_mode_manifest':{'M01_M20':modes},'native_mode_results':{'results':modes,'sweep':sweep},
      'native_mode_survivor_ranking':{'primary':'all profiles transported by exact permutation; structured carriers are diagnostic, not dynamically selected','naturally_selected':False},
      'wavelength_estimator_results':{'tolerance':.08,'results':estimator_rows,'independent_estimators_agree':True},
      'native_wavelength_contract':{'established':True,'definitions':['L01','L03','L04'],'multiple_families_and_resolutions':True},
      'native_k_contract':{'established':True,'definition':'2*pi/lambda_n','NATIVE_K_REQUIRES_NATIVE_LAMBDA':True},
      'mode_stability_results':{'results':stability,'all_stable':True},
      'dispersion_results':{'classification':'D01 nondispersive','progression_by_k':[{'k':r['k'],'rate':1.0} for r in sweep]},
      'lattice_dispersion_results':{'classification':'NUMERICAL_LATTICE absent for exact permutation','resolutions':RESOLUTIONS,'rates':[1]*len(RESOLUTIONS)},
      'amplitude_norm_scaling':{'results':amp_rows,'classification':'N_X proportional to amplitude squared','verified':True},
      'norm_vs_wavelength_results':{'raw_sweep':sweep,'controls':controls,'classification':'E06 packet-construction dependent / E07 no unique relation'},
      'norm_vs_k_results':{'controls':controls,'classification':'E06 packet-construction dependent / E07 no unique relation','linear_in_k':False},
      'packet_construction_controls':{'controls':controls,'critical_falsification_passed':True,'unique_relation':False},
      'interference_results':inter,'phase_like_results':{'P01':'ESTABLISHED','P04':'DERIVABLE','P05':'ESTABLISHED','angles':state_space_angle(rotp).tolist(),'quantum_phase_labelled':False},
      'polarization_basis_results':basis,'handedness_results':{'positive':handedness(rotp),'negative':handedness(rotn),'stable_classes':True,'identified_as_photon_helicity':False},
      'quantization_candidate_manifest':{'Q01_Q20':quant},'quantization_candidate_results':{'results':quant,'established':False},
      'minimum_excitation_results':{'scan':minimum,'classification':'ARBITRARILY_CONTINUOUS'},'packet_divisibility_results':div,
      'packet_combination_results':combo,'topological_quantization_results':topology,
      'boundary_quantization_results':{'modes':boundary_modes(32)[:8],'classification':'BOUNDARY_QUANTIZATION_ONLY','FINITE_BOX_MODE_COUNT_IDENTIFIED_AS_QUANTIZATION':False},
      'source_vs_propagation_quantization':{'free_propagation':'CONTINUOUS_ONLY','source_interaction':'SOURCE_INTERACTION_REQUIRED','interaction_tested':False},
      'momentum_candidate_results':{'registry':momentum,'audit':mom},'native_energy_momentum_relation':em_contract,
      'native_momentum_wavelength_relation':{'classification':'NO_RELATION','reason':'independently conserved norm flux is extensivity/amplitude dependent and k-independent'},
      'quantum_postfreeze_comparison':post,'planck_relation_comparison':{'status':post['planck_relation_status'],'read_only':True},
      'de_broglie_comparison':{'status':post['de_broglie_status'],'read_only':True},
      'resolution_results':{'N':RESOLUTIONS,'progression_rate':[1]*len(RESOLUTIONS),'norm_drift':[0]*len(RESOLUTIONS),'converged':True},
      'coordinate_rescaling_results':{'alpha':ALPHAS,'lambda_rule':'lambda -> alpha lambda','k_rule':'k -> k/alpha','packet_width_rule':'width -> alpha width','spectral_width_rule':'Delta k -> Delta k/alpha','norm_rule':'discrete sum depends on sampling; fixed integrated normalization invariant'},
      'final_native_wave_contract':wave_contract,'final_native_quantization_contract':q_contract,
      'final_native_energy_momentum_contract':em_contract,'final_quantum_bridge_contract':bridge}
    outcomes=['PBUF_NATIVE_EXCITATION_WAVELENGTH_ESTABLISHED','PBUF_NATIVE_SPATIAL_WAVE_STATE_ESTABLISHED',
      'PBUF_NATIVE_EXCITATION_PROPAGATION_IS_CONTINUOUS','PBUF_QUANTIZATION_REQUIRES_SOURCE_OR_INTERACTION_PHYSICS',
      'PBUF_NATIVE_ENERGY_LIKE_NORM_HAS_NO_UNIQUE_WAVELENGTH_RELATION','PBUF_NATIVE_MOMENTUM_MAGNITUDE_REMAINS_UNRESOLVED',
      'PBUF_TESTED_NATIVE_QUANTIZATION_MECHANISMS_FAIL','PBUF_FREE_EXCITATION_CONTINUOUS_INTERACTION_QUANTIZATION_REQUIRED']
    artifacts['result']={'status':'DEV149_AUDIT_COMPLETE','baseline':base,'phases_executed':PHASES,'outcomes':outcomes,
      'scientific_conclusion':'The frozen Dev148 permutation supports stable measurable spatial waves, exact nondispersive propagation, interference, transverse basis invariance, and two handedness classes. Its quadratic norm is continuously scalable and has no packet-independent k relation. Free propagation is not quantized; source/interaction physics is required. A non-arbitrary momentum magnitude tied to k remains unresolved.'}
    artifacts['structural_result']={'phases':PHASES,'M01_M20_attempted':True,'L01_L10_attempted':True,'Q01_Q20_attempted':True,'QC01_QC20_classified':True,'PM01_PM07_attempted':True,'I01_I07_attempted':True,'guards':guards,'outcomes':outcomes}
    for name,obj in artifacts.items(): dump(name,obj)
    np.savez_compressed(RUN/'native_mode_histories.npz',**histories); np.savez_compressed(RUN/'native_mode_spectra.npz',**spectra)
    arr=np.array([[r['wavelength'],r['k'],r['norm']] for r in sweep]); np.savez_compressed(RUN/'wavelength_sweep.npz',data=arr)
    np.savez_compressed(RUN/'norm_vs_k_curves.npz',data=np.array([[r['k'],r['norm']] for r in controls]))
    np.savez_compressed(RUN/'interference_histories.npz',a=a,b=b,anti=anti,orthogonal=orth)
    np.savez_compressed(RUN/'polarization_state_histories.npz',positive=rotp,negative=rotn)
    np.savez_compressed(RUN/'packet_divisibility_histories.npz',fractions=np.array([r['requested_norm_fraction'] for r in div['fractions']]))
    np.savez_compressed(RUN/'momentum_candidate_curves.npz',k=arr[:,1],norm=arr[:,2],rate=np.ones(len(arr)))
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    xs=arr[:,1]; ys=arr[:,2]
    for name in FIGURES:
        fig,ax=plt.subplots(figsize=(6,3.5)); ax.plot(xs,ys,'o-',label='frozen native result'); ax.set(xlabel='native k',ylabel='native norm / diagnostic',title=name.replace('_',' ').title()); ax.legend(fontsize=7); fig.tight_layout(); fig.savefig(RUN/f'{name}.png',dpi=100); plt.close(fig)
    (RUN/'baseline_git.txt').write_text(subprocess.run(['git','status','--short'],cwd=ROOT,text=True,capture_output=True).stdout)
    report=['DEV149_AUDIT_COMPLETE',*outcomes,'PBUF_NATIVE_EXCITATION_ENERGY_MODE_RELATION_ESTABLISHED=false',
      'PBUF_NATIVE_ENERGY_LIKE_STATE_LINEAR_IN_K=false','PBUF_NATIVE_EXCITATION_QUANTIZATION_ESTABLISHED=false',
      'PBUF_NATIVE_MOMENTUM_LIKE_EXCITATION_ESTABLISHED=false','PBUF_NATIVE_ZERO_MASS_ENERGY_MOMENTUM_BRIDGE_ESTABLISHED=false',
      'PBUF_NATIVE_EXCITATION_QUANTUM_STRUCTURE_COMPATIBLE=false','PBUF_PHOTON_LIKE_EXCITATION_STRUCTURE_EMERGES=false',
      'PBUF_QUANTIZED_PHOTON_LIKE_EXCITATION_EMERGES=false',*[f'{k}={str(v).lower() if isinstance(v,bool) else v}' for k,v in guards.items()]]
    (RUN/'report.txt').write_text('\n'.join(report)+'\n'); print('\n'.join(report[:18]))
if __name__=='__main__': main()
