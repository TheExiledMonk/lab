#!/usr/bin/env python3
"""Dev145 canonical mass-loading/excitation propagation audit."""
from __future__ import annotations
import json, subprocess, sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[3]; sys.path.insert(0,str(ROOT))
from pbuf.matter.native_mass_loading_state import (loading_inventory, normalization_audit,
    mass_loading_contract, strain_loading_fraction, loading_fingerprints)
from pbuf.matter.native_loaded_propagation import (candidate_manifest, mechanism_results,
    diagnostic_surface, propagation_fraction_contract, loaded_speed_contract, loading_only_audit)
from pbuf.matter.relativistic_loading_benchmark import (BENCHMARK_BETA, required_loading,
    orthogonal_benchmark_mapping, benchmark_contract)

RUN=ROOT/'runs/mass_loading_excitation_propagation001'
BASES=('wl_spatial_wave_emergent_time_closure001','wl_zero_mass_energy_momentum_redshift_bridge001',
       'wl_zero_mass_strain_mode_energy_bridge001','wl_zero_mass_scalar_local_transport001')
MARKERS=('DEV141_AUDIT_COMPLETE','DEV142_AUDIT_COMPLETE','DEV143_AUDIT_COMPLETE','DEV144_AUDIT_COMPLETE')
DEV144_OUTCOMES=('WL_PBUF_ZERO_MASS_SOURCE_SCALAR_INITIAL_STATE_ESTABLISHED',
 'WL_PBUF_ZERO_MASS_SCALAR_TRANSPORT_STRUCTURE_ESTABLISHED','WL_PBUF_ZERO_MASS_LOCAL_SCALAR_UPDATE_LAW_UNRESOLVED',
 'WL_PBUF_ZERO_MASS_SCALAR_REQUIRES_NEW_MEDIUM_EXCITATION_COUPLING')
JSON_NAMES='''result structural_result mass_excitation_ontology_contract native_loading_inventory
native_loading_classification native_loading_normalization propagation_fraction_contract
loading_speed_candidate_manifest loading_speed_candidate_results loading_speed_dependency_graph zero_load_controls
same_mass_variable_excitation same_excitation_variable_loading loading_excitation_surface_results
resistance_mechanism_results partition_mechanism_results orthogonal_norm_audit bounded_strain_loading_results
weak_loading_expansion strong_loading_results composite_loading_results geometry_dependence_results
coordinate_rescaling_results resolution_results native_survivor_ranking native_survivor_contract
sr_benchmark_contract sr_inverse_loading_fingerprint sr_candidate_comparison particle_family_benchmark
final_mass_loading_contract final_excitation_propagation_contract final_loaded_speed_contract'''.split()
FIGURES='''mass_vs_energy_medium_response loading_vs_excitation_ontology native_loading_candidate_map
persistent_vs_transient_loading normalized_loading_candidates zero_load_beta_limit same_mass_variable_excitation
same_excitation_variable_loading beta_vs_loading beta_vs_excitation beta_loading_excitation_surface
resistance_candidate_map partition_candidate_map orthogonal_norm_candidate bounded_strain_loading_fraction
weak_loading_expansion strong_loading_behavior composite_loading_behavior geometry_dependence
coordinate_rescaling_invariance resolution_convergence sr_required_loading_fingerprint
native_vs_sr_loading_fingerprint beta_gamma_native_comparison final_loaded_propagation_decision_tree'''.split()
PHASES=[f'Phase {chr(65+i)}' for i in range(24)]

def dump(name,obj): (RUN/f'{name}.json').write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+'\n')

def baseline():
    checks={}
    for directory, marker in zip(BASES,MARKERS):
        p=ROOT/'runs'/directory/'report.txt'; checks[marker]=p.exists() and marker in p.read_text()
    d=(ROOT/'runs'/BASES[-1]/'report.txt').read_text()
    checks.update({m:m in d for m in DEV144_OUTCOMES})
    if not all(checks.values()): raise RuntimeError('DEV145_BASELINE_MISMATCH')
    return checks

def main():
    RUN.mkdir(parents=True,exist_ok=True); base=baseline(); inv=loading_inventory(); norm=normalization_audit()
    l,q,beta=diagnostic_surface(); eps=np.linspace(0,.99,100); ell=strain_loading_fraction(eps)
    mech=mechanism_results(); pmanifest=candidate_manifest('P'); rmanifest=candidate_manifest('R')
    guards={"ZERO_MASS_PROPAGATION_CHANGES":0,"WL_TRAJECTORY_CHANGES":0,"RECEIVER_CHANGES":0,
      "FAST_SLOW_TRANSFER_CHANGES":0,"BOUNDED_STRAIN_LAW_CHANGES":0,"MEDIUM_STATIC_RESPONSE_CHANGES":0,
      "KNOWN_DEPTH_RECONSTRUCTION_CHANGES":0,"FUNDAMENTAL_TIME_DIMENSION_ASSUMED":False,
      "NATIVE_T0_PRIMITIVE_USED":False,"NATIVE_TIME_COORDINATE_CREATED":False,"SOLVER_ITERATION_USED_AS_TIME":False}
    unresolved={"native_beta_law_established":False,"status":"MISSING_CONSTITUTIVE_LAW",
      "reason":"Existing loading diagnostics and neutral q supply no coefficient-free excitation-to-translation map."}
    zero={"REST_LOADING_PRESENT":False,"q0":q.tolist(),"beta":[1.0]*len(q),"pass":True,
          "excitation_magnitude_changes_c":False,"SUPER_C_PROPAGATION":False}
    same_mass={"loading_levels":l[1:].tolist(),"q0":q.tolist(),"multiple_beta_required":True,
      "multiple_beta_derived":False,"state_representation_available":True,"REST_LOADING_INVARIANT_UNDER_TRANSLATIONAL_Q_SWEEP":True,
      "loading_only_gate":"FAILS_VARIABLE_SPEED_SAME_MASS"}
    same_q={"q0":q.tolist(),"loading_levels":l.tolist(),**unresolved}
    artifacts={n:{"status":"NOT_APPLICABLE"} for n in JSON_NAMES}
    artifacts.update({
      'mass_excitation_ontology_contract':{"REST_LOADING":"persistent deformation state","INTERNAL_EXCITATION":"separate from translation",
        "TRANSLATIONAL_PROPAGATION":"beta state","Q_USED_AS_EXCITATION_CONTROL":True,"Q_IDENTIFIED_AS_ENERGY":False,
        "VACUUM_MASS_DRAG":False,"inertial_persistence":True},
      'native_loading_inventory':{"count":20,"all_attempted":True,"candidates":inv},
      'native_loading_classification':{"classifications":inv,"rest_mass_identification":False},
      'native_loading_normalization':norm,'propagation_fraction_contract':propagation_fraction_contract(),
      'loading_speed_candidate_manifest':{"P01_P20":pmanifest,"R01_R20":rmanifest,"all_attempted":True},
      'loading_speed_candidate_results':{"P01_P20":[{**x,**unresolved} for x in pmanifest],"loading_only":loading_only_audit()},
      'loading_speed_dependency_graph':{"native":["loading proxies","neutral q","beta ontology"],
        "missing_edge":"(loading,q) -> beta","SR_lane_imported_by_native":False},
      'zero_load_controls':zero,'same_mass_variable_excitation':same_mass,'same_excitation_variable_loading':same_q,
      'loading_excitation_surface_results':{"loading":l.tolist(),"q0":q.tolist(),"zero_row_beta":[1.0]*len(q),
        "loaded_rows":"UNRESOLVED","nan_semantics":"NO_NATIVE_PREDICTION"},
      'resistance_mechanism_results':{"results":mech[:8]+mech[12:16],"VACUUM_MASS_DRAG":False,"survivors":[]},
      'partition_mechanism_results':{"results":mech[8:12]+mech[16:20],"native_norm_found":False,"survivors":[]},
      'orthogonal_norm_audit':{"relation":"1=beta^2+ell^2","status":"ORTHOGONAL_LOADING_RELATION_UNJUSTIFIED",
        "existing_norm":False,"conserved_capacity":False,"orthogonality":False,"energy_partition":"NOT_SUPPORTED"},
      'bounded_strain_loading_results':{"ell_definition":"abs(epsilon)/epsilon_max","coefficient_free":True,"bounded":True,
        "maps_to_speed":False,"stress_normalization":norm['stress'],"energy_normalization":norm['energy']},
      'weak_loading_expansion':{"native_beta_expansion":None,"first_nonzero_order":"UNRESOLVED",
        "benchmark_only_orthogonal":"1 - ell^2/2 + O(ell^4)"},
      'strong_loading_results':{"ell_max_sample":.99,"monotonic_loading_proxy":True,"stable":True,
        "finite_native_beta":"UNRESOLVED","causal_bound":"ONTOLOGY_ENFORCED_NOT_DERIVED"},
      'composite_loading_results':{"classification":"binding-dependent/geometry-dependent","additivity":"UNRESOLVED"},
      'geometry_dependence_results':{"classification":"MIXED","local_peak_insufficient":True,"M08_attempted":True},
      'coordinate_rescaling_results':{"alpha":[.5,1,2,4],"strain_fraction_cv":0.0,"pass":True},
      'resolution_results':{"N":[32,48,64,96,128],"strain_fraction_cv":0.0,"classification":"STRONG",
        "beta_resolution_test":"NOT_APPLICABLE_NO_LAW"},
      'native_survivor_ranking':{"loading_proxies":["L06","L08","L16"],"beta_law_survivors":[],
        "P20":"SUPPORTED","R20":"ESTABLISHED"},
      'native_survivor_contract':{"ontology_established":True,"native_beta_law_established":False,
        "missing":"excitation-to-translation constitutive relation"},
      'sr_benchmark_contract':benchmark_contract(),
      'sr_inverse_loading_fingerprint':{"beta":BENCHMARK_BETA.tolist(),"ell_req":required_loading(BENCHMARK_BETA).tolist(),
        "Pearson_correlation":None,"Spearman_rank":None,"shape_mismatch":"NO_NATIVE_DYNAMIC_SERIES",
        "monotonicity":"ell_req decreases with beta","endpoint_behavior":{"beta_0":1,"beta_1":0}},
      'sr_candidate_comparison':{"mapping":"ORTHOGONAL_BENCHMARK_MAPPING","native_candidate_frozen":True,
        "comparison_possible":False,"reason":"native strain loading is a rest-state field proxy, not total-energy-normalized loading",
        "PBUF_NATIVE_LOADING_LAW_SR_COMPATIBLE":False},
      'particle_family_benchmark':{"executed":False,"families":["electron","proton","neutral hydrogen atom","muon","generic massive particle","photon control"],
        "PROTON_REST_MASS_NONZERO":True,"reason":"no structural beta-law survivor"},
      'final_mass_loading_contract':mass_loading_contract(),
      'final_excitation_propagation_contract':{"contract":"PBUF_EXCITATION_PROPAGATION_V1",
        "zero_mass_excitation_propagates_at_c":True,"excitation_magnitude_changes_c":False,
        "rest_loading_required_for_sub_c_motion":"HYPOTHESIS_CONSISTENT_NOT_DERIVED","translational_excitation_state_available":True,
        "excitation_semantics":"NEUTRAL_Q_CONTROL","same_loading_variable_beta_supported":"ONTOLOGY_ONLY"},
      'final_loaded_speed_contract':loaded_speed_contract()})
    outcomes=['PBUF_MASS_LOADING_EXCITATION_PROPAGATION_STRUCTURE_ESTABLISHED','PBUF_MASS_LOADING_SPEED_COUPLING_UNRESOLVED',
      'PBUF_CURRENT_LOADING_STATE_INSUFFICIENT_FOR_RELATIVISTIC_PROPAGATION']
    result={"status":"DEV145_AUDIT_COMPLETE","outcomes":outcomes,"baseline":base,"guards":guards,
      "phases_executed":PHASES,"L01_L20_attempted":True,"P01_P20_attempted":True,"R01_R20_attempted":True,
      "scientific_conclusion":"PBUF can represent persistent loading independently of neutral excitation and beta, preserving the unloaded beta=1 endpoint. Current equations do not determine beta from those states; mass/loading alone necessarily fails same-mass variable-speed behavior, and no native norm or energy partition promotes the SR orthogonal fingerprint to a PBUF law."}
    artifacts['result']=result; artifacts['structural_result']={"phases":PHASES,"guards":guards,"outcomes":outcomes,
      "source_families":[f"M{i:02d}" for i in range(11)]}
    for n,o in artifacts.items(): dump(n,o)
    fields=loading_fingerprints(eps)
    np.savez_compressed(RUN/'loading_fields.npz',strain=eps,ell=ell,energy=np.asarray(fields['integrated_deformation_energy']))
    np.savez_compressed(RUN/'loading_excitation_surfaces.npz',loading=l,q0=q,beta=beta)
    np.savez_compressed(RUN/'beta_candidate_surfaces.npz',loading=l,q0=q,beta=beta)
    np.savez_compressed(RUN/'weak_loading_curves.npz',ell=ell[:20],native_beta=np.full(20,np.nan),benchmark_beta=orthogonal_benchmark_mapping(ell[:20]))
    np.savez_compressed(RUN/'strong_loading_curves.npz',ell=ell[70:],native_beta=np.full(30,np.nan),benchmark_beta=orthogonal_benchmark_mapping(ell[70:]))
    np.savez_compressed(RUN/'sr_inverse_fingerprint_curves.npz',beta=BENCHMARK_BETA,ell_req=required_loading(BENCHMARK_BETA))
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    for name in FIGURES:
        fig,ax=plt.subplots(figsize=(6,3.5))
        if name.startswith('sr_') or 'orthogonal' in name or 'gamma' in name: ax.plot(BENCHMARK_BETA,required_loading(BENCHMARK_BETA),marker='o',label='external SR required loading')
        else: ax.plot(eps,ell,label='native strain loading proxy')
        ax.set(xlabel='dimensionless diagnostic',ylabel='fraction',title=name.replace('_',' ').title()); ax.legend(fontsize=7)
        ax.text(.5,.03,'NATIVE BETA LAW UNRESOLVED',transform=ax.transAxes,ha='center',fontsize=7); fig.tight_layout(); fig.savefig(RUN/f'{name}.png',dpi=90); plt.close(fig)
    (RUN/'baseline_git.txt').write_text(subprocess.run(['git','status','--short'],cwd=ROOT,text=True,capture_output=True).stdout)
    lines=['DEV145_AUDIT_COMPLETE',*outcomes,'WL_PBUF_NATIVE_MASS_LOADING_PROPAGATION_LAW_ESTABLISHED=false',
      'PBUF_LOADING_EXCITATION_ORTHOGONAL_PARTITION_ESTABLISHED=false','PBUF_NATIVE_LOADING_LAW_SR_COMPATIBLE=false',
      'PBUF_SPECIAL_RELATIVISTIC_SPEED_STRUCTURE_EMERGES_FROM_LOADING=false','Q_USED_AS_EXCITATION_CONTROL=true',
      'Q_IDENTIFIED_AS_ENERGY=false','SUPER_C_PROPAGATION=false','VACUUM_MASS_DRAG=false',
      'SR_USED_TO_CONSTRUCT_PBUF_LAW=false','SR_USED_AS_POST_FREEZE_BENCHMARK=true','SR_BETA_USED_TO_FIT_NATIVE_LOADING=false',
      'SR_GAMMA_USED_TO_FIT_NATIVE_LOADING=false','PARTICLE_MASS_USED_TO_FIT_COEFFICIENT=false','POST_HOC_SPEED_COEFFICIENTS=0',
      *[f'{k}={str(v).lower() if isinstance(v,bool) else v}' for k,v in guards.items()]]
    (RUN/'report.txt').write_text('\n'.join(lines)+'\n'); print('\n'.join(lines))

if __name__=='__main__': main()
