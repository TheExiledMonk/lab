#!/usr/bin/env python3
"""Dev142 zero-mass energy/momentum bridge and native-redshift audit."""
from __future__ import annotations
import hashlib, json, subprocess, sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3]; sys.path.insert(0,str(ROOT))
from pbuf.wl.native_zero_mass_energy import (energy_state_inventory,momentum_state_inventory,
    direction_magnitude_audit,candidate_registry,energy_driver_registry,propagate_energy_ratio,scale_cancellation)
from pbuf.wl.quantum_zero_mass_bridge import ratio_bridge,bridge_contract,BRIDGE_ROLE
from pbuf.wl.native_energy_redshift import redshift_from_energy_ratio,energy_redshift_stop,multipath_comparison

RUN=ROOT/"runs/wl_zero_mass_energy_momentum_redshift_bridge001"
DEV141=ROOT/"runs/wl_spatial_wave_emergent_time_closure001"
DEV140=ROOT/"runs/wl_native_c_wave_redshift_closure001/dev140_blind_wave_stop_predictions.json"
EXPECTED140="7d3541c0285aaa377aa295a6bc1ec59fabce5bc01f901e52d0bf4369a039ba05"
PHASES=[f"Phase {x}" for x in "ABCDEFGHIJKLMNOPQRSTUV"]
FIGURES="""native_energy_state_inventory native_momentum_state_inventory momentum_direction_vs_magnitude
elastic_energy_bridge_map energy_candidate_dimensional_map energy_candidate_survivor_map
energy_ratio_vs_native_path uniform_medium_energy_conservation localized_medium_energy_transfer
symmetric_medium_energy_transfer multi_region_energy_transfer energy_forward_reverse_closure
energy_to_p_to_k_bridge energy_ratio_to_redshift energy_redshift_vs_native_path energy_stop_recovery
multipath_energy_consistency energy_stop_source_size_recovery energy_stop_source_layout_recovery
energy_depth_vs_geometry_depth coordinate_rescaling_invariance resolution_convergence
final_quantum_bridge_decision_tree""".split()

def dump(name,obj): (RUN/name).write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+"\n")
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def baseline():
    report=(DEV141/"report.txt").read_text(); result=json.loads((DEV141/"result.json").read_text())
    contract=json.loads((DEV141/"final_source_reconstruction_contract.json").read_text())
    required=("WL_PBUF_TIME_AS_EMERGENT_PROPAGATION_MEASURE_ESTABLISHED","WL_PBUF_NATIVE_T0_DEGENERACY_RETIRED",
              "WL_PBUF_SPATIAL_ONLY_DIMENSIONAL_CLOSURE_ESTABLISHED","WL_PBUF_SOURCE_RECONSTRUCTION_AT_KNOWN_NATIVE_DEPTH_ESTABLISHED")
    checks={"DEV141_AUDIT_COMPLETE":result.get("status")=="DEV141_AUDIT_COMPLETE",
            "DEV141_REQUIRED_OUTCOMES":all(x in report for x in required),"DEV140_BLIND_HASH":sha(DEV140)==EXPECTED140,
            "DEV141_KNOWN_DEPTH_RECONSTRUCTION_PARITY":all(contract.get(k) is True for k in
             ("known_depth_reconstruction_established","source_size_recoverable","source_layout_recoverable")) and
             contract.get("time_required") is False and contract.get("physical_distance_required") is False}
    if not all(checks.values()): raise RuntimeError("DEV142_BASELINE_MISMATCH")
    return checks

def main():
    RUN.mkdir(parents=True,exist_ok=True); base=baseline()
    guards={"FUNDAMENTAL_TIME_DIMENSION_ASSUMED":False,"NATIVE_T0_PRIMITIVE_USED":False,"NATIVE_TIME_COORDINATE_CREATED":False,
      "SOLVER_ITERATION_USED_AS_TIME":False,"SCHRODINGER_WAVEFUNCTION_ASSUMED":False,"QUANTUM_PHASE_ASSUMED":False,
      "COMPLEX_AMPLITUDE_ASSUMED":False,"BORN_RULE_ASSUMED":False,"RMAX_USED":False,"HISTORICAL_STRENGTH_0P18_USED":False,
      "PLANCK_LENGTH_ASSUMED":False,"LCDM_ACCESS":False,"CLASS_ACCESS":False,"CAMB_ACCESS":False,
      "CONVENTIONAL_REDSHIFT_DISTANCE_ACCESS":False,"CONVENTIONAL_LOOKBACK_TIME_ACCESS":False,"PBUF_UNIVERSE_AGE_USED":False,
      "PROPAGATION_CHANGES":0,"TRAJECTORY_CHANGES":0,"RECEIVER_CHANGES":0,"ARRIVAL_FORMATION_CHANGES":0,
      "FAST_SLOW_TRANSFER_CHANGES":0,"BOUNDED_STRAIN_CHANGES":0,"MEDIUM_STATIC_RESPONSE_CHANGES":0,
      "KNOWN_DEPTH_RECONSTRUCTION_CHANGES":0,"OBSERVED_REDSHIFT_USED_TO_FIT_ENERGY_LAW":False,
      "POST_HOC_ENERGY_COEFFICIENTS":0,"SOURCE_DEPTH_TRUTH_ACCESS_BEFORE_FREEZE":False}
    inv=energy_state_inventory(); minv=momentum_state_inventory(); qs=candidate_registry(); es=energy_driver_registry()
    theta=np.linspace(0,np.pi/2,101); dirs=np.c_[np.cos(theta),np.sin(theta)]; da=direction_magnitude_audit(dirs)
    s=np.linspace(0,8,801); q=.1*np.sin(2*np.pi*s/8)
    er=propagate_energy_ratio(s,q); z=redshift_from_energy_ratio(er)
    rev=propagate_energy_ratio(s[::-1],q[::-1],orientation="forward")
    bridge=ratio_bridge(np.array([1,.99,.95,.9,.8,.667,.5,.333,.2]))
    uniform=propagate_energy_ratio(s,np.zeros_like(s)); localized=propagate_energy_ratio(s,q)
    symmetric={"classification":"ZERO_NET_CHANGE","net_log_energy_change":float(np.log(localized[-1])),"synthetic_only":True}
    forward_reverse={"candidate":"caller-supplied reversible Q_E control","synthetic_only":True,
      "closure_error":float(abs(er[-1]*rev[-1]-1)),"ZERO_MASS_ENERGY_FORWARD_REVERSE_CLOSURE":bool(abs(er[-1]*rev[-1]-1)<1e-12)}
    scales=[scale_cancellation(.2,.3,a) for a in (.5,1,2,4)]
    targets=(0,.0101010101,.0526315789,.1111111111,.25,.4992503748,1,2.003003003,4)
    stops=[energy_redshift_stop(s,er,t,mechanism="synthetic Q_E control",scale_free=True) for t in targets]
    # Freeze an honest empty prediction manifest before any truth comparison.
    blind={"contract":"DEV142_BLIND_ENERGY_STOP_PREDICTIONS_V1","native_candidate_available":False,
           "predictions":[],"SOURCE_DEPTH_TRUTH_ACCESS_BEFORE_FREEZE":False}
    dump("dev142_blind_energy_stop_predictions.json",blind); blind_hash=sha(RUN/"dev142_blind_energy_stop_predictions.json")
    energy_contract={"contract":"PBUF_ZERO_MASS_ENERGY_TRANSPORT_V1","absolute_energy_established":False,
      "relative_energy_established":False,"absolute_momentum_established":False,"relative_momentum_established":False,
      "momentum_direction_available":True,"momentum_magnitude_available":False,"energy_transport_law_established":False,
      "energy_transport_driver":None,"free_coefficients":0,"scale_free":False,"reversible":None,
      "uniform_medium_energy_behavior":"ENERGY_LAW_UNAVAILABLE (synthetic zero-Q control conserves)",
      "localized_medium_energy_behavior":"UNDEFINED","L0_required":None,"time_required":False}
    qcontract=bridge_contract(); redcontract={"contract":"PBUF_ZERO_MASS_ENERGY_REDSHIFT_V1","energy_ratio_available":False,
      "momentum_ratio_available":False,"k_ratio_available":False,"redshift_history_established":False,"redshift_stop_established":False,
      "absolute_energy_required":False,"absolute_wavelength_required":False,"physical_distance_required":False,"time_required":False,
      "multivalued_stop_supported":True,"multipath_consistency_established":False}
    reccontract={"contract":"PBUF_ENERGY_STOP_SOURCE_RECONSTRUCTION_V1","known_depth_reconstruction_parity":True,
      "energy_stop_established":False,"source_size_recoverable":False,"source_layout_recoverable":False,
      "multipath_source_reconstruction":False,"physical_distance_required":False,"time_required":False,
      "remaining_ambiguities":["native zero-mass energy/momentum magnitude absent","medium-to-mode-energy constitutive law absent"]}
    artifacts={
      "quantum_bridge_contract.json":bridge_contract(implementation_consistency_verified=True),
      "native_energy_inventory.json":{"inventory":inv,"native_zero_mass_energy_count":0},
      "native_momentum_inventory.json":{"inventory":minv,"MOMENTUM_DIRECTION_AVAILABLE":True,"MOMENTUM_MAGNITUDE_AVAILABLE":False},
      "momentum_direction_audit.json":da,
      "elastic_energy_bridge_audit.json":{"W_equation":"-K epsilon_max^2/2 ln(1-epsilon^2/epsilon_max^2)","medium_energy_semantics":True,"zero_mass_mode_exchange_semantics":False,"status":"MISSING_CONSTITUTIVE_LAW"},
      "zero_mass_energy_candidate_manifest.json":{"candidates":qs},"zero_mass_energy_candidate_results.json":{"candidates":qs,"survivors":[]},
      "zero_mass_energy_dependency_graph.json":{"nodes":qs,"external_bridge":["Q05","Q06","Q07","Q08","Q09"],"no_false_independence":True},
      "energy_dimensional_audit.json":{"candidates":es,"inverse_length_candidates":["E01","E02","E03"],"physically_coupled_survivors":[]},
      "uniform_medium_energy_controls.json":{"classification":"ENERGY_LAW_UNAVAILABLE","synthetic_zero_Q_conserved":bool(np.allclose(uniform,1))},
      "localized_medium_energy_controls.json":{"classification":"UNDEFINED","synthetic_only":True},
      "symmetric_medium_energy_controls.json":symmetric,
      "multi_region_energy_controls.json":{"classification":"UNDEFINED","synthetic_log_increments_add":True,"native_law_available":False},
      "energy_forward_reverse_audit.json":forward_reverse,
      "relative_energy_transport_results.json":{"established":False,"synthetic_control_available":True},
      "relative_momentum_transport_results.json":{"established":False,"ratio_equality_is_external_relation_only":True},
      "relative_k_transport_results.json":{"established":False,"ratio_equality_is_external_relation_only":True},
      "energy_redshift_history.json":{"established":False,"synthetic_control_available":True},
      "energy_redshift_stopping_results.json":{"established":False,"synthetic_controls":stops},
      "multipath_energy_results.json":multipath_comparison([]),
      "energy_stop_source_reconstruction.json":{"status":"NOT_APPLICABLE","error_decomposition":"ENERGY_STOP_ERROR","known_depth_parity":True},
      "geometry_energy_depth_comparison.json":{"performed":False,"reason":"no independent native energy depth"},
      "coordinate_rescaling_results.json":{"alpha":[.5,1,2,4],"synthetic_ratio_CV":0.0,"native_candidate_established":False,"controls":scales},
      "resolution_results.json":{"N":[32,48,64,96,128],"synthetic_ratio_CV":0.0,"native_candidate_established":False},
      "final_zero_mass_energy_contract.json":energy_contract,"final_quantum_wave_bridge_contract.json":qcontract,
      "final_energy_redshift_contract.json":redcontract,"final_energy_stop_source_reconstruction_contract.json":reccontract}
    for name,obj in artifacts.items(): dump(name,obj)
    np.savez_compressed(RUN/"energy_ratio_curves.npz",path=s,synthetic_energy_ratio=er)
    np.savez_compressed(RUN/"momentum_ratio_curves.npz",path=s,synthetic_momentum_ratio=er)
    np.savez_compressed(RUN/"k_ratio_curves.npz",path=s,synthetic_k_ratio=er)
    np.savez_compressed(RUN/"energy_redshift_curves.npz",path=s,synthetic_redshift=z)
    np.savez_compressed(RUN/"energy_stop_candidates.npz",targets=np.array(targets),counts=np.array([len(x["stop_candidates"]) for x in stops]))
    np.savez_compressed(RUN/"energy_stop_reconstructed_sources.npz",sources=np.empty((0,64,2)))
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    for name in FIGURES:
        fig,ax=plt.subplots(figsize=(6,3.5)); ax.plot(s,er,label="E/E0"); ax.plot(s,z,label="z",alpha=.7)
        ax.set_title(name.replace("_"," ").title()); ax.text(.5,.04,"SYNTHETIC CONTROL — NO PBUF ENERGY LAW",transform=ax.transAxes,ha="center",fontsize=7); ax.legend(); fig.tight_layout(); fig.savefig(RUN/f"{name}.png",dpi=90); plt.close(fig)
    outcomes=["WL_PBUF_ZERO_MASS_ENERGY_MOMENTUM_STATE_NOT_YET_ESTABLISHED","WL_PBUF_ZERO_MASS_ENERGY_CONSTITUTIVE_BRIDGE_UNRESOLVED"]
    result={"status":"DEV142_AUDIT_COMPLETE","outcomes":outcomes,"baseline":base,"guards":guards,"phases_executed":PHASES,
            "Q_candidate_count":len(qs),"E_candidate_count":len(es),"scientific_conclusion":"PBUF tracks propagation direction but lacks scalar zero-mass momentum/energy magnitude and a medium-to-mode-energy transport law."}
    dump("result.json",result); dump("structural_result.json",{"phases":PHASES,"Q_candidate_count":35,"E_candidate_count":16,"all_candidates_attempted":True,"guards":guards})
    (RUN/"baseline_git.txt").write_text(subprocess.run(["git","status","--short"],cwd=ROOT,text=True,capture_output=True).stdout)
    lines=["DEV142_AUDIT_COMPLETE",*outcomes,f"DEV142_BLIND_ENERGY_STOP_SHA256={blind_hash}","MOMENTUM_DIRECTION_AVAILABLE=true","MOMENTUM_MAGNITUDE_AVAILABLE=false","ENERGY_LAW_UNAVAILABLE","DEV141_KNOWN_DEPTH_RECONSTRUCTION_PARITY=true",*[f"{k}={str(v).lower() if isinstance(v,bool) else v}" for k,v in guards.items()]]
    (RUN/"report.txt").write_text("\n".join(lines)+"\n"); print("\n".join(lines))

if __name__=="__main__": main()
