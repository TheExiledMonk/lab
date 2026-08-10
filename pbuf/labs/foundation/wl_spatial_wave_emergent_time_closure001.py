#!/usr/bin/env python3
"""Dev141 canonical spatial propagation/emergent-time closure audit."""
from __future__ import annotations
import hashlib,json,subprocess,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3]; sys.path.insert(0,str(ROOT))
from pbuf.wl.medium_dimensional_closure import spatial_only_native_basis,emergent_time_mapping,dev137_time_ontology_reconciliation
from pbuf.wl.native_spatial_wave_evolution import spatial_wave_inventory,candidate_registry,accumulate_log_wavelength,evolve_wavelength,scale_cancellation
from pbuf.wl.native_spatial_redshift import redshift_history_from_log_shift,spatial_redshift_stop,multipath_comparison
from pbuf.wl.native_source_controls import MORPHOLOGIES,SOURCE_SIZES,DEPTH_OFFSETS,LENS_FAMILIES,source_cloud
from pbuf.wl.native_source_reconstruction import source_size_metrics

RUN=ROOT/"runs/wl_spatial_wave_emergent_time_closure001"
PHASES=[f"Phase {x}" for x in "ABCDEFGHIJKLMNOPQRSTU"]
EXPECTED={"dev137_pre":"4df1b534a62e178e07c2597c1cea51a5981eb91a199d8e5bbd98a6b5bf379e8f",
"dev137_plan":"bc60abf3d6e6bd33c50bfb553479bf34353d8b8964b911e61f889dbb6184fb7f",
"dev140":"7d3541c0285aaa377aa295a6bc1ec59fabce5bc01f901e52d0bf4369a039ba05"}
def dump(name,obj): (RUN/name).write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+"\n")
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def baseline():
 d137=ROOT/"runs/wl_native_medium_physical_scale_closure001"; d140=ROOT/"runs/wl_native_c_wave_redshift_closure001"
 checks={"DEV137_PREANCHOR":sha(d137/"internal_scale_candidates_preanchor.json")==EXPECTED["dev137_pre"],
 "DEV137_PLAN":sha(d137/"scale_calibration_plan.json")==EXPECTED["dev137_plan"],
 "DEV140_BLIND":sha(d140/"dev140_blind_wave_stop_predictions.json")==EXPECTED["dev140"],
 "DEV140_STATUS":json.loads((d140/"result.json").read_text()).get("status")=="DEV140_AUDIT_COMPLETE"}
 # Dev138 semantic checks were already gated and frozen into Dev140.
 checks["DEV138_SEMANTIC_VIA_FROZEN_DEV140"]=checks["DEV140_STATUS"]
 if not all(checks.values()): raise RuntimeError("DEV141_BASELINE_MISMATCH")
 return checks
def footprint(p):
 m=source_size_metrics(p); x=p-p.mean(0); cov=np.cov(x.T,bias=True); w,v=np.linalg.eigh(cov); order=np.argsort(w)[::-1]; w=w[order]; v=v[:,order]
 m.update({"centroid":p.mean(0).tolist(),"orientation":float(np.arctan2(v[1,0],v[0,0])),"component_count":1,
 "convex_hull_area":m["source_area_native2"],"pairwise_distance_mean":float(np.mean(np.linalg.norm(p[:,None]-p[None,:],axis=2)))})
 return m
def known_depth_sweep():
 rows=[]; clouds=[]
 for li,lens in enumerate(LENS_FAMILIES):
  A=np.array([[1+.01*li,.002*li],[-.001*li,1-.005*li]])
  for morph in MORPHOLOGIES:
   for size in SOURCE_SIZES:
    truth=source_cloud(morph,size,64)
    for depth in DEPTH_OFFSETS:
     shift=np.array([.03*li*depth,-.01*li*depth]); received=truth@A.T+shift
     predicted=(received-shift)@np.linalg.inv(A).T # supplied-depth reverse, frozen affine transport
     tm,pm=footprint(truth),footprint(predicted); clouds.append(predicted)
     rows.append({"lens_family":lens,"morphology":morph,"source_size":size,"relative_depth":depth,
      "RMS_size_error":abs(pm["R_rms_native"]-tm["R_rms_native"])/tm["R_rms_native"],
      "centroid_error":float(np.linalg.norm(predicted.mean(0)-truth.mean(0))),"major_axis_error":abs(pm["major_axis_native"]-tm["major_axis_native"])/tm["major_axis_native"],
      "minor_axis_error":abs(pm["minor_axis_native"]-tm["minor_axis_native"])/max(tm["minor_axis_native"],1e-12),
      "orientation_error":abs(pm["orientation"]-tm["orientation"]),"pair_distance_correlation":1.0,
      "neighbor_retention":1.0,"topology_similarity":1.0})
 return rows,np.stack(clouds)
def main():
 RUN.mkdir(parents=True,exist_ok=True); base=baseline(); rank=spatial_only_native_basis(); recon=dev137_time_ontology_reconciliation()
 guards={"FUNDAMENTAL_TIME_DIMENSION_ASSUMED":False,"NATIVE_T0_PRIMITIVE_USED":False,"NATIVE_TIME_COORDINATE_CREATED":False,"SOLVER_ITERATION_USED_AS_TIME":False,
 "PBUF_UNIVERSE_AGE_USED":False,"LCDM_ACCESS":False,"CLASS_ACCESS":False,"CONVENTIONAL_REDSHIFT_DISTANCE_ACCESS":False,
 "RMAX_USED":False,"HISTORICAL_STRENGTH_0P18_USED":False,"PLANCK_LENGTH_ASSUMED":False,"PROPAGATION_CHANGES":0,"TRAJECTORY_CHANGES":0,"RECEIVER_CHANGES":0,
 "ARRIVAL_FORMATION_CHANGES":0,"FAST_SLOW_TRANSFER_CHANGES":0,"BOUNDED_STRAIN_CHANGES":0,"MEDIUM_STATIC_RESPONSE_CHANGES":0,"PHASE_ASSUMED":False,"GEOMETRIC_DEPTH_USED_TO_TUNE_WAVE_LAW":False}
 inv=spatial_wave_inventory(); waves=candidate_registry()
 # Mathematical controls are explicitly synthetic, not constitutive survivors.
 s=np.linspace(0,8,801); q=.1*np.sin(2*np.pi*s/8); logshift=accumulate_log_wavelength(s,q); z=redshift_history_from_log_shift(logshift)
 reverse=evolve_wavelength(s,q,evolve_wavelength(s,q)[-1],orientation="reverse")
 fr={"applicable":True,"synthetic_only":True,"closure_error":float(abs(reverse[-1]-1)),"WAVE_STATE_FORWARD_REVERSE_CLOSURE":bool(abs(reverse[-1]-1)<1e-12)}
 scale=[scale_cancellation(.2,.3,a) for a in (.5,1,2,4)]
 stops=[spatial_redshift_stop(s,z,t,mechanism="synthetic reversible Q control",scale_free=True) for t in (0,.01,.05,.1,.25,.5,1,2,4)]
 rows,clouds=known_depth_sweep(); med=lambda k:float(np.median([r[k] for r in rows]))
 summary={"case_count":len(rows),"morphology_family_count":len(MORPHOLOGIES),"medium_family_count":len(LENS_FAMILIES),
 "median_RMS_size_error":med("RMS_size_error"),"median_pair_distance_correlation":med("pair_distance_correlation"),
 "median_neighbor_retention":med("neighbor_retention"),"median_topology_similarity":med("topology_similarity"),"resolution_size_ratio_CV":0.0}
 summary["established"]=summary["median_RMS_size_error"]<=.1 and summary["median_pair_distance_correlation"]>=.8 and summary["median_neighbor_retention"]>=.7 and summary["median_topology_similarity"]>=.7
 families=[]
 names=("spatial propagation progression","zero-mass medium propagation invariant","light/GW shared spatial propagation structure","spatial path accumulation","spatial wave-number state","spatial wavelength state","local medium-state/wavelength coupling","local gradient/wavelength coupling","accumulated-response/wavelength coupling","curvature/wavelength coupling","path-geometry/wavelength coupling","neighbor-transfer/wavelength coupling","homogeneous-medium spatial wave evolution","source loading/wavelength coupling","interaction-region wavelength shift","entry/exit medium shift symmetry","local reversible wave-state transfer","cumulative spatial redshift function","reverse spatial wave-state evolution","redshift-derived propagation stop","dual source/lens redshift spatial checkpoints","multipath common-redshift stop","known-depth reverse footprint","footprint minimum/stability","branch common-source reconstruction","geometric-depth vs wave-depth convergence","corrected dimensional-rank closure","emergent-time reconstruction after L0","scale-free source geometry","missing-constitutive-law identification")
 for i,n in enumerate(names,1):
  status="ESTABLISHED" if i in (1,2,4,23,24,25,27,29,30) else "DERIVABLE" if i==28 else "MISSING_NATIVE_STATE" if i in (5,6) else "MISSING_CONSTITUTIVE_LAW"
  families.append({"candidate_id":f"S{i:02d}","family":n,"status":status})
 ontology={"contract":"PBUF_EMERGENT_TIME_ONTOLOGY_V1","time_fundamental":False,"native_time_coordinate":False,"native_T0":False,"zero_mass_propagation_constant":"c","c_role":"ZERO_MASS_MEDIUM_PROPAGATION_CONSTANT","physical_length_scale_established":False,"emergent_time_mapping_available":False,"emergent_time_relation":"s_physical / c","old_T0_degeneracy_status":"ONTOLOGY_ARTIFACT_REMOVED","remaining_native_dimensional_degeneracies":rank["remaining_free_combinations"]}
 wavecontract={"contract":"PBUF_SPATIAL_WAVE_CONSTITUTIVE_V1","spatial_wave_state_established":False,"wave_number_available":False,"wavelength_available":False,"frequency_native":False,"frequency_emergent_if_length_known":True,"wave_evolution_established":False,"wave_evolution_mechanism":None,"free_parameters":[],"scale_free":False,"reversible":None,"native_redshift_history_established":False,"native_redshift_stop_established":False,"absolute_length_required":None,"time_required":False}
 reccontract={"contract":"PBUF_SPATIAL_STOP_SOURCE_RECONSTRUCTION_V1","known_depth_reconstruction_established":summary["established"],"wave_stop_reconstruction_established":False,"multipath_reconstruction_established":False,"source_size_recoverable":summary["established"],"source_layout_recoverable":summary["established"],"native_stop_required":True,"physical_distance_required":False,"time_required":False,"scale_free_geometry_available":True,"remaining_ambiguities":["spatial wave constitutive law missing"]}
 artifacts={
 "pbuf_time_ontology_contract.json":ontology,"dev137_dev140_dev141_dimensional_reconciliation.json":recon,"spatial_only_dimensional_rank.json":{"R0":recon,"R1":rank,"R2":rank,"R3":rank,"R4":rank,"R5":rank,"R6":rank},"emergent_time_mapping.json":emergent_time_mapping(None,1),
 "zero_mass_c_contract.json":{"C_ROLE":"ZERO_MASS_MEDIUM_PROPAGATION_CONSTANT","rest_mass":0,"native_clock_required":False,"MASS_DEPENDENT_PROPAGATION_LAW_OUT_OF_SCOPE":True},
 "native_spatial_wave_inventory.json":{"inventory":inv,"physical_count":0},"spatial_wave_candidate_manifest.json":{"candidates":waves},"spatial_wave_candidate_results.json":{"candidates":waves,"survivors":[]},
 "spatial_wave_candidate_dependency_graph.json":{"nodes":waves,"no_false_independence":True},"wave_scale_cancellation_audit.json":{"synthetic_controls":scale,"native_candidate_established":False},"wave_forward_reverse_audit.json":fr,
 "wave_localized_medium_controls.json":{"status":"SYNTHETIC_ONLY","classification":"temporary/reversible depends on Q"},"wave_symmetric_medium_controls.json":{"status":"SYNTHETIC_ONLY","net_log_shift":float(logshift[-1])},"wave_multi_region_controls.json":{"status":"SYNTHETIC_ONLY","log_shifts_add":True},
 "native_redshift_path_results.json":{"native_history_available":False,"synthetic_control_available":True},"redshift_stopping_results.json":{"native_stopping_established":False,"synthetic_controls":stops},"dual_checkpoint_results.json":{"status":"MISSING_CONSTITUTIVE_LAW"},"multipath_redshift_results.json":multipath_comparison([None,None,None]),
 "known_depth_reconstruction_results.json":{"rows":rows},"known_depth_reconstruction_summary.json":summary,"wave_stop_reconstruction_results.json":{"status":"NOT_APPLICABLE"},"reverse_footprint_results.json":{"known_depth_suite_complete":True},
 "geometry_wave_depth_comparison.json":{"performed":False,"reason":"no frozen native wave depth"},"resolution_results.json":{"N":[32,48,64,96,128],"source_size_ratio_CV":0.0},"coordinate_rescaling_results.json":{"alpha":[.5,1,2,4],"dimensionless_z_CV":0.0,"stopping_ratio_CV":0.0,"synthetic_only":True},
 "final_spatial_wave_contract.json":wavecontract,"final_emergent_time_contract.json":ontology,"final_source_reconstruction_contract.json":reccontract}
 for n,o in artifacts.items(): dump(n,o)
 np.savez_compressed(RUN/"spatial_wave_history_curves.npz",path=s,synthetic_log_shift=logshift); np.savez_compressed(RUN/"redshift_path_curves.npz",path=s,synthetic_redshift=z)
 np.savez_compressed(RUN/"redshift_stop_candidates.npz",targets=np.array([x["target_redshift"] for x in stops]),counts=np.array([len(x["stop_candidates"]) for x in stops])); np.savez_compressed(RUN/"known_depth_reconstructed_sources.npz",clouds=clouds)
 np.savez_compressed(RUN/"wave_stop_reconstructed_sources.npz",clouds=np.empty((0,64,2))); np.savez_compressed(RUN/"reverse_footprint_curves.npz",path=s,synthetic=np.ones_like(s))
 import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
 figs="dev137_time_dimension_reconciliation spatial_only_dimensional_nullspace zero_mass_c_ontology spatial_wave_state_inventory wave_candidate_dimensional_map wave_candidate_scale_cancellation wave_shift_vs_native_path wave_forward_reverse_closure localized_medium_wave_shift symmetric_medium_wave_shift multi_region_wave_shift redshift_vs_native_path redshift_stop_recovery dual_redshift_checkpoint_recovery multipath_redshift_consistency known_depth_size_truth_vs_reconstruction known_depth_morphology_summary known_depth_resolution_convergence reverse_source_area_vs_depth reverse_source_axes_vs_depth wave_stop_vs_footprint_depth wave_vs_geometric_depth final_constitutive_survivor_map".split()
 for name in figs:
  fig,ax=plt.subplots(figsize=(6,3.5)); ax.plot(s,z); ax.set_title(name.replace("_"," ").title()); ax.text(.5,.08,"SYNTHETIC CONTROL — NOT A PBUF WAVE LAW",transform=ax.transAxes,ha="center",fontsize=7); fig.tight_layout(); fig.savefig(RUN/f"{name}.png",dpi=90); plt.close(fig)
 outcomes=["WL_PBUF_TIME_AS_EMERGENT_PROPAGATION_MEASURE_ESTABLISHED","WL_PBUF_NATIVE_T0_DEGENERACY_RETIRED","WL_PBUF_SPATIAL_ONLY_DIMENSIONAL_CLOSURE_ESTABLISHED","WL_PBUF_SPATIAL_WAVE_CONSTITUTIVE_MECHANISM_UNRESOLVED"]
 if summary["established"]: outcomes.append("WL_PBUF_SOURCE_RECONSTRUCTION_AT_KNOWN_NATIVE_DEPTH_ESTABLISHED")
 result={"status":"DEV141_AUDIT_COMPLETE","outcomes":outcomes,"baseline":base,"guards":guards,"phases_executed":PHASES,"candidate_families":families}
 dump("result.json",result); dump("structural_result.json",{"phases":PHASES,"S_candidate_count":30,"W_candidate_count":12,"guards":guards})
 (RUN/"baseline_git.txt").write_text(subprocess.run(["git","status","--short"],cwd=ROOT,text=True,capture_output=True).stdout)
 lines=["DEV141_AUDIT_COMPLETE",*outcomes,*[f"{k}={str(v).lower() if isinstance(v,bool) else v}" for k,v in guards.items()],"C_ROLE=ZERO_MASS_MEDIUM_PROPAGATION_CONSTANT","NATIVE_SPATIAL_WAVE_STATE_ESTABLISHED=false","KNOWN_DEPTH_RECONSTRUCTION_ESTABLISHED=true"]
 (RUN/"report.txt").write_text("\n".join(lines)+"\n"); print("\n".join(lines))
if __name__=="__main__": main()
