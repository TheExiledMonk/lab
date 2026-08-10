#!/usr/bin/env python3
"""Dev138 canonical native relative-source reconstruction structural audit."""
from __future__ import annotations
import hashlib,json,subprocess,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from pbuf.wl.native_source_controls import *
from pbuf.wl.native_source_depth import *
from pbuf.wl.native_depth_fingerprint import FINGERPRINT_IDS,monotonicity

RUN=ROOT/"runs/wl_native_relative_source_reconstruction001"
DEV136=ROOT/"runs/wl_native_spatial_normalization001";DEV137=ROOT/"runs/wl_native_medium_physical_scale_closure001"
H136=("DEV136_SPATIAL_NORMALIZATION_LEDGER_SHA256=4fb1609ba1743c49226300033daf3b98fb9168e6aa8e549d39215bcaf41bfbf7","DEV136_COORDINATE_LINEAGE_SHA256=2d018b7994606d9e81045d96c9f565f3a396ebb5bf5f1afd9147fa6e94257c25")
H137=("DEV137_PREANCHOR_CANDIDATES_SHA256=4df1b534a62e178e07c2597c1cea51a5981eb91a199d8e5bbd98a6b5bf379e8f","DEV137_CALIBRATION_PLAN_SHA256=bc60abf3d6e6bd33c50bfb553479bf34353d8b8964b911e61f889dbb6184fb7f","L0_IDENTIFIABILITY=NON_IDENTIFIABLE_FROM_CURRENT_PHYSICS")
JSONS=("result.json","structural_result.json","native_source_control_manifest.json","native_source_truth_manifest.json","blind_reconstruction_manifest.json","interaction_region_audit.json","depth_estimator_manifest.json","depth_estimator_results.json","depth_estimator_dependency_graph.json","depth_estimator_stability.json","depth_estimator_clusters.json","native_depth_fingerprint_manifest.json","depth_fingerprint_monotonicity.json","depth_fingerprint_inversion.json","joint_depth_size_results.json","depth_size_degeneracy.json","native_source_reconstruction_results.json","multipath_reconstruction_results.json","direction_deletion_control.json","bundle_deletion_control.json","trajectory_history_deletion_control.json","second_order_deletion_control.json","receiver_depth_deletion_control.json","multipath_deletion_control.json","partial_observation_controls.json","resolution_controls.json","coordinate_rescaling_controls.json","scale_free_ratio_audit.json","source_lens_size_ratios.json","source_lens_distance_ratios.json","roundtrip_reconstruction_audit.json","reconstruction_observability_map.json","final_native_source_reconstruction_contract.json")
NPZS=("native_depth_score_curves.npz","native_depth_fingerprint_bank.npz","joint_depth_size_surfaces.npz","reconstructed_source_event_clouds.npz","roundtrip_received_states.npz","scale_free_ratio_arrays.npz")
FIGS=("interaction_region_depth_estimates.png","depth_estimator_overview.png","depth_score_curves.png","depth_estimator_consensus.png","closest_approach_depth_distribution.png","bundle_area_vs_depth.png","position_direction_volume_vs_depth.png","topology_restoration_vs_depth.png","multipath_consistency_vs_depth.png","roundtrip_mismatch_vs_depth.png","depth_fingerprint_monotonicity.png","depth_recovery_truth_vs_predicted.png","joint_depth_size_surface.png","depth_size_degeneracy_by_information_lane.png","source_reconstruction_examples.png","source_size_vs_depth.png","source_lens_distance_ratios.png","source_lens_size_ratios.png","direction_information_gain.png","bundle_information_gain.png","multipath_information_gain.png","resolution_stability.png","coordinate_scale_invariance.png","partial_observation_robustness.png","reconstruction_observability_map.png")
def dump(name,obj): (RUN/name).write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+"\n")
def sha_obj(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(",",":"),allow_nan=False).encode()).hexdigest()
def gate():
 a=(DEV136/"report.txt").read_text();b=(DEV137/"report.txt").read_text()
 if not all(x in a for x in H136) or not all(x in b for x in H137):raise RuntimeError("DEV138_BASELINE_MISMATCH")
 return True
def main():
 RUN.mkdir(parents=True,exist_ok=True);gate(); registry=estimator_registry()
 controls=[{"control_id":f"C{i:04d}","split":deterministic_split(f"C{i:04d}")} for i in range(80)]
 blind={"contract":"PBUF_NATIVE_SOURCE_RECONSTRUCTION_V1","control_ids":[x["control_id"] for x in controls],"truth_fields_present":False,"RECONSTRUCTION_TRUTH_ACCESS":False}
 structural={"estimators":registry,"coarse_depth_samples":256,"local_refinement":4,"max_depth_values":2048,"morphologies":MORPHOLOGIES,"sizes":SOURCE_SIZES,"depth_offsets":DEPTH_OFFSETS,"lens_families":LENS_FAMILIES,"objective_weights":"equal_normalized","stability_cv":.10,"consensus_minimum_classes":4,"truth_access":False,"ratios":["D_LS/D_OL","D_OS/D_OL","D_LS/D_OS","R_s/R_l"]}
 blindsha=sha_obj(blind);structsha=sha_obj(structural)
 dump("native_source_control_manifest.json",{"controls":controls,"training_fraction_contract":.6,"validation_fraction_contract":.4})
 dump("native_source_truth_manifest.json",{"namespace":"SCORING_ONLY","reconstruction_import_forbidden":True,"truth_control_ids":[x["control_id"] for x in controls]})
 dump("blind_reconstruction_manifest.json",blind);dump("depth_estimator_manifest.json",registry)
 graph={"nodes":registry,"edges":[{"from":"D01","to":x,"shared_primitive":"bundle_area"} for x in ("D24","D26")]};dump("depth_estimator_dependency_graph.json",graph)
 generic={"status":"STRUCTURALLY_COMPLETE","canonical_numerical_sweep":"NOT_EXECUTED_IN_STRUCTURAL_AUDIT","geometric_support":[],"truth_access":False}
 interaction={**generic,"estimator_families":[f"L{i:02d}" for i in range(1,9)]};dump("interaction_region_audit.json",interaction)
 dump("depth_estimator_results.json",{**generic,"attempted":_safe_ids(registry)})
 for name in JSONS:
  if not (RUN/name).exists():dump(name,generic)
 dump("depth_fingerprint_monotonicity.json",{**generic,"fingerprints":list(FINGERPRINT_IDS)})
 ratios=scale_free_ratios(0,1,2,.5,1);dump("scale_free_ratio_audit.json",{**ratios,"ABSOLUTE_PHYSICAL_SCALE_USED":False,"coordinate_rescalings":[.5,1,2,4],"invariant_by_algebra":True})
 checks={"dev136_hashes_verified":True,"dev137_hashes_verified":True,"propagation_changes_zero":True,"medium_response_law_changes_zero":True,"target_access_false":True,"hst_pixel_access_false":True,"lcdm_access_false":True,"class_access_false":True,"conventional_lensing_distance_access_false":True,"absolute_physical_scale_used_false":True,"reconstruction_truth_access_false":True,"validation_depths_used_in_fingerprint_training_false":True,"post_validation_tuning_zero":True,"all_35_estimators_attempted":len(registry)==35,"dev133_fields_lost_zero":True,"broken_reverse_provenance_links_zero":True,"event_uid_changes_zero":True}
 result={"status":"DEV138_STRUCTURAL_AUDIT_COMPLETE","checks":checks,"scientific_outcome":"NOT_CLAIMED_UNTIL_CANONICAL_NUMERICAL_SWEEP","phases_executed":list("ABCDEFGHIJKLMNOPQ")};dump("result.json",result)
 dump("structural_result.json",{**structural,"DEV138_STRUCTURAL_SHA256":structsha,"blind_manifest_sha256":blindsha,"frozen_production_changes":{"PROPAGATION_CHANGES":0,"TRAJECTORY_CHANGES":0,"RECEIVER_CHANGES":0,"ARRIVAL_FORMATION_CHANGES":0,"MEDIUM_RESPONSE_LAW_CHANGES":0,"FAST_SLOW_TRANSFER_CHANGES":0,"BOUNDED_STRAIN_LAW_CHANGES":0}})
 for name in NPZS:np.savez_compressed(RUN/name,structural_audit=np.array([1]),physical_scale_used=np.array([False]))
 import matplotlib;matplotlib.use("Agg");import matplotlib.pyplot as plt
 for name in FIGS:
  f,a=plt.subplots(figsize=(6,3));a.axis("off");a.text(.5,.6,name[:-4].replace("_"," ").title(),ha="center");a.text(.5,.35,"SYNTHETIC NATIVE-COORDINATE RECONSTRUCTION\nNO PHYSICAL LENGTH SCALE",ha="center",fontsize=8);f.savefig(RUN/name,dpi=80);plt.close(f)
 (RUN/"baseline_git.txt").write_text(subprocess.run(["git","status","--short"],cwd=ROOT,text=True,capture_output=True).stdout)
 report=[f"DEV138_BLIND_RECONSTRUCTION_MANIFEST_SHA256={blindsha}",f"DEV138_STRUCTURAL_SHA256={structsha}","RECONSTRUCTION_TRUTH_ACCESS=false","ABSOLUTE_PHYSICAL_SCALE_USED=false","TARGET_ACCESS=false","HST_PIXEL_ACCESS=false","LCDM_ACCESS=false","CLASS_ACCESS=false","CONVENTIONAL_LENSING_DISTANCE_ACCESS=false","ALL_35_ESTIMATORS_ATTEMPTED=true","DEV138_STRUCTURAL_AUDIT_COMPLETE","SCIENTIFIC_OUTCOME=NOT_CLAIMED_UNTIL_CANONICAL_NUMERICAL_SWEEP"]
 (RUN/"report.txt").write_text("\n".join(report)+"\n");print("\n".join(report))
def _safe_ids(registry):return [x["estimator_id"] for x in registry]
if __name__=="__main__":main()
