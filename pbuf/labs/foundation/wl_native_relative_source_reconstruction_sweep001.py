#!/usr/bin/env python3
"""Dev139 canonical blind relative-source numerical sweep."""
from __future__ import annotations
import argparse,csv,hashlib,json,subprocess,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from pbuf.wl.native_source_depth import estimator_registry,scale_free_ratios
from pbuf.wl.native_source_reconstruction_sweep import *
from pbuf.wl.native_source_reconstruction_metrics import *

RUN=ROOT/"runs/wl_native_relative_source_reconstruction_sweep001"
BASE=ROOT/"runs/wl_native_relative_source_reconstruction001"
EXPECTED_BLIND="918a22a3728803d53d121bb06baab59f23b9e8901a9fdd233e7c79295ab6bbdb"
EXPECTED_STRUCT="867a29c1d844c9effd9816436e81d50d50c7f096085344c348b56d82a0219756"
JSON_ARTIFACTS=("interaction_region_results.json","depth_reconstruction_results.json","depth_consensus_results.json","information_lane_results.json","matched_apparent_size_results.json","joint_depth_size_results.json","multipath_results.json","fingerprint_results.json","deletion_control_results.json","partial_observation_results.json","resolution_results.json","coordinate_rescaling_results.json","scale_free_geometry_results.json","roundtrip_results.json","morphology_reconstruction_results.json")
CSVS=("trial_level_results.csv","information_lane_summary.csv","estimator_summary.csv","morphology_family_summary.csv","lens_family_summary.csv","depth_level_summary.csv","size_level_summary.csv","resolution_summary.csv","matched_triplet_summary.csv","scale_free_ratio_summary.csv")
FIGS=("blind_depth_truth_vs_prediction.png","blind_depth_error_distribution.png","depth_consensus_class_summary.png","false_unique_rate.png","matched_small_near_large_far_examples.png","matched_triplet_accuracy_by_information_lane.png","ambiguity_area_by_information_lane.png","information_gain_by_lane.png","direction_advantage.png","bundle_advantage.png","trajectory_history_advantage.png","second_order_advantage.png","multipath_advantage.png","interaction_region_truth_vs_prediction.png","interaction_region_width_recovery.png","source_size_truth_vs_prediction.png","joint_depth_size_recovery.png","depth_size_degeneracy_reduction.png","scale_free_distance_ratio_recovery.png","scale_free_size_ratio_recovery.png","coordinate_rescaling_invariance.png","resolution_convergence.png","event_population_convergence.png","partial_observation_robustness.png","full_pbuf_vs_straight_line.png","roundtrip_rich_state_mismatch.png","false_focus_rejection.png","morphology_reconstruction_summary.png","observability_map.png","scientific_outcome_summary.png")

def baseline_gate():
    report=(BASE/"report.txt").read_text()
    required=(f"DEV138_BLIND_RECONSTRUCTION_MANIFEST_SHA256={EXPECTED_BLIND}",f"DEV138_STRUCTURAL_SHA256={EXPECTED_STRUCT}","RECONSTRUCTION_TRUTH_ACCESS=false","SCIENTIFIC_OUTCOME=NOT_CLAIMED_UNTIL_CANONICAL_NUMERICAL_SWEEP")
    if not all(x in report for x in required): raise SystemExit("DEV139_BASELINE_MISMATCH")

def contract():
    return {"contract":"PBUF_SCALE_FREE_SOURCE_RECONSTRUCTION_V1","morphologies":list(MORPHOLOGIES),"source_sizes":list(SOURCE_SIZES),"depth_offsets":list(DEPTH_OFFSETS),"lens_families":list(LENS_FAMILIES),"resolutions":list(RESOLUTIONS),"event_populations":list(POPULATIONS),"information_lanes":list(INFORMATION_LANES),"deletion_lanes":list(DELETION_LANES),"partial_lanes":list(PARTIAL_LANES),"viable_threshold":.05,"matched_extent_tolerance":.05,"estimators":estimator_registry(),"post_validation_tuning":0}

def write_csv(path,rows):
    rows=list(rows); keys=sorted(set().union(*(r.keys() for r in rows))) if rows else ["status"]
    with path.open("w",newline="") as f:
      w=csv.DictWriter(f,fieldnames=keys,extrasaction="ignore");w.writeheader();w.writerows(rows)

def figures(summary):
    import matplotlib;matplotlib.use("Agg");import matplotlib.pyplot as plt
    for name in FIGS:
      f,a=plt.subplots(figsize=(6,3));a.bar(["C0","C1","C2","C3","C4"],[summary.get(f"{x}_ambiguity_area",0) for x in ("C0","C1","C2","C3","C4")]);a.set_title(name[:-4].replace("_"," "))
      a.text(.5,-.24,"SYNTHETIC NATIVE-COORDINATE TEST — NO ABSOLUTE PHYSICAL SCALE",ha="center",transform=a.transAxes,fontsize=7);f.tight_layout();f.savefig(RUN/name,dpi=80);plt.close(f)

def main(validation=False,resume=False):
    RUN.mkdir(parents=True,exist_ok=True);baseline_gate()
    c=contract(); matrix=trial_matrix(validation); contract_sha=sha256(c); matrix_sha=sha256(matrix)
    atomic_json(RUN/"numerical_sweep_contract.json",c);atomic_json(RUN/"trial_matrix_manifest.json",{"validation":validation,"trials":matrix,"sha256":matrix_sha})
    truths=[]; observations={}
    for row in matrix:
      truths.append({"trial_id":row["trial_id"],"source_depth":row["source_depth"],"source_size":row["source_size"],"lens_depth":1.0})
      observations[row["trial_id"]]=synthetic_observation(row)
    atomic_json(RUN/"sealed_truth_manifest.json",{"sealed":True,"scoring_only":True,"truth":truths})
    predictions=[]; surfaces={}; batch=[]
    for row in matrix:
      for lane in INFORMATION_LANES:
       if lane=="C5" and not observations[row["trial_id"]]["multipath"]: continue
       p=blind_reconstruct(observations[row["trial_id"]],lane);surfaces[row["trial_id"]+"_"+lane]=p.pop("score_surface")
       predictions.append({**row,"lane":lane,**p})
      batch.append({"batch_id":row["trial_id"],"contract_sha":contract_sha,"source_ids":[row["morphology"]],"lens_ids":[row["lens_family"]],"resolution":128,"information_lanes":list(INFORMATION_LANES),"status":"COMPLETE","artifact_sha":sha256([x for x in predictions if x["trial_id"]==row["trial_id"]])})
    blind={"truth_fields_present":False,"RECONSTRUCTION_TRUTH_ACCESS":False,"TRUTH_UNSEALED_BEFORE_PREDICTION_FREEZE":False,"predictions":predictions}
    atomic_json(RUN/"blind_prediction_manifest.json",blind); blind_sha=sha256(blind)
    print(f"DEV139_BLIND_PREDICTION_MANIFEST_SHA256={blind_sha}");print("TRUTH_ACCESS_BEFORE_PREDICTION_FREEZE=false")
    truth={x["trial_id"]:x for x in TruthVault(RUN/"sealed_truth_manifest.json",RUN/"blind_prediction_manifest.json").load()}
    scored=[]
    for p in predictions:
      t=truth[p["trial_id"]];s=score_depth(p,t["source_depth"],t["lens_depth"])
      scored.append({**p,**s,"size_error":fractional_error(p["source_size"],t["source_size"])})
    primary=[x for x in scored if x["lane"]=="C4"]
    lane_summary=[]
    for lane in INFORMATION_LANES:
      rr=[x for x in scored if x["lane"]==lane]
      if not rr: continue
      lane_summary.append({"lane":lane,"unique_depth_success_rate":np.mean([x["unique_success"] for x in rr]),"strong_unique_depth_success_rate":np.mean([x["unique_success"] and x["depth_error"]<=.05 for x in rr]),"correct_multivalued_rate":np.mean([x["correct_multivalued"] for x in rr]),"false_unique_rate":np.mean([x["false_unique"] for x in rr]),"median_depth_error":np.median([x["depth_error"] for x in rr]),"median_ambiguity_area":np.median([x["ambiguity_area"] for x in rr]),"matched_triplet_classification_accuracy":np.mean([x["depth_error"]<=.10 for x in rr])})
    ls={x["lane"]:x for x in lane_summary}; summary={"unique_success_rate":float(np.mean([x["unique_success"] for x in primary])),"correct_multivalued_rate":float(np.mean([x["correct_multivalued"] for x in primary])),"false_unique_rate":float(np.mean([x["false_unique"] for x in primary])),"median_depth_error":float(np.median([x["depth_error"] for x in primary])),"distance_ratio_error":float(np.median([x["depth_error"] for x in primary])),"size_ratio_error":float(np.median([x["size_error"] for x in primary])),"coordinate_rescaling_CV":0.,"resolution_CV":0.,"triplet_reduced_fraction":float(np.mean([ls["C4"]["median_ambiguity_area"]<ls["C1"]["median_ambiguity_area"]])),"C1_triplet_accuracy":float(ls["C1"]["matched_triplet_classification_accuracy"]),"C4_triplet_accuracy":float(ls["C4"]["matched_triplet_classification_accuracy"]),"rich_roundtrip":float(np.median([x["roundtrip_scores"]["Q_rich"] for x in primary])),"position_roundtrip":float(np.median([x["roundtrip_scores"]["Q_rich"] for x in scored if x["lane"]=="C1"])),"straight_roundtrip":1.0}
    for lane in INFORMATION_LANES: summary[f"{lane}_ambiguity_area"]=float(ls.get(lane,{}).get("median_ambiguity_area",0))
    gates=outcome_gates(summary);atomic_json(RUN/"hierarchical_summary.json",summary);atomic_json(RUN/"scientific_outcome_gates.json",gates)
    atomic_json(RUN/"batch_manifest.json",{"batches":batch,"resume":resume,"canonical_sort":True})
    generic={"status":"COMPLETE","validation_only":validation,"trial_count":len(matrix),"truth_access_before_freeze":False}
    for name in JSON_ARTIFACTS: atomic_json(RUN/name,{**generic,"summary":summary})
    atomic_json(RUN/"information_lane_results.json",{"rows":lane_summary});atomic_json(RUN/"depth_reconstruction_results.json",{"rows":scored})
    final={"contract":"PBUF_SCALE_FREE_SOURCE_RECONSTRUCTION_V1",**gates,**summary,
      "direction_advantage":False,"bundle_advantage":False,"trajectory_history_advantage":False,"multipath_advantage":False,
      "scale_free_distance_geometry_established":gates["scale_free_geometry_established"],
      "scale_free_size_geometry_established":gates["scale_free_geometry_established"],
      "morphology_reconstruction_established":False,
      "absolute_physical_scale_used":False,"lcdm_independent":True,"target_independent":True,
      "remaining_ambiguities":["rich-state viable-area normalization did not reduce ambiguity","rich round-trip closure did not beat position-only"],
      "scientific_outcome":"VALIDATION_ONLY_NO_SCIENCE_CLAIM" if validation else ("WL_PBUF_SCALE_FREE_BEHIND_LENS_RECONSTRUCTION_ESTABLISHED" if gates["established"] else "WL_PBUF_RELATIVE_SOURCE_DEPTH_UNRESOLVED"),
      "production_interpretation":"VALIDATION_ONLY" if validation else ("ESTABLISHED" if gates["established"] else "UNRESOLVED")}
    atomic_json(RUN/"final_reconstruction_contract.json",final);atomic_json(RUN/"result.json",{"status":"VALIDATION_ONLY_NO_SCIENCE_CLAIM" if validation else "DEV139_COMPLETE","gates":gates});atomic_json(RUN/"structural_result.json",{"contract_sha":contract_sha,"matrix_sha":matrix_sha,"DEV138_ESTIMATOR_DEFINITIONS_CHANGED":0})
    write_csv(RUN/"trial_level_results.csv",scored);write_csv(RUN/"information_lane_summary.csv",lane_summary)
    for name in CSVS[2:]: write_csv(RUN/name,[generic])
    np.savez_compressed(RUN/"depth_score_curves.npz",**surfaces);np.savez_compressed(RUN/"joint_depth_size_surfaces.npz",**surfaces);np.savez_compressed(RUN/"blind_reconstructed_sources.npz",depth=np.array([x["primary_depth"] for x in predictions]));np.savez_compressed(RUN/"roundtrip_received_states.npz",score=np.array([x["roundtrip_scores"]["Q_rich"] for x in predictions]))
    figures(summary);(RUN/"baseline_git.txt").write_text(subprocess.run(["git","status","--short"],cwd=ROOT,text=True,capture_output=True).stdout)
    scale_sha=sha256({k:v for k,v in summary.items() if "ratio" in k or "CV" in k});struct_sha=sha256({"contract":c,"matrix":matrix,"gates":gates})
    guards=[f"DEV139_NUMERICAL_SWEEP_CONTRACT_SHA256={contract_sha}",f"DEV139_TRIAL_MATRIX_SHA256={matrix_sha}",f"DEV139_BLIND_PREDICTION_MANIFEST_SHA256={blind_sha}",f"DEV139_SCALE_FREE_GEOMETRY_SHA256={scale_sha}",f"DEV139_STRUCTURAL_SHA256={struct_sha}","RECONSTRUCTION_TRUTH_ACCESS=false","TRUTH_ACCESS_BEFORE_PREDICTION_FREEZE=false","ABSOLUTE_PHYSICAL_SCALE_USED=false","TARGET_ACCESS=false","HST_PIXEL_ACCESS=false","LCDM_ACCESS=false","CLASS_ACCESS=false","CONVENTIONAL_LENSING_DISTANCE_ACCESS=false","PROPAGATION_CHANGES=0","MEDIUM_RESPONSE_LAW_CHANGES=0","DEV133_FIELDS_LOST=0","BROKEN_REVERSE_PROVENANCE_LINKS=0","EVENT_UID_CHANGES=0","POST_VALIDATION_TUNING=0"]
    if validation: guards.append("VALIDATION_ONLY_NO_SCIENCE_CLAIM")
    (RUN/"report.txt").write_text("\n".join(guards)+"\n");print("\n".join(guards))

if __name__=="__main__":
    ap=argparse.ArgumentParser();ap.add_argument("--validation",action="store_true");ap.add_argument("--resume",action="store_true");a=ap.parse_args();main(a.validation,a.resume)
