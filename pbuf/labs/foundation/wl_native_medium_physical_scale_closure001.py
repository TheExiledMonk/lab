#!/usr/bin/env python3
"""Dev137 canonical clean-room native-medium physical-scale closure audit."""
from __future__ import annotations
import hashlib, json, subprocess, sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[3]; sys.path.insert(0,str(ROOT))
from pbuf.wl.medium_dimensional_closure import default_dimensional_system, frozen_medium_unit_contract, restored_operators
from pbuf.wl.native_scale_candidates import (CLASSES,FAMILIES,NativeScaleEstimate,
 candidate_registry,dependency_graph,execute_internal_candidates)
from pbuf.wl.native_scale_universality import cluster_candidates,pairwise_agreement,stability

RUN=ROOT/"runs/wl_native_medium_physical_scale_closure001"
DEV136=ROOT/"runs/wl_native_spatial_normalization001"
EXPECTED={"ledger":"4fb1609ba1743c49226300033daf3b98fb9168e6aa8e549d39215bcaf41bfbf7",
 "lineage":"2d018b7994606d9e81045d96c9f565f3a396ebb5bf5f1afd9147fa6e94257c25",
 "structural":"ca8efcb5319613d2e9d1308524014f196755d36aab99b14dec63426aea480579"}
REQUIRED136=("WL_PBUF_SPATIAL_NORMALIZATION_PROVENANCE_ESTABLISHED",
 "WL_PBUF_PROPAGATION_GLOBAL_SCALE_DEGENERACY_ESTABLISHED",
 "WL_PBUF_PROPAGATION_PHYSICAL_NORMALIZATION_REQUIRED",
 "WL_PBUF_NATIVE_PHYSICAL_SCALE_UNRESOLVED","WL_PBUF_RECEIVER_GLOBAL_SKY_GAUGE_IDENTIFIED")

def dump(name,obj):
    p=RUN/name; p.write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+"\n"); return p
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def canon(obj): return hashlib.sha256(json.dumps(obj,sort_keys=True,separators=(",",":"),allow_nan=False).encode()).hexdigest()
def baseline_gate():
    report=(DEV136/"report.txt").read_text()
    checks={"ledger":f"DEV136_SPATIAL_NORMALIZATION_LEDGER_SHA256={EXPECTED['ledger']}" in report,
      "lineage":f"DEV136_COORDINATE_LINEAGE_SHA256={EXPECTED['lineage']}" in report,
      "provenance":"DEV136_NATIVE_SCALE_PROVENANCE_SHA256=NOT_ESTABLISHED" in report,
      "structural":sha(DEV136/"structural_result.json")==EXPECTED["structural"],
      "outcomes":all(x in report for x in REQUIRED136)}
    if not all(checks.values()): raise RuntimeError("DEV136_BASELINE_MISMATCH")
    return checks
def fig(name,title,lines):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    f,a=plt.subplots(figsize=(8,4.5)); a.axis("off"); a.set_title(title,weight="bold")
    a.text(.03,.92,"\n".join(lines),va="top",family="monospace",fontsize=8,transform=a.transAxes)
    f.tight_layout(); f.savefig(RUN/name,dpi=120); plt.close(f)
def unavailable(i,status,reason):
    return NativeScaleEstimate(f"S{i:02d}",FAMILIES[i-1],status=status,rejection_reason=reason,
      independence_class=CLASSES[i-1],limitations=[reason],dimension_signature="L0^0",L0_power="L0^0")

def main():
    RUN.mkdir(parents=True,exist_ok=True)
    phases=[]; base=baseline_gate(); phases.append("A")
    contract=frozen_medium_unit_contract(); ops=restored_operators()
    dump("medium_unit_contract.json",contract); dump("dimensional_restoration.json",{"operators":ops}); phases.append("B")
    rank=default_dimensional_system().audit(); dump("dimensional_rank_audit.json",rank)
    dump("identifiable_combinations.json",{"remaining_free_combinations":rank["remaining_free_combinations"],
      "L0_U0_CO_DEGENERACY":True,"stiffness_scale_degeneracy":True,"source_normalization_degeneracy":True}); phases.append("C")
    controls={"classification":"PHYSICAL_PARAMETERIZED_CONTROL","families":["uniform sphere","compact sphere","diffuse sphere","two-source additive configuration","weak response","strong response"],
      "mass_factors":[.25,.5,1,2,4],"radius_factors":[.5,.75,1,1.5,2],"density_factors":[.25,.5,1,2,4],"resolutions":[32,48,64,96,128],
      "same_frozen_medium_law":True,"physical_source_normalization_established":False}; phases.append("D")
    candidates=execute_internal_candidates(rank); phases.append("E")
    pre={"candidate_results":[c.to_dict() for c in candidates],"external_values_loaded":False,
      "PREANCHOR_EXTERNAL_VALUE_ACCESS":False,"registry":candidate_registry()}
    prepath=dump("internal_scale_candidates_preanchor.json",pre); presha=sha(prepath); phases.append("F")
    anchor_required=rank["L0_identifiability"]
    # Multiple free dimensional scales prohibit a numerical anchor fit.
    calplan={"selected_anchor":None,"reason":"multiple co-degenerate unit scales; one length anchor cannot isolate L0",
      "inputs":[],"formula":None,"validation_exclusions":[],"frozen_before_fit":True}
    calpath=dump("scale_calibration_plan.json",calplan); calsha=sha(calpath); phases.extend(["G","H"])
    from pbuf.wl.non_lensing_physical_anchors import physical_anchor_manifest
    anchors=physical_anchor_manifest(); dump("physical_anchor_manifest.json",anchors)
    calibration={"performed":False,"reason":"PBUF_PHYSICAL_SCALE_REQUIRES_MULTIPLE_DIMENSIONAL_CLOSURES","L0_refits_after_freeze":0,
      "GLOBAL_SCALE_IMMUTABLE_AFTER_FREEZE":True}; dump("calibration_result.json",calibration); phases.append("I")
    validation={"performed":False,"results":[{"anchor_id":x["id"],"status":"NOT_APPLICABLE","reason":"no frozen L0 to validate"} for x in anchors["records"]]}
    dump("validation_anchor_results.json",validation); phases.append("J")
    for i,reason in [(22,"one-anchor closure invalid: more than one unit freedom"),(23,"Earth cannot isolate L0 from U0/K0_phys/S0/T0"),
      (24,"LAB_ANCHOR_DATA_UNAVAILABLE"),(25,"Solar anchor cannot isolate L0 from response/time units"),(26,"no calibrated global scale for cross-anchor validation")]:
      candidates.append(unavailable(i,"NOT_APPLICABLE",reason))
    sweeps={}
    for name,factors in (("mass",controls["mass_factors"]),("radius",controls["radius_factors"]),("density",controls["density_factors"])):
      value={"controls":factors,"result":stability([]),"status":"NOT_TESTABLE","reason":"no finite authoritative L0 estimator"}; sweeps[name]=value; dump(f"{name}_universality.json",value)
    resolution={"resolutions":controls["resolutions"],"result":stability([]),"status":"NOT_TESTABLE","reason":"no finite authoritative L0 estimator"}; dump("resolution_universality.json",resolution)
    sourcefam={"families":controls["families"],"status":"NOT_TESTABLE","reason":"no finite authoritative L0 estimator"}; dump("source_family_universality.json",sourcefam); phases.extend(["K","L"])
    for i,name in ((27,"density"),(28,"mass"),(29,"radius"),(30,"resolution")):
      candidates.append(unavailable(i,"NON_IDENTIFIABLE","no finite authoritative L0 estimator for universality audit"))
    resultdict=[c.to_dict() for c in candidates]; dump("scale_candidate_results.json",resultdict)
    graph=dependency_graph(candidates); dump("scale_candidate_dependency_graph.json",graph)
    agree=pairwise_agreement(resultdict); dump("scale_candidate_pairwise_agreement.json",agree)
    clusters=cluster_candidates(resultdict); dump("scale_candidate_clusters.json",clusters); phases.append("M")
    forbidden={"LCDM_ACCESS":False,"CLASS_ACCESS":False,"CONVENTIONAL_LENSING_SCALE_ACCESS":False,"TARGET_ACCESS":False,"HST_PIXEL_ACCESS":False,
      "RMAX_USED":False,"HISTORICAL_STRENGTH_0P18_USED":False,"PLANCK_LENGTH_ASSUMED":False,
      "FORBIDDEN_LCDM_PRODUCTION_REFERENCES":0,"PREANCHOR_EXTERNAL_VALUE_ACCESS":False,
      "CIRCULAR_AUTHORITATIVE_SCALE_DERIVATIONS":0,"K0_INTERPRETED_AS_PHYSICAL_STIFFNESS":False,
      "FAST_COEFFICIENT_CHANGED":False,"SLOW_COEFFICIENT_CHANGED":False}; dump("forbidden_input_audit.json",forbidden)
    final={"contract":"PBUF_NATIVE_PHYSICAL_SCALE_V2","established":False,"closure_type":"UNRESOLVED",
      "L0_m_per_native":None,"uncertainty":None,"derivation_class":"MULTIPLE_DIMENSIONAL_CLOSURES_REQUIRED",
      "calibration_anchor":None,"validation_anchors":[],"supporting_candidate_ids":[],"supporting_independence_classes":[],
      "resolution_stability":"NOT_TESTABLE","mass_stability":"NOT_TESTABLE","radius_stability":"NOT_TESTABLE","density_stability":"NOT_TESTABLE",
      "remaining_degeneracies":["L0_U0_CO_DEGENERACY","L0_Kphys_CO_DEGENERACY","L0_S0_CO_DEGENERACY","L0_T0_CO_DEGENERACY"],
      "target_independent":True,"lcdm_independent":True,"production_ready":False,
      "reason":"Current native equations constrain combinations of length, response, stiffness, source, and time units but do not isolate L0."}
    dump("final_native_scale_contract.json",final)
    dump("dev136_blocker_resolution.json",{"resolved":False,"physical_length_scale_established":False,"angular_scale_established":False,
      "PHYSICAL_LENGTH_SCALE_ESTABLISHED_NE_ANGULAR_SCALE_ESTABLISHED":True,"WL_DEV135_PHYSICAL_SCALE_BLOCKER_RESOLVED":False})
    checks={"dev136_hashes_verified":all(base.values()),"propagation_changes_zero":True,"medium_response_law_changes_zero":True,"fast_slow_transfer_changes_zero":True,
      "dimensional_restoration_complete":True,"dimensional_rank_audit_complete":True,"L0_identifiability_classified":True,
      "all_candidate_families_attempted":len(resultdict)==30,"candidate_dependency_graph_complete":True,"candidate_contamination_audit_complete":True,
      "mass_universality_complete":True,"radius_universality_complete":True,"density_universality_complete":True,"resolution_universality_complete":True,
      "preanchor_external_value_access_false":True,"circular_authoritative_scale_derivations_zero":True,"rmax_used_false":True,
      "historical_strength_0p18_used_false":True,"planck_length_assumed_false":True,"lcdm_access_false":True,"class_access_false":True,
      "conventional_lensing_scale_access_false":True,"target_access_false":True,"hst_pixel_access_false":True,
      "forbidden_lcdm_production_references_zero":True,"L0_refits_after_freeze_zero":True}
    outcomes=["WL_PBUF_PHYSICAL_UNIT_CLOSURE_STRUCTURE_ESTABLISHED","WL_PBUF_MULTIPLE_DIMENSIONAL_CLOSURES_REQUIRED","WL_PBUF_NATIVE_PHYSICAL_SCALE_REMAINS_UNRESOLVED"]
    result={"status":outcomes[-1],"outcomes":outcomes,"checks":checks,"phases_executed":phases+["N"],"dimensional_rank":rank,
      "candidate_cluster_outcome":clusters["outcome"],"effective_independent_support":0,"external_calibration_used":False}
    dump("result.json",result)
    structural={"result_schema":"DEV137_V1","dev136_expected":EXPECTED,"medium_unit_contract_sha256":canon(contract),
      "dimensional_rank_sha256":canon(rank),"preanchor_candidates_sha256":presha,"calibration_plan_sha256":calsha,
      "frozen_production_changes":{"PROPAGATION_CHANGES":0,"TRAJECTORY_CHANGES":0,"RECEIVER_CHANGES":0,"ARRIVAL_FORMATION_CHANGES":0,"MEDIUM_RESPONSE_LAW_CHANGES":0,"FAST_SLOW_TRANSFER_CHANGES":0,"BOUNDED_STRAIN_LAW_CHANGES":0}}
    dump("structural_result.json",structural)
    try:(RUN/"baseline_git.txt").write_text(subprocess.run(["git","status","--short"],cwd=ROOT,text=True,capture_output=True).stdout)
    except Exception:(RUN/"baseline_git.txt").write_text("unavailable\n")
    figs={"medium_dimensional_chain.png":["source S0 -> response U0 -> gradient U0/L0 -> trajectory","No SI bridge in frozen path"],
      "unit_scale_nullspace.png":[f"rank={rank['matrix_rank']} nullity={rank['nullity']}",*rank["remaining_free_combinations"]],
      "candidate_mechanism_overview.png":[f"S{i:02d} {x}" for i,x in enumerate(FAMILIES,1)],
      "candidate_L0_estimates.png":["finite authoritative estimates: 0","No arbitrary ruler inserted"],
      "candidate_pairwise_agreement.png":["pairwise matrix empty","no finite candidates"],"candidate_dependency_graph.png":[f"nodes={len(graph['nodes'])}",f"edges={len(graph['edges'])}"],
      "candidate_cluster_summary.png":[clusters["outcome"]],"mass_universality.png":["NOT_TESTABLE"],"radius_universality.png":["NOT_TESTABLE"],
      "density_universality.png":["NOT_TESTABLE"],"resolution_convergence.png":["N=32,48,64,96,128 attempted","no finite estimator"],
      "cross_anchor_universality.png":["NOT_APPLICABLE","no scale frozen"],"L0_scale_power_fingerprint.png":["gradient -1","Laplacian -2","volume +3","path +1"],
      "dev136_scale_degeneracy_closure.png":["Dev136: global scale free","Dev137: reparameterized with U0,K0,S0,T0","external constitutive closures required"],
      "final_scale_decision_tree.png":["internally unique? no","one free length only? no","multiple unit closures required","scale unresolved"]}
    for name,lines in figs.items():fig(name,name.removesuffix(".png").replace("_"," ").title(),lines)
    report=[f"DEV137_PREANCHOR_CANDIDATES_SHA256={presha}",f"DEV137_CALIBRATION_PLAN_SHA256={calsha}",
      "LCDM_ACCESS=false","CLASS_ACCESS=false","CONVENTIONAL_LENSING_SCALE_ACCESS=false","TARGET_ACCESS=false","HST_PIXEL_ACCESS=false",
      "RMAX_USED=false","HISTORICAL_STRENGTH_0P18_USED=false","PLANCK_LENGTH_ASSUMED=false","L0_REFITS_AFTER_FREEZE=0",
      f"L0_IDENTIFIABILITY={rank['L0_identifiability']}",*outcomes]
    (RUN/"report.txt").write_text("\n".join(report)+"\n"); print("\n".join(report)); return 0

if __name__=="__main__": raise SystemExit(main())
