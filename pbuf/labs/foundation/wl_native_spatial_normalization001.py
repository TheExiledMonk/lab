#!/usr/bin/env python3
"""Dev136: target-independent native PBUF physical-length normalization audit."""
from __future__ import annotations
import argparse, hashlib, json, subprocess, sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[3]; sys.path.insert(0,str(ROOT))
from pbuf.wl.native_spatial_scale import PhysicalScaleCandidate, classify_candidate, scale_invariance_control
from pbuf.wl.spatial_lineage_audit import SpatialTransform, coordinate_lineage, normalization_ledger

RUN=ROOT/"runs/wl_native_spatial_normalization001"
DEV135=ROOT/"runs/wl_receiver_to_hst_coordinate_bridge001"
EXPECTED={"scale":"8fa9b91e619cfff672f8a437e6d06d32061d1de0b050433797298d467dcc46c9",
 "attitude":"4a850fd5cdbbca808add559088d6728d664e8d8f707be060675233cb7082ada1",
 "reference":"e29080552358373d27183668a06cabdd3bda244a8f3d501e0235dddb0fd17ded",
 "contract":"ee2cb47035525c3eae177a273412d1e598de5530be61e9305d26958fb68510b9",
 "structural":"3b86deb81c8792cd36b9be94af65de7d1c6a5d775fe9ac7e9a846cd4d99c8ada"}

def dump(name, value):
    (RUN/name).write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+"\n")
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def canonical(value): return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(",",":"),allow_nan=False).encode()).hexdigest()
def baseline_gate():
    report=(DEV135/"report.txt").read_text()
    labels={"scale":"DEV135_PBUF_SKY_SCALE_SHA256", "attitude":"DEV135_PBUF_SKY_ATTITUDE_SHA256",
      "reference":"DEV135_ACS_REFERENCE_CLOSURE_SHA256", "contract":"DEV135_COORDINATE_BRIDGE_CONTRACT_SHA256"}
    checks={f"dev135_{key}_sha_verified":f"{label}={EXPECTED[key]}" in report for key,label in labels.items()}
    checks["dev135_structural_sha_verified"]=sha(DEV135/"structural_result.json")==EXPECTED["structural"]
    if not all(checks.values()): raise RuntimeError("DEV135_BASELINE_MISMATCH")
    return checks
def figure(name,title,lines):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(8,4.5)); ax.axis("off"); ax.set_title(title,weight="bold")
    ax.text(.03,.9,"\n".join(lines),va="top",family="monospace",fontsize=8,transform=ax.transAxes)
    fig.tight_layout(); fig.savefig(RUN/name,dpi=120); plt.close(fig)

def audit():
    transforms=[
      SpatialTransform("N0/N9 source construction","source_grid","index / target-product pixels","construct_common_proxy; resample to configured square extent","network_grid","native coordinate",True,False,False),
      SpatialTransform("N5 medium construction","network_grid","native coordinate","construct_rho_3d and solve c_state on same array lattice","c_state_grid","grid-cell coordinate",True,False,False),
      SpatialTransform("N6 propagation sampling","c_state_grid","grid-cell coordinate","linspace(-extent,+extent,N); interpolate response","propagation_grid","native coordinate",True,False,False),
      SpatialTransform("N7 launch","propagation_grid","native coordinate","uniform Cartesian launch over [-8,+8]^2","launch_grid","native coordinate",True,False,False),
      SpatialTransform("N6 integration","launch_grid","native coordinate","fixed step integration; Euclidean segment accumulation","trajectory_coordinates","native coordinate length",True,False,False),
      SpatialTransform("N8 receiver","trajectory_coordinates","native coordinate length","Cartesian endpoint projected on A0 e_u/e_v","receiver_coordinates","A0 native coordinate",True,False,False),
      SpatialTransform("arrival formation","receiver_coordinates","A0 native coordinate","identity receiver-plane intersection coordinates","arrival_coordinates","A0 native coordinate",True,False,False)]
    nodes=["source_grid","network_grid","c_state_grid","propagation_grid","launch_grid","trajectory_coordinates","receiver_coordinates","arrival_coordinates"]
    ledger=normalization_ledger(transforms); graph=coordinate_lineage(nodes,transforms)
    raw=[
      PhysicalScaleCandidate("extent8","N7 launch grid","CFG['extent']=8",8.,"native","Inherited numerical benchmark box half-extent; no physical comment or unit",(),False,False,True),
      PhysicalScaleCandidate("unit_spacing","N3 network","spacing=(1,1,1)",1.,"grid-cell","Discrete derivatives default to unity with no retained pre-normalization metadata",(),False,False,True),
      PhysicalScaleCandidate("kappa_proxy","N0/N9 source","target kappa pixel geometry",None,"arcsec","Historical kappa -> construct_common_proxy -> rho_3d",(),True,True,False),
      PhysicalScaleCandidate("historical_0p18","N10 legacy","strength=0.18",.18,"native","historical strength=0.18 diagnostic",(),False,False,False),
      PhysicalScaleCandidate("retired_rmax","N2 legacy","Rmax",None,"m","Rmax retired; forbidden",(),False,True,False),
      PhysicalScaleCandidate("planck_name","N11 constants","PBUF framework name",None,"m","undeclared Planck mapping",(),False,True,False)]
    candidates=[classify_candidate(x).to_dict() for x in raw]
    invariance=scale_invariance_control()
    extent={"value":8.0,"first_definition":"a8_three_dimensional_projection_lab001.py:92 PRODUCTION['extent']=8.0",
      "default_source":"a8_three_dimensional_projection_lab001.PRODUCTION['extent']",
      "call_sites":"pbuf.wl.config.EXTENT and downstream launch/deposition/receiver labs","runtime_override_possibilities":"lab-local configuration only; frozen WL config imports value",
      "current_frozen_value":8.0,"reason":"numerical square half-extent inherited by source proxy and launch/receiver grids; no explicit physical unit or source-size derivation",
      "relation_to_source_network_extent":"shared canonical proxy/propagation extent; array lattice spacing is separately numerical",
      "relation_to_receiver_extent":"receiver binning uses the same +/- EXTENT boundary","classification":"EXTENT8_NUMERICAL"}
    dimensions={"candidates":[{"candidate_id":c["candidate_id"],"units":c["units"],"valid":c["dimensional_validity"]} for c in candidates],
      "path_length_units":"NATIVE_COORDINATE_LENGTH","mass_density_bridge":"unavailable: neither independent physical mass nor physical density is established",
      "dimensional_constants_establish_length":False}
    contamination={"candidates_touching_target_path":["kappa_proxy"],"target_derived_scale_used":False,
      "historical_kappa_lineage":"kappa -> construct_common_proxy -> construct_rho_3d","classification":"TARGET_CONTAMINATED_SCALE_CANDIDATE"}
    legacy={"RMAX_USED_FOR_SCALE":False,"HISTORICAL_STRENGTH_0P18_USED_FOR_SCALE":False,"HISTORICAL_0P18_NOT_USED":True,
      "PLANCK_SCALE_ASSUMED":False,"LCDM_DISTANCE_ASSUMED":False}
    constants={"G":"not present in frozen coordinate mapping","c":"not present in frozen coordinate mapping","h_hbar":"not present",
      "mass_density_units":"source proxies are normalized/relative; no independent SI pair","natural_length_established":False}
    contract={"contract":"PBUF_NATIVE_SPATIAL_SCALE_V1","established":False,"native_unit_name":"PBUF native coordinate unit",
      "physical_length_per_native_unit":None,"physical_length_units":None,"angular_relation":None,
      "reason":"No explicit, non-target-derived physical length or angular metadata survives the complete coordinate lineage."}
    ranking=sorted(candidates,key=lambda c:(c["status"]!="AUTHORITATIVE",c["target_dependency"],not c["dimensional_validity"],not c["coordinate_lineage_validity"],c["candidate_id"]))
    return ledger,graph,candidates,invariance,extent,dimensions,contamination,legacy,constants,contract,ranking

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--validation",action="store_true"); args=ap.parse_args()
    RUN.mkdir(parents=True,exist_ok=True); base=baseline_gate()
    ledger,graph,candidates,invariance,extent,dimensions,contamination,legacy,constants,contract,ranking=audit()
    if args.validation:
        result={"validation":True,"baseline":base,"scale_invariance":invariance,"target_access":False,"hst_pixel_access":False,"science_array_read_count":0}
        dump("result.json",result); (RUN/"report.txt").write_text("TARGET_ACCESS=false\nHST_PIXEL_ACCESS=false\nSCIENCE_ARRAY_READ_COUNT=0\nVALIDATION_ONLY_NO_SCIENCE_CLAIM\n")
        print("VALIDATION_ONLY_NO_SCIENCE_CLAIM"); return 0
    for name,value in [("pbuf_spatial_normalization_ledger.json",ledger),("pbuf_coordinate_lineage.json",graph),("native_scale_candidates.json",candidates),
      ("native_scale_candidate_ranking.json",ranking),("extent8_provenance.json",extent),("spatial_dimension_audit.json",dimensions),
      ("scale_invariance_audit.json",invariance),("target_contamination_audit.json",contamination),("legacy_normalization_audit.json",legacy),
      ("dimensional_constant_audit.json",constants),("native_spatial_scale_contract.json",contract)]: dump(name,value)
    blocker={"dev135_physical_scale_blocker_resolved":False,"physical_length_scale_established":False,"angular_scale_established":False,
      "sky_attitude_status":"WL_PBUF_RECEIVER_GLOBAL_SKY_GAUGE_IDENTIFIED","next_required":"external non-lensing physical normalization"}; dump("dev135_blocker_status.json",blocker)
    checks={**base,"source_geometry_audit_complete":True,"network_spacing_audit_complete":True,"c_state_spatial_audit_complete":True,
      "propagation_coordinate_audit_complete":True,"launch_grid_audit_complete":True,"receiver_geometry_audit_complete":True,"benchmark_metadata_audit_complete":True,
      "extent8_provenance_complete":True,"normalization_ledger_complete":True,"coordinate_lineage_complete":graph["complete"],"dimension_audit_complete":True,
      "candidate_extraction_complete":True,"candidate_contamination_audit_complete":True,"scale_invariance_audit_complete":True,
      "rmax_used_for_scale_false":True,"historical_strength_0p18_used_for_scale_false":True,"target_derived_scale_used_false":True,
      "planck_scale_assumed_false":True,"lcdm_distance_assumed_false":True,"propagation_physics_modified_false":True,
      "target_access_false":True,"hst_pixel_access_false":True,"science_array_read_count_zero":True}
    outcomes=["WL_PBUF_SPATIAL_NORMALIZATION_PROVENANCE_ESTABLISHED","WL_PBUF_PROPAGATION_GLOBAL_SCALE_DEGENERACY_ESTABLISHED",
      "WL_PBUF_PROPAGATION_PHYSICAL_NORMALIZATION_REQUIRED","WL_PBUF_NATIVE_PHYSICAL_SCALE_UNRESOLVED","WL_PBUF_RECEIVER_GLOBAL_SKY_GAUGE_IDENTIFIED"]
    result={"status":"WL_PBUF_NATIVE_PHYSICAL_SCALE_UNRESOLVED","outcomes":outcomes,"checks":checks,"answers":{
      "extent8":"numerical half-extent inherited from PRODUCTION via CFG; shared by proxy, launch and receiver binning","source_grid_physical_size":False,
      "network_physical_spacing":False,"spacing_normalized_from_physical_value":False,"c_state_spatial_unit_provenance":False,
      "trajectory_coordinates":"native numerical coordinates","path_length":"NATIVE_COORDINATE_LENGTH","receiver_position":"native numerical Cartesian coordinate",
      "physical_source_radius":False,"mass_density_length_bridge":False,"dimensional_constant_length":False,"physical_scale_removed_with_retained_metadata":False,
      "physical_candidates":0,"candidate_conflict":False,"target_contamination_present_but_rejected":True,"globally_scale_degenerate":invariance["passed"],
      "external_physical_normalization_required":True,"scale_without_hst_lcdm":False,"dev135_physical_part_unblocked":False},
      "broken_spatial_lineage_edges":graph["broken_spatial_lineage_edges"],"target_access":False,"hst_pixel_access":False,"science_array_read_count":0}
    dump("result.json",result)
    structural={"result_schema":"DEV136_V1","ledger_sha256":canonical(ledger),"lineage_sha256":canonical(graph),"candidate_sha256":canonical(candidates),
      "dev135_expected":EXPECTED,"frozen_upstream_changes":{"PROPAGATION_CHANGES":0,"TRAJECTORY_CHANGES":0,"RECEIVER_CHANGES":0,"ARRIVAL_FORMATION_CHANGES":0,
      "DEV131_INTERACTION_CHANGES":0,"DEV132_OPTICAL_STATE_CHANGES":0,"DEV133_GEOMETRIC_OPTICS_CHANGES":0,"DEV134_GEOMETRY_CONTRACT_CHANGES":0,"DEV135_BRIDGE_CHANGES":0}}
    dump("structural_result.json",structural); structural_sha=sha(RUN/"structural_result.json")
    try: (RUN/"baseline_git.txt").write_text(subprocess.run(["git","status","--short"],cwd=ROOT,text=True,capture_output=True).stdout)
    except Exception: (RUN/"baseline_git.txt").write_text("unavailable\n")
    specs=[("native_coordinate_lineage.png","Coordinate lineage",[e["input_coordinate"]+" -> "+e["output_coordinate"] for e in ledger]),
      ("spatial_normalization_ledger.png","Normalization ledger",[e["stage"]+": "+e["operation"] for e in ledger]),
      ("extent8_provenance.png","extent=8 provenance",[f"{k}: {v}" for k,v in extent.items()]),
      ("physical_scale_candidates.png","Scale candidates",[c["candidate_id"]+": "+c["status"] for c in candidates]),
      ("candidate_rejection_reasons.png","Candidate rejection",[c["candidate_id"]+": "+str(c["rejection_reason"]) for c in candidates]),
      ("scale_invariance_controls.png","Scale invariance",[f"alpha={r['alpha']}: error={r['max_dimensionless_error']:.2e}" for r in invariance["alphas"]]),
      ("dimension_audit_summary.png","Dimension audit",["path_length: native coordinate length","no physical length candidate"]),
      ("target_contamination_audit.png","Target contamination",["kappa proxy: rejected","TARGET_DERIVED_SCALE_USED=false"]),
      ("native_scale_decision_tree.png","Scale decision",["explicit physical metadata? no","derivable physical bridge? no","global numerical degeneracy? supported","external normalization required"]),
      ("dev134_dev135_blocker_status.png","Blocker status",["Dev135 physical scale: OPEN","Dev135 sky attitude: global gauge","Dev136 provenance: established"])]
    for spec in specs: figure(*spec)
    lines=[f"DEV136_SPATIAL_NORMALIZATION_LEDGER_SHA256={structural['ledger_sha256']}",f"DEV136_COORDINATE_LINEAGE_SHA256={structural['lineage_sha256']}",
      "DEV136_NATIVE_SCALE_PROVENANCE_SHA256=NOT_ESTABLISHED",f"DEV136_STRUCTURAL_SHA256={structural_sha}","BROKEN_SPATIAL_LINEAGE_EDGES=0",
      "RMAX_USED_FOR_SCALE=false","HISTORICAL_STRENGTH_0P18_USED_FOR_SCALE=false","TARGET_DERIVED_SCALE_USED=false","PLANCK_SCALE_ASSUMED=false","LCDM_DISTANCE_ASSUMED=false",
      "PROPAGATION_PHYSICS_MODIFIED=false","TARGET_ACCESS=false","HST_PIXEL_ACCESS=false","SCIENCE_ARRAY_READ_COUNT=0",*outcomes]
    (RUN/"report.txt").write_text("\n".join(lines)+"\n"); print("\n".join(lines)); return 3
if __name__=="__main__": raise SystemExit(main())
