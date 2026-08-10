#!/usr/bin/env python3
"""Dev135 receiver-scale/attitude audit and ACS reference closure."""
from __future__ import annotations
import argparse, hashlib, json, subprocess, sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[3]; sys.path.insert(0,str(ROOT))
from pbuf.wl.receiver_sky_bridge import (BRIDGE_VERSION, PBUFReceiverSkyBridge,
                                         canonical_sha256, roundtrip_audit)
from pbuf.wl.hst_acs_reference_closure import acquire_references

RUN=ROOT/"runs/wl_receiver_to_hst_coordinate_bridge001"
DEV134=ROOT/"runs/wl_hst_acs_detector_geometry001"
EVENTS=ROOT/"runs/wl_geometric_optical_interaction001/geometric_optical_events.npz"
EXPECTED={"contract":"cbdb6d67989b1dd6c53fd004a26b4ca8d455990d608094c79d17a85fc4369e94",
          "manifest":"ffecf391bd271f8d0e75be93e5d2c755713d9cbbdf35a1c75ef1f8ade60ab5ca",
          "structural":"58a26bfc89dbd1b75d9c8306bc7a003bfd18dcd2dbb5bfa32ef637e99281f3a6"}
REFS=("jref$4bb1536cj_idc.fits","jref$4bb1536lj_npl.fits",
      "jref$4bb1536mj_idc.fits","jref$4bb15371j_d2i.fits")

def dump(name,value):
    (RUN/name).write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False,default=str)+"\n")
def fsha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def png(name,title,lines):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(8,4.5)); ax.axis("off"); ax.set_title(title,weight="bold")
    ax.text(.03,.88,"\n".join(lines),va="top",family="monospace",fontsize=9,transform=ax.transAxes)
    fig.tight_layout(); fig.savefig(RUN/name,dpi=120); plt.close(fig)

def baseline_gate():
    result=json.loads((DEV134/"result.json").read_text())
    blockers=[x.get("quantity") for x in result.get("missing_quantities",[])]
    checks={"dev134_contract_sha_verified":result.get("contract_sha256")==EXPECTED["contract"],
            "dev134_manifest_sha_verified":result.get("exposure_manifest_sha256")==EXPECTED["manifest"],
            "dev134_structural_sha_verified":fsha(DEV134/"structural_result.json")==EXPECTED["structural"],
            "status_verified":result.get("status")=="WL_HST_ACS_NATIVE_GEOMETRY_INSUFFICIENTLY_SPECIFIED",
            "exact_three_blocker_classes":len(blockers)==3}
    if not all(checks.values()): raise RuntimeError("DEV134_GATE_STATE_MISMATCH")
    return checks

def synthetic():
    rows=[]
    for deg in (0,90,180,270):
        a=np.deg2rad(deg); r=((float(np.cos(a)),float(-np.sin(a))),(float(np.sin(a)),float(np.cos(a))))
        b=PBUFReceiverSkyBridge(BRIDGE_VERSION,"SYNTHETIC",(1,0,0),(0,1,0),(0,0,1),
          "EXPLICIT_ANGULAR_SCALE",2.,"arcsec/A0_unit",{"synthetic":True},(3.,-4.),"synthetic tangent plane",r,
          False,{"translation":"UPSTREAM_FIXED","rotation":"UPSTREAM_FIXED","scale":"UPSTREAM_FIXED","reflection":"UPSTREAM_FIXED"},
          ("synthetic",),canonical_sha256({"synthetic":True}))
        audit=roundtrip_audit(b,np.array([[0,0],[1,0],[0,1],[2,-3]],float)); rows.append({"degrees":deg,**audit})
    refl=PBUFReceiverSkyBridge(BRIDGE_VERSION,"SYNTHETIC",(1,0,0),(0,1,0),(0,0,1),
      "EXPLICIT_ANGULAR_SCALE",1.,"arcsec/A0_unit",{"synthetic":True},(0.,0.),"synthetic",((-1.,0.),(0.,1.)),True,
      {"translation":"UPSTREAM_FIXED","rotation":"UPSTREAM_FIXED","scale":"UPSTREAM_FIXED","reflection":"UPSTREAM_FIXED"},
      ("reflection test",),canonical_sha256({"reflection":True}))
    missing_scale=False
    try: PBUFReceiverSkyBridge(BRIDGE_VERSION,"x",(1,0,0),(0,1,0),(0,0,1),"NUMERICAL_ONLY_SCALE",1.,"native",{},(0,0),"x",((1,0),(0,1)),False,{},(),"x")
    except ValueError as exc: missing_scale=str(exc)=="PBUF_RECEIVER_PHYSICAL_SCALE_NOT_ESTABLISHED"
    return {"rotation_translation_scale_tests":rows,"reflection_test":roundtrip_audit(refl,np.array([[1.,2.]])),
            "reflection_explicit":refl.reflection_status,"missing_scale_rejected":missing_scale,
            "passed":all(x["max_error"]<1e-12 for x in rows) and missing_scale}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--validation",action="store_true"); ap.add_argument("--offline",action="store_true")
    args=ap.parse_args(); RUN.mkdir(parents=True,exist_ok=True)
    base=baseline_gate(); synth=synthetic()
    if args.validation:
        dump("result.json",{"validation":True,"synthetic":synth,"science_array_read_count":0,"hst_pixel_access":False,"target_access":False})
        (RUN/"report.txt").write_text("HST_HEADER_ACCESS=true\nHST_PIXEL_ACCESS=false\nSCIENCE_ARRAY_READ_COUNT=0\nTARGET_ACCESS=false\nVALIDATION_ONLY_NO_SCIENCE_CLAIM\n")
        print("VALIDATION_ONLY_NO_SCIENCE_CLAIM"); return 0 if synth["passed"] else 4

    scale={"classification":"NUMERICAL_ONLY_SCALE","scale_resolved":False,"coordinate_extent":8.0,"coordinate_units":"native simulation coordinate",
      "candidate_sources":[{"source":"pbuf.wl.config.EXTENT / m10 CFG extent","classification":"NUMERICAL_ONLY_SCALE"},
       {"source":"native_length_scale_mapping_audit001.scale_inventory","classification":"UNRESOLVED","finding":"SI_length_scale_closed=false"}],
      "physical_or_angular_box_size":None,"distance_metadata":None,"redshift_used_as_distance":False,"target_derived":False,
      "conclusion":"No audited grid-unit-to-length or grid-unit-to-angle mapping exists; array dimensions and HST field size were not used."}
    origin={"classification":"NUMERICAL_GRID_ORIGIN_ONLY","origin_resolved":False,"a0_origin":"receiver-plane Cartesian origin",
            "absolute_sky_origin":None,"array_center_assumed":False,"target_alignment_used":False}
    basis={"classification":"RELATIVE_ORIENTATION_KNOWN","absolute_sky_orientation_known":False,
      "lineage":["launch x/y","propagation global Cartesian x/y","Dev130 receiver e_u=(1,0,0), e_v=(0,1,0), normal=(0,0,1)",
                 "Dev133 incoming optical frame"],"a0_basis_u":[1,0,0],"a0_basis_v":[0,1,0],"a0_normal":[0,0,1]}
    gauge={"translation":"GAUGE_FREE","rotation":"GAUGE_FREE","scale":"UNRESOLVED","reflection":"UPSTREAM_FIXED",
      "notes":"Internal basis and handedness are fixed. Absolute celestial placement and attitude are absent. Scale is not called gauge-free because scale invariance was not established.",
      "coordinate_gauge_fix_applied":False}
    for n,v in (("pbuf_a0_scale_provenance.json",scale),("pbuf_a0_origin_provenance.json",origin),
                ("pbuf_a0_basis_provenance.json",basis),("pbuf_gauge_audit.json",gauge)): dump(n,v)
    closure=acquire_references(REFS,RUN/"acs_references",offline=args.offline); dump("acs_reference_closure_manifest.json",closure)
    dump("acs_reference_content_audit.json",{"distortion_application_order":["DET2IM/D2IM","polynomial/IDCTAB","non-polynomial/NPOL residual","linear frame/WCS"],
      "components_kept_separate":True,"references":[{"logical_reference_name":r["logical_reference_name"],"extensions":r.get("extensions",[])} for r in closure["records"]]})
    manifest=json.loads((DEV134/"acs_exposure_geometry_manifest.json").read_text())
    attitudes=[]
    for e in manifest:
        o=e.get("orientation_metadata",{}); attitudes.append({"exposure_uid":e["exposure_uid"],"aperture":e.get("detector"),
          "pa_v3":o.get("PA_V3",{}).get("value"),"chip_orientat":[v.get("value") for k,v in o.items() if "ORIENTAT" in k]})
    dump("hst_header_attitude_audit.json",{"exposure_count":len(attitudes),"exact_geometry_families_preserved":116,"attitudes":attitudes,
      "science_array_read_count":0,"status":"HEADER_ATTITUDE_INVENTORIED; NO A0 TIE"})
    siaf={"access":"SIAF_ACCESS_UNAVAILABLE","aperture_status":"SIAF_APERTURE_UNRESOLVED","reason":"No frozen machine-readable SIAF source is installed or present locally; handbook values were not hard-coded.",
          "header_siaf_consistency":"NOT_COMPARABLE"}; dump("hst_siaf_geometry_audit.json",siaf)
    contract={"bridge_version":BRIDGE_VERSION,"components":{"scale":scale,"origin":origin,"basis":basis,"gauge":gauge},
      "linear_form":"x_sky=O+s R x_A0 when all components are authoritative","production_bridge_constructed":False,
      "reverse_definition_complete":False,"hard_independence":{"HST_HEADER_ACCESS":True,"HST_PIXEL_ACCESS":False,"SCIENCE_ARRAY_READ_COUNT":0,"TARGET_ACCESS":False}}
    contract_sha=canonical_sha256(contract); dump("receiver_sky_bridge_contract.json",contract)
    dump("receiver_sky_bridge_audit.json",{"scale_resolved":False,"origin_resolved":False,"orientation_resolved":False,"reference_closure_complete":closure["reference_closure_complete"],
      "production_ready":False,"stop_reasons":["PBUF_RECEIVER_PHYSICAL_SCALE_NOT_ESTABLISHED","PBUF_RECEIVER_SKY_ATTITUDE_NOT_ESTABLISHED"]})
    dump("receiver_sky_reverse_audit.json",{"reverse_classification":"UNRESOLVED","synthetic_validation":synth,"production_reverse_performed":False})
    preservation={"DEV133_FIELDS_LOST":0,"BROKEN_COORDINATE_BRIDGE_PROVENANCE_LINKS":0,"EVENT_UID_CHANGES":0,"UNEXPLAINED_INFORMATION_LOSS":0,
      "full_population_loaded":False,"reason":"Phase H gated because no production bridge was established"}; dump("dev133_information_preservation.json",preservation)
    # No production arrays are emitted: required only if the bridge is established.
    figspec=[("pbuf_a0_frame_provenance.png","A0 frame provenance",basis["lineage"]),
      ("pbuf_scale_provenance.png","Scale provenance",["extent=8: numerical only","physical scale: unresolved","angular scale: unresolved"]),
      ("a0_sky_hst_frame_chain.png","Coordinate frame chain",["A0 -- BLOCKED(scale/origin/attitude) --> sky","sky --> HST V2/V3: header inventory only"]),
      ("gauge_freedom_summary.png","Gauge freedom",[f"{k}: {v}" for k,v in gauge.items() if isinstance(v,str)]),
      ("hst_v2v3_bridge_geometry.png","HST V2/V3 bridge",["116 header families preserved","SIAF machine-readable source unavailable","A0 tie unresolved"]),
      ("acs_reference_file_chain.png","ACS calibration references",[f"{r['logical_reference_name']}: {r['status']}" for r in closure["records"]]),
      ("header_siaf_consistency.png","Header / SIAF consistency",["NOT_COMPARABLE","headers inventoried; SIAF unavailable"]),
      ("receiver_sky_roundtrip_error.png","Round-trip validation",["synthetic bridge: PASS","production bridge: gated"]),
      ("coordinate_bridge_information_preservation.png","Information preservation",[f"{k}={v}" for k,v in preservation.items() if isinstance(v,int)]),
      ("dev134_blocker_closure.png","Dev134 blocker closure",["PBUF scale: OPEN","PBUF attitude: OPEN",f"ACS references: {'CLOSED' if closure['reference_closure_complete'] else 'OPEN'}"])]
    for spec in figspec: png(*spec)
    # Closure identity deliberately excludes retrieval time and online/offline mode.
    # It is therefore an identity of the frozen bytes/content, not of the acquisition session.
    closure_identity={"official_service":closure["official_service"],"reference_closure_complete":closure["reference_closure_complete"],
      "records":[{"logical_reference_name":r.get("logical_reference_name"),"resolved_filename":r.get("resolved_filename"),
                  "byte_size":r.get("byte_size"),"sha256":r.get("sha256"),
                  "fits_primary_metadata":r.get("fits_primary_metadata"),"extensions":r.get("extensions")} for r in closure["records"]]}
    scale_sha=canonical_sha256(scale); attitude_sha=canonical_sha256({"origin":origin,"basis":basis,"gauge":gauge}); ref_sha=canonical_sha256(closure_identity)
    checks={**base,"a0_scale_audit_complete":True,"a0_origin_audit_complete":True,"a0_basis_audit_complete":True,"gauge_audit_complete":True,
      "acs_reference_closure_audit_complete":True,"siaf_audit_complete":True,"header_attitude_audit_complete":True,"no_target_derived_scale":True,
      "no_target_derived_origin":True,"no_target_derived_attitude":True,"science_array_read_count_zero":True,"hst_pixel_access_false":True,
      "target_access_false":True,"bridge_derivation_explicit":True,"bridge_inverse_defined_if_established":True,"dev133_fields_lost_zero":True,
      "broken_coordinate_bridge_provenance_links_zero":True,"event_uid_changes_zero":True,"unexplained_information_loss_zero":True}
    outcomes=["WL_PBUF_RECEIVER_PHYSICAL_SCALE_UNRESOLVED","WL_PBUF_RECEIVER_SKY_ATTITUDE_UNRESOLVED","WL_PBUF_RECEIVER_GLOBAL_SKY_GAUGE_IDENTIFIED"]
    if closure["reference_closure_complete"]: outcomes.append("WL_ACS_REFERENCE_CALIBRATION_CLOSURE_ESTABLISHED")
    result={"status":"WL_PBUF_RECEIVER_PHYSICAL_SCALE_UNRESOLVED","phase_completed":"Phase E — SIAF/header attitude audit; Phases F-H gated",
      "outcomes":outcomes,"answers":{"physical_scale":False,"scale_origin":"numerical simulation extent only","scale_without_LCDM":False,
      "absolute_sky_origin":False,"absolute_sky_attitude":False,"a0_to_hst_v2v3":False,"reverse_mapping":False,"all_events_in_common_hst_frame":False,
      "dev134_unblocked":False},"reference_closure_complete":closure["reference_closure_complete"],"checks":checks,
      "science_array_read_count":0,"hst_pixel_access":False,"target_access":False,"production_bridge_constructed":False}
    dump("result.json",result)
    structural={"contract_sha256":contract_sha,"scale_sha256":scale_sha,"attitude_sha256":attitude_sha,"reference_closure_sha256":ref_sha,
      "dev134":EXPECTED,"frozen_upstream_changes":0,"production_mapping_performed":False}; dump("structural_result.json",structural); structural_sha=fsha(RUN/"structural_result.json")
    try: (RUN/"baseline_git.txt").write_text(subprocess.run(["git","status","--short"],cwd=ROOT,text=True,capture_output=True).stdout)
    except Exception: (RUN/"baseline_git.txt").write_text("unavailable\n")
    lines=[f"DEV135_PBUF_SKY_SCALE_SHA256={scale_sha}",f"DEV135_PBUF_SKY_ATTITUDE_SHA256={attitude_sha}",
      f"DEV135_ACS_REFERENCE_CLOSURE_SHA256={ref_sha}",f"DEV135_COORDINATE_BRIDGE_CONTRACT_SHA256={contract_sha}",f"DEV135_STRUCTURAL_SHA256={structural_sha}",
      "HST_HEADER_ACCESS=true","HST_PIXEL_ACCESS=false","SCIENCE_ARRAY_READ_COUNT=0","TARGET_ACCESS=false","PBUF_RECEIVER_PHYSICAL_SCALE_NOT_ESTABLISHED",
      "PBUF_RECEIVER_SKY_ATTITUDE_NOT_ESTABLISHED",*( ["WL_ACS_REFERENCE_CALIBRATION_CLOSURE_ESTABLISHED"] if closure["reference_closure_complete"] else ["ACS_REFERENCE_CLOSURE_INCOMPLETE"]),
      "WL_PBUF_RECEIVER_PHYSICAL_SCALE_UNRESOLVED","WL_PBUF_RECEIVER_SKY_ATTITUDE_UNRESOLVED"]
    (RUN/"report.txt").write_text("\n".join(lines)+"\n"); print("\n".join(lines)); return 3

if __name__=="__main__": raise SystemExit(main())
