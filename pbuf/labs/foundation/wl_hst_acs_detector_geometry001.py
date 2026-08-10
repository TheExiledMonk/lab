#!/usr/bin/env python3
"""Dev134 ACS/WFC native-geometry gate and synthetic validation lane."""
from __future__ import annotations
import argparse, hashlib, json, subprocess, sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from pbuf.data.hst_acs_calibration_audit import canonical_sha256
from pbuf.wl.hst_acs_header_geometry import exposure_manifest_sha256, inventory
from pbuf.wl.hst_acs_geometry import (CONTRACT_VERSION, AffineDetectorTransform, RectangularChip,
    basis_audit, classify_detector_points)

RUN=ROOT/"runs/wl_hst_acs_detector_geometry001"
DATA=ROOT/"PBUF_raw_benchmark/WLRAW-001_Abell2744"
UP=ROOT/"runs/wl_geometric_optical_interaction001"
EXPECTED={"optical_system_sha256":"0a4afbdc18259113b69e8e2a82cfcfb03a65c6d4dec488e9efa18c7c60ec4cff",
          "contract_sha256":"5ab9f1e490154ebb90bc3a5431e2c406331274508d96a8983577044ee1006ddd",
          "structural_sha256":"580249b502680478d80a16bafdd8035b8c624ee4f8512c506ad0fdd7389dfd98"}

def dump(name, value):
    (RUN/name).write_text(json.dumps(value,indent=2,sort_keys=True,default=lambda x:x.item() if isinstance(x,np.generic) else str(x))+"\n")
def file_sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def synthetic_audit():
    points=np.array([[0.,0.],[1.,0.],[0.,1.],[-1.,0.],[0.,-1.]])
    rows=[]
    for deg in (0,90,180,270):
        a=np.deg2rad(deg);m=np.array([[np.cos(a),-np.sin(a)],[np.sin(a),np.cos(a)]])
        t=AffineDetectorTransform(m,np.array([10.,20.]),f"ROT{deg}");out=t.forward(points);back=t.reverse(out)
        rows.append({"orientation_degrees":deg,"max_roundtrip_error":float(np.max(np.abs(back-points))),"determinant":float(np.linalg.det(m))})
    reflection=AffineDetectorTransform(np.diag([-1.,1.]),np.zeros(2),"REFLECTION")
    chips=(RectangularChip("CHIP_1",0,10,0,4),RectangularChip("CHIP_2",0,10,6,10))
    probe=np.array([[5,2],[5,8],[5,5],[-1,2],[0,2],[np.nan,0]],float)
    status,ids=classify_detector_points(probe,chips)
    singular=False
    try: AffineDetectorTransform(np.array([[1.,0.],[0.,0.]]),np.zeros(2)).reverse([[0,0]])
    except ValueError: singular=True
    return {"rotation_tests":rows,"reflection_determinant":float(np.linalg.det(reflection.matrix)),
      "two_chip_status":status.tolist(),"two_chip_ids":ids.tolist(),"singular_rejected":singular,
      "basis":basis_audit([1,0,0],[0,1,0],[0,0,1]),"passed":all(r["max_roundtrip_error"]<1e-14 for r in rows) and singular and
      status.tolist()==["ACTIVE_CHIP","ACTIVE_CHIP","INTER_CHIP_GAP","OUTSIDE_DETECTOR","BOUNDARY","INVALID"]}

def contract():
    return {"version":CONTRACT_VERSION,"frames":{"A0":"PBUF incoming optical frame","A1":"telescope/instrument reference frame",
      "A2":"ACS/WFC ideal focal/detector reference frame","A3":"chip-native physical frame","A4":"calibrated FLT/FLC array frame","A5":"diagnostic sky/WCS only"},
      "source_hierarchy":["LOCAL_FITS_HEADER","LOCAL_REFERENCE_FILE","OFFICIAL_INSTRUMENT_DOCUMENT","DERIVED_FROM_FROZEN_SOURCE"],
      "boundary_rule":"lower inclusive; upper exclusive; exact edges labeled BOUNDARY","coordinates":"continuous float64; no pixel index or accumulation",
      "fits_origin":"FITS 1-based metadata retained separately from continuous physical and numpy coordinates",
      "losses":["DETECTOR_GEOMETRIC_GAP_LOSS","OUTSIDE_ACTIVE_DETECTOR","DISTORTION_DOMAIN_FAILURE"],
      "reverse_lanes":["FULL_METADATA_REVERSE","OUTPUT_GEOMETRY_ONLY_REVERSE"],"science_array_access":False,"image_generation":False}

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--validation",action="store_true");args=ap.parse_args();RUN.mkdir(parents=True,exist_ok=True)
    prior=json.loads((UP/"result.json").read_text())
    if any(prior.get(k)!=v for k,v in EXPECTED.items()) or file_sha(UP/"structural_result.json")!=EXPECTED["structural_sha256"]:
        print("DEV133_STRUCTURAL_BASELINE_MISMATCH");return 2
    required={"WL_GEOMETRIC_OPTICAL_INTERACTION_ESTABLISHED","WL_GEOMETRIC_OPTICAL_EVENT_PRESERVATION_ESTABLISHED","WL_GEOMETRIC_OPTICAL_REVERSE_TRANSPORT_ESTABLISHED","WL_COMBINED_PBUF_OPTICAL_TRANSPORT_STRUCTURE_ESTABLISHED","WL_OPTICAL_APERTURE_INFORMATION_LOSS_BOUNDARY_ESTABLISHED","WL_BIDIRECTIONAL_GEOMETRIC_OPTICAL_ARCHITECTURE_ESTABLISHED"}
    if not required.issubset(prior.get("outcomes",[])) or prior.get("DEV132_FIELDS_LOST") or prior.get("BROKEN_OPTICAL_PROVENANCE_LINKS") or prior.get("EVENT_UID_CHANGES"): return 2
    c=contract();dump("acs_wfc_geometry_contract.json",c);contract_sha=canonical_sha256(c)
    synth=synthetic_audit();dump("acs_frame_audit.json",synth)
    if args.validation:
        dump("result.json",{"validation":True,"synthetic_tests_passed":synth["passed"],"science_array_read_count":0,"hst_header_access":True,"hst_pixel_access":False})
        (RUN/"report.txt").write_text(f"DEV134_ACS_GEOMETRY_CONTRACT_SHA256={contract_sha}\nHST_HEADER_ACCESS=true\nHST_PIXEL_ACCESS=false\nSCIENCE_ARRAY_READ_COUNT=0\nVALIDATION_ONLY_NO_SCIENCE_CLAIM\n")
        print("VALIDATION_ONLY_NO_SCIENCE_CLAIM");return 0 if synth["passed"] else 4
    exposures=inventory(DATA);manifest=[e.manifest() for e in exposures];dump("acs_exposure_geometry_manifest.json",manifest)
    manifest_sha=exposure_manifest_sha256(exposures)
    refs=sorted({str(v.value) for e in exposures for v in e.reference_file_metadata.values()})
    missing_refs=[r for r in refs if "$" in r and not any(DATA.rglob(Path(r.split("$",1)[1]).name))]
    missing=[{"quantity":"A0 PBUF receiver-coordinate angular/physical scale","required_for":"A0->A1/A2 transform","reason":"no units or authoritative scale in Dev132/133 state"},
      {"quantity":"A0 origin/attitude tie to an exposure sky or V2/V3 reference","required_for":"exposure pointing transform","reason":"event state contains receiver coordinates but no RA/Dec or instrument-reference tie"}]
    if missing_refs: missing.append({"quantity":"referenced distortion/calibration files","identifiers":missing_refs,"required_for":"complete ideal-to-physical distortion audit"})
    source={"dataset":str(DATA.resolve()),"exposure_count":len(exposures),"HST_HEADER_ACCESS":True,"HST_PIXEL_ACCESS":False,
      "SCIENCE_ARRAY_READ_COUNT":0,"reference_identifiers":refs,"unresolved_quantities":missing}
    dump("acs_geometry_source_manifest.json",source)
    # Exact-family grouping uses header geometry only.
    groups={}
    for e in exposures:
        key=canonical_sha256({"chips":e.manifest()["chips"],"orientation":e.manifest()["orientation_metadata"],"references":e.manifest()["reference_file_metadata"]})
        groups.setdefault(key,[]).append(e.exposure_uid)
    dump("acs_geometry_families.json",{"number_of_unique_geometry_families":len(groups),"families":[{"geometry_sha256":k,"members":v} for k,v in sorted(groups.items())]})
    dump("acs_header_triplet_audit.json",{"exposure_count":len(exposures),"status":"HEADER_INVENTORY_COMPLETE; semantic RAW/FLT/FLC transform comparison pending missing bridge resolution","science_array_read_count":0})
    dump("acs_distortion_audit.json",{"status":"DISTORTION_MODEL_PARTIALLY_AVAILABLE","header_sip_available":True,"referenced_files_locally_available":not missing_refs,"unresolved_reference_identifiers":missing_refs})
    structural={"dev133":EXPECTED,"dataset_identity":str(DATA.resolve()),"contract_sha256":contract_sha,"exposure_manifest_sha256":manifest_sha,"source_manifest":source,
      "mapping_frozen":False,"stop_reason":"ACS_NATIVE_GEOMETRY_INCOMPLETE_FROM_LOCAL_SOURCES"};dump("structural_result.json",structural);structural_sha=file_sha(RUN/"structural_result.json")
    result={"validation":False,"phase_completed":"Phase A — header/reference inventory","exposure_count":len(exposures),"logical_mapping_count":285156*len(exposures),
      "production_mapping_performed":False,"status":"WL_HST_ACS_NATIVE_GEOMETRY_INSUFFICIENTLY_SPECIFIED","missing_quantities":missing,"contract_sha256":contract_sha,
      "exposure_manifest_sha256":manifest_sha,"structural_sha256":structural_sha,"science_array_read_count":0,"hst_header_access":True,"hst_pixel_access":False}
    dump("result.json",result)
    try:(RUN/"baseline_git.txt").write_text(subprocess.run(["git","status","--short"],cwd=ROOT,text=True,capture_output=True).stdout)
    except Exception:(RUN/"baseline_git.txt").write_text("unavailable\n")
    lines=[f"DEV134_ACS_GEOMETRY_CONTRACT_SHA256={contract_sha}",f"DEV134_EXPOSURE_GEOMETRY_MANIFEST_SHA256={manifest_sha}",f"DEV134_STRUCTURAL_SHA256={structural_sha}","HST_HEADER_ACCESS=true","HST_PIXEL_ACCESS=false","SCIENCE_ARRAY_READ_COUNT=0","ACS_NATIVE_GEOMETRY_INCOMPLETE_FROM_LOCAL_SOURCES","WL_HST_ACS_NATIVE_GEOMETRY_INSUFFICIENTLY_SPECIFIED"]
    (RUN/"report.txt").write_text("\n".join(lines)+"\n");print("\n".join(lines));return 3

if __name__=="__main__":raise SystemExit(main())

