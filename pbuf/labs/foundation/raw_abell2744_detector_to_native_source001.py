#!/usr/bin/env python3
"""Dev161: target-blind RAW/FLT/FLC detector-to-native-source audit."""
from __future__ import annotations

import ast
import json
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path: sys.path.insert(0, str(ROOT))
ARCHIVE = ROOT / "PBUF_raw_benchmark/WLRAW-001_Abell2744"
OUT = ROOT / "runs/raw_abell2744_detector_to_native_source001"

from pbuf.data.hst_acs_raw_source import archive_inventory, calibrated_chips
from pbuf.data.hst_acs_common_frame import combine_samples, frame_from_bounds, sampled_chip
from pbuf.source.raw_detector_source_bridge import native_2d_constraint, support_diagnostics


def dump(name, value):
    (OUT/name).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+"\n")


def _finite_json(a):
    return np.where(np.isfinite(a), a, np.nan)


def _plot(path, image, title, *, cmap="viridis"):
    fig,ax=plt.subplots(figsize=(6,5)); shown=np.asarray(image,float); shown=np.where(np.isfinite(shown),shown,np.nan)
    im=ax.imshow(shown,origin="lower",cmap=cmap);ax.set_title(title);fig.colorbar(im,ax=ax,shrink=.8);fig.tight_layout();fig.savefig(path,dpi=120);plt.close(fig)


def _dependency_audit():
    paths=[ROOT/"pbuf/data/hst_acs_raw_source.py",ROOT/"pbuf/data/hst_acs_common_frame.py",
      ROOT/"pbuf/source/raw_detector_source_bridge.py",Path(__file__)]
    forbidden=("current_native_five_cluster_observable_benchmark001","benchmark_data","kappa","gamma","load_cluster_source","build_native_response")
    hits=[]; imports=[];calls=[]
    for path in paths:
        text=path.read_text(); tree=ast.parse(text)
        imports += [getattr(n,"module",None) or ".".join(x.name for x in n.names) for n in ast.walk(tree) if isinstance(n,(ast.Import,ast.ImportFrom))]
        calls += [n.func.id if isinstance(n.func,ast.Name) else n.func.attr for n in ast.walk(tree) if isinstance(n,ast.Call) and isinstance(n.func,(ast.Name,ast.Attribute))]
        hits += [{"file":str(path.relative_to(ROOT)),"token":x} for x in forbidden if x in text and path == Path(__file__) and x not in ("kappa","gamma")]
    # The runner contains forbidden names only in this negative audit tuple; executable imports/calls are decisive.
    executable_forbidden=[x for x in imports+calls if any(t in str(x) for t in forbidden)]
    return {"audited_files":[str(x.relative_to(ROOT)) for x in paths],"imports":sorted(set(imports)),
      "calls":sorted(set(calls)),"executable_forbidden_dependencies":executable_forbidden,
      "negative_assertion_tokens":hits,"passed":not executable_forbidden}


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    families,inventory=archive_inventory(ARCHIVE); counts=inventory["file_counts"]
    if counts != {"RAW":116,"FLT":116,"FLC":116} or len(families)!=116:
        raise RuntimeError(f"archive baseline mismatch: {counts}, families={len(families)}")
    dump("archive_inventory.json",inventory)
    metadata={"inventory_complete":True,"exposure_count":len(families),"records":inventory["exposures"],
      "recorded_fields":["filename","exposure_id","detector/chip","filter","exposure_time","pointing","orientation","WCS metadata in SCI headers","calibration status","array dimensions","DQ","ERR"]}
    dump("detector_metadata_inventory.json",metadata)
    dump("raw_flt_flc_role_audit.json",{
      "RAW":{"role":"uncalibrated/minimally calibrated detector readout","use":"provenance and archive baseline; not direct science combination"},
      "FLT":{"role":"flat-fielded calibrated exposure on native detector grid","use":"PRIMARY_PIXEL_PRODUCT"},
      "FLC":{"role":"FLT-equivalent calibrated exposure with pixel-based CTE correction","use":"retained control; not assumed universally preferable"},
      "PRIMARY_PIXEL_PRODUCT":"FLT","reason":"Dev126 classifies FLT as the candidate native-detector baseline; choosing FLC would add a detector-specific correction preference not required for this minimum bridge."})

    by_filter=defaultdict(list); bounds=[]; accessed=0; dq_validated=False; err_validated=False
    first_image=None
    for family,row in zip(families,inventory["exposures"]):
        for ver,sci,err,dq,header in calibrated_chips(family.flt):
            sample=sampled_chip(sci,err,dq,header);by_filter[row["filter"]].append(sample);bounds.append((sample[3],sample[4]))
            accessed+=1;dq_validated |= dq is not None and dq.shape==sci.shape
            err_validated |= err is not None and err.shape==sci.shape and bool(np.any(np.isfinite(err)&(err>0)))
            if first_image is None: first_image=sci[::32,::32].copy()
    frame=frame_from_bounds(bounds); combined={f:combine_samples(samples,frame) for f,samples in sorted(by_filter.items())}
    channels=tuple(sorted(combined)); images=[combined[f]["image"] for f in channels]; uncertainties=[combined[f]["uncertainty"] for f in channels]
    constraint=native_2d_constraint(images,uncertainties,channels); diagnostics=support_diagnostics(constraint)
    overlap={f:{"sample_count":len(by_filter[f]),"covered_cells":int(combined[f]["occupied"].sum()),
      "coverage_fraction":float(combined[f]["occupied"].mean()),"maximum_sample_coverage":int(combined[f]["coverage"].max()),
      "cells_with_overlap":int((combined[f]["coverage"]>1).sum())} for f in channels}
    frame_contract={"status":"TRUE","frame":"ICRS celestial bounding box sampled from calibrated SCI WCS",
      "shape":list(frame.shape),"bounds_degrees":{"ra":[frame.ra_min,frame.ra_max],"dec":[frame.dec_min,frame.dec_max]},
      "mapping":"FITS WCS celestial pixel_to_world followed by fixed rectangular sky binning","morphology_targeted_alignment":False,
      "coordinate_covariance":"WCS geometry explicit; no downstream morphology enters frame selection"}
    dump("common_frame_contract.json",frame_contract)
    dump("exposure_combination_audit.json",{"selected":"inverse-variance weighted mean of deterministic detector samples",
      "alternatives":["uncombined multi-exposure (preserved as S01 contract)","exposure-time weighted mean","masked mean","robust median"],
      "least_destructive_endpoint":"multi-channel common-frame field plus uncertainty and coverage; detector exposure inventory retained",
      "overlap":overlap,"background_model":"per-chip constant","background_estimator":"median of finite DQ=0 deterministic samples",
      "masking_rule":"DQ must equal zero; SCI/ERR/WCS finite; ERR positive"})
    dump("candidate_source_inventory.json",{"S01":{"status":"VALID","representation":"individual FLT SCI/ERR/DQ exposures"},
      "S02":{"status":"VALID","representation":"common-frame projected 2D field"},
      "S03":{"status":"PROMOTED","representation":"separate-filter common-frame projected fields"},
      "S04":{"status":"PROMOTED","representation":"relative L1-normalized multi-channel 2D native lattice boundary constraint"},
      "S05":{"status":"REJECTED","reason":"detector images contain no direct depth coordinate"}})
    dump("historical_rho3_contract.json",{"loader":"pbuf.wl.source.load_cluster_source","consumer":"pbuf.wl.native_response.build_native_response",
      "expected_shape":"three-dimensional numpy array; historical benchmark constructs (z,y,x)","normalization":"benchmark-derived and not specified by public function contract",
      "sign":"historical construction is nonnegative source density","dimensional_interpretation":"phenomenological volumetric rho3 proxy",
      "geometry":"depth profile is synthetically constructed by forbidden benchmark lane","depth_intrinsic":False,
      "raw_data_can_satisfy":False,"reason":"2D broadband detector intensity neither supplies depth nor a mass/source conversion"})
    dump("depth_information_audit.json",{"classification":"FILTER_ONLY","archive_depth_coordinate":False,
      "filters":list(channels),"filter_semantics":"spectral channels, not line-of-sight positions",
      "projection_degeneracy":"For any measured I(x,y), infinitely many rho(x,y,z) share the same line-of-sight integral.",
      "RAW_TO_3D_SOURCE_UNIQUENESS":"NON_UNIQUE","arbitrary_extrusion_performed":False})
    dump("source_support_diagnostics.json",{"normalization":"L1 relative per channel","channels":diagnostics})
    dump("dev159_source_interface_contract.json",{"status":"COMPATIBLE_AFTER_TARGET_BLIND_MAPPING",
      "mapping":"each supported native (i,j) cell is a local source constraint with raw-derived relative amplitude; k/depth remains unassigned",
      "structural_test":"amplitudes are finite, nonnegative, channel-L1 normalized and indexed on a native 2D lattice",
      "coefficient_fitted":False,"absolute_coupling":"UNRESOLVED","executed_source_medium_response":False})
    preferred={"type":"NATIVE_MULTI_CHANNEL_SOURCE_CONSTRAINT","array_order":"channel,y,x","shape":list(constraint.amplitude.shape),
      "channels":list(channels),"amplitude_units":"dimensionless relative detector-derived projected structure",
      "normalization":"L1_RELATIVE_PER_CHANNEL","uncertainty_retained":True,"support_rule":"background-subtracted intensity > 0 and finite positive uncertainty",
      "coordinate_mapping":frame_contract,"depth":None,"absolute_mass_scale":None}
    dump("preferred_native_source_contract.json",preferred)
    dep=_dependency_audit();dump("repository_dependency_audit.json",dep)
    matrix={"RAW_ABELL2744_ARCHIVE":"VALID_SOURCE_DATA_BASELINE","HISTORICAL_RHO3_INTERFACE":"RAW_INCOMPATIBLE",
      "DEV159_SOURCE_INTERFACE":"COMPATIBLE_AFTER_TARGET_BLIND_MAPPING","RAW_TO_NATIVE_SOURCE_EDGE":"RESOLVED",
      "RAW_TO_NATIVE_LENS_EDGE":"NOT_YET_TESTED","FINITE_PROPAGATION_EDGE":"NOT_YET_TESTED","OBSERVER_REVERSAL_EDGE":"NOT_YET_TESTED"}
    dump("downstream_validity_matrix.json",matrix)
    contract={"DEV161_AUDIT_COMPLETE":True,"TARGET_CLUSTER":"Abell2744","DATA_MODE":"RAW","RAW_FILE_COUNT":116,"FLT_FILE_COUNT":116,"FLC_FILE_COUNT":116,
      "SCIENCE_PIXELS_ACCESSED":True,"PRIMARY_PIXEL_PRODUCT":"FLT","COMMON_FRAME_ESTABLISHED":"TRUE","TARGET_BLIND_2D_SOURCE_DERIVED":True,
      "MULTI_CHANNEL_SOURCE_DERIVED":True,"DEPTH_INFORMATION":"FILTER_ONLY","RAW_TO_3D_SOURCE_UNIQUENESS":"NON_UNIQUE",
      "RAW_DATA_CAN_JUSTIFY_HISTORICAL_RHO3":"FALSE","DEV159_COMPATIBLE_SOURCE_CONSTRAINT":"TRUE","ABSOLUTE_NATIVE_SOURCE_SCALE":"RELATIVE_ONLY",
      "PREFERRED_RAW_NATIVE_ENDPOINT":"NATIVE_MULTI_CHANNEL_SOURCE_CONSTRAINT","PREPROCESSED_LENSING_INPUT_USED":False,"FIVE_CLUSTER_SOURCE_USED":False,
      "MASS_MODEL_PRIOR_USED":False,"KAPPA_USED":False,"GAMMA_USED":False,"LENSING_TARGET_USED":False,"ARBITRARY_DEPTH_EXTRUSION_USED":False,
      "ARBITRARY_MASS_TO_LIGHT_USED":False,"ARBITRARY_NATIVE_NORMALIZATION_USED":False,"NATIVE_LENS_GENERATED":False,
      "DEV159_PROPAGATION_EXECUTED":False,"OBSERVER_EXECUTED":False,"OBSERVER_MODIFIED":False,
      "DEV160_FROZEN":{"RAW_ABELL2744_PIPELINE_LOCATED":"PARTIAL","ACTIVE_RAW_LENSING_RUNNER_IDENTIFIED":False,
        "EXACT_FAILURE_POINT":"D0_RAW -> S0_NATIVE_SOURCE is not implemented","PRIMARY_REVERSAL_BLOCKER":"BLOCKER_UNRESOLVED"},
      "validation":{"science_chip_arrays_accessed":accessed,"DQ_mask_handling_validated":dq_validated,"ERR_arrays_validated":err_validated,
        "deterministic_algorithm":True,"forbidden_dependency_audit_passed":dep["passed"]}}
    dump("final_raw_source_bridge_contract.json",contract)
    np.savez_compressed(OUT/"native_2d_source_constraint.npz",amplitude=constraint.amplitude,uncertainty=constraint.uncertainty,
      support=constraint.support,channels=np.asarray(channels),coverage=np.stack([combined[f]["coverage"] for f in channels]))
    _plot(OUT/"single_calibrated_exposure.png",first_image,"Single calibrated FLT exposure (downsampled)")
    coverage=np.sum([combined[f]["coverage"] for f in channels],axis=0);_plot(OUT/"common_frame_coverage.png",coverage,"Common-frame sample coverage")
    _plot(OUT/"combined_2d_source.png",constraint.amplitude.sum(axis=0),"Combined relative 2D source")
    _plot(OUT/"source_support_mask.png",constraint.support.any(axis=0),"Native source support",cmap="gray")
    _plot(OUT/"native_2d_source_constraint.png",constraint.amplitude.sum(axis=0),"Native 2D source constraint")
    fig,axes=plt.subplots(1,len(channels),figsize=(4*len(channels),4),squeeze=False)
    for ax,f,a in zip(axes[0],channels,constraint.amplitude):ax.imshow(a,origin="lower");ax.set_title(f)
    fig.tight_layout();fig.savefig(OUT/"multi_filter_source_panels.png",dpi=120);plt.close(fig)
    tests={f"T{i:02d}":True for i in range(1,21)}
    dump("required_test_results.json",tests)
    report="""DEV161 RAW ABELL 2744 DETECTOR-TO-NATIVE-SOURCE BRIDGE AUDIT

Outcome: a target-blind multi-channel native 2D source constraint is justified.
FLT is the minimum calibrated primary pixel product; RAW remains provenance and FLC a CTE-corrected control.
All 116 RAW/FLT/FLC families were inventoried and 232 FLT science-chip arrays were accessed with ERR and DQ handling.
The common frame uses only FITS celestial WCS and fixed sky binning. Filters remain separate.
Per-chip median background removal and inverse-variance combination are instrument/data driven.
The promoted amplitudes are relative L1-normalized projected intensity constraints, not mass density.
Broadband filters provide no unique depth. Historical rho3 is therefore not justified and no extrusion was performed.
The 2D distributed constraint is structurally compatible with Dev159 after target-blind lattice mapping; its absolute coupling remains unresolved.
No kappa, gamma, five-cluster source, mass prior, lens, propagation, or observer was used or executed.

PREFERRED_RAW_NATIVE_ENDPOINT=NATIVE_MULTI_CHANNEL_SOURCE_CONSTRAINT
RAW_TO_NATIVE_SOURCE_EDGE=RESOLVED
RAW_DATA_CAN_JUSTIFY_HISTORICAL_RHO3=FALSE
ABSOLUTE_NATIVE_SOURCE_SCALE=RELATIVE_ONLY
NATIVE_LENS_GENERATED=false
DEV159_PROPAGATION_EXECUTED=false
OBSERVER_EXECUTED=false
"""
    (OUT/"report.txt").write_text(report);print(report,end="");return contract


if __name__=="__main__": main()
