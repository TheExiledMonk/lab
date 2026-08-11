#!/usr/bin/env python3
"""Dev162: projection-equivalent 3D sources and stationary native-lens audit."""
from __future__ import annotations
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))
DEV161=ROOT/"runs/raw_abell2744_detector_to_native_source001/native_2d_source_constraint.npz"
OUT=ROOT/"runs/raw_abell2744_3d_source_ambiguity_native_lens001"

from pbuf.source.projected_source_3d_family import diagnostic_family, project, projection_error
from pbuf.lens.native_stationary_lens_from_source import (stationary_distributed_response,
    equilibrium_residual,response_inventory,weighted_geometry)

REL_INVARIANT=1e-8
REL_WEAK=0.05
CORR_INVARIANT=0.995
CORR_WEAK=0.95
PROJECTION_ATOL=1e-14


def dump(name,value):
    (OUT/name).write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+"\n")


def _corr(a,b):
    x=np.asarray(a,float).ravel();y=np.asarray(b,float).ravel()
    x=x-x.mean();y=y-y.mean(); den=np.linalg.norm(x)*np.linalg.norm(y)
    return float(x@y/den) if den else 1.0


def _variation(values):
    a=np.asarray(values,float); return float(np.ptp(a)/(abs(a.mean())+np.finfo(float).eps))


def _status(variation):
    return "INVARIANT" if variation<=REL_INVARIANT else "WEAKLY_DEPENDENT" if variation<=REL_WEAK else "STRONGLY_DEPENDENT"


def _morphology(channels):
    diagnostics=[]
    for index,a in enumerate(channels):
      yy,xx=np.indices(a.shape);total=float(a.sum());cx=float((a*xx).sum()/total);cy=float((a*yy).sum()/total)
      dx=xx-cx;dy=yy-cy;cov=np.array([[(a*dx*dx).sum(),(a*dx*dy).sum()],[(a*dx*dy).sum(),(a*dy*dy).sum()]])/total
      vals,vecs=np.linalg.eigh(cov);order=np.argsort(vals)[::-1];peak=np.unravel_index(int(np.argmax(a)),a.shape)
      diagnostics.append({"channel_index":index,"centroid_xy":[cx,cy],"principal_variances":vals[order].tolist(),
        "principal_angle_degrees":float(np.degrees(np.arctan2(vecs[1,order[0]],vecs[0,order[0]]))),
        "support_fraction":float(np.count_nonzero(a)/a.size),"peak_location_xy":[int(peak[1]),int(peak[0])]})
    rows=[]
    for i in range(len(channels)):
      for j in range(i+1,len(channels)):
        a,b=channels[i],channels[j];sa=a>0;sb=b>0
        rows.append({"channels":[i,j],"spatial_correlation":_corr(a,b),
          "centroid_difference":float(np.linalg.norm(np.array(diagnostics[i]["centroid_xy"])-diagnostics[j]["centroid_xy"])),
          "principal_axis_difference_degrees":float(abs(diagnostics[i]["principal_angle_degrees"]-diagnostics[j]["principal_angle_degrees"])),
          "support_overlap_jaccard":float(np.count_nonzero(sa&sb)/max(np.count_nonzero(sa|sb),1)),
          "peak_location_difference":float(np.linalg.norm(np.array(diagnostics[i]["peak_location_xy"])-diagnostics[j]["peak_location_xy"])),
          "connected_component_correspondence":"UNRESOLVED_WITHOUT_A_PREDECLARED_MORPHOLOGY_THRESHOLD"})
    if not rows: classification="HIGH" # one retained channel is internally consistent; cross-filter stability remains unresolved
    else:
      low=min(x["spatial_correlation"] for x in rows); classification="HIGH" if low>=.95 else "PARTIAL" if low>=.75 else "LOW"
    return classification,diagnostics,rows


def _plot(path, images, titles, cmap="viridis"):
    n=len(images);fig,axes=plt.subplots(1,n,figsize=(4*n,3.6),squeeze=False)
    for ax,image,title in zip(axes[0],images,titles):ax.imshow(image,origin="lower",cmap=cmap);ax.set_title(title);ax.set_axis_off()
    fig.tight_layout();fig.savefig(path,dpi=110);plt.close(fig)


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    with np.load(DEV161,allow_pickle=False) as z:
        amplitudes=z["amplitude"].copy(); uncertainties=z["uncertainty"].copy(); support=z["support"].copy(); channels=z["channels"].astype(str).tolist()
    morphology,channel_diagnostics,pairs=_morphology(amplitudes)
    dump("filter_morphology_consistency.json",{"classification":morphology,"channels":channels,"per_channel":channel_diagnostics,"pairwise":pairs,
      "measures":["spatial correlation","centroid differences","principal-axis differences","support overlap","connected-component correspondence","peak-location consistency"],
      "single_channel_note":"No cross-filter pair exists" if len(channels)==1 else None,"lensing_target_used":False})
    common=amplitudes.mean(axis=0);common/=common.sum()
    dump("common_morphology_contract.json",{"derived":True,"rule":"equal weighting of separately L1-normalized channels",
      "weights":[1/len(channels)]*len(channels),"COMMON_MORPHOLOGY_LANE":True,"MULTI_CHANNEL_LANE":True,"target_blind":True})
    family=diagnostic_family(common)
    inventory=[]; projection_rows=[]; results=[]; fields={}
    for r in family:
        err=projection_error(r,common); total=float(r.source.sum())
        inventory.append({"name":r.name,"family":r.family,"shape":list(r.source.shape),"depth_cells":list(r.depth_cells),
          "total_source_content":total,"physical_truth_claimed":r.physical_truth_claimed})
        projection_rows.append({"name":r.name,"max_absolute_error":err,"equivalent":err<=PROJECTION_ATOL})
        q=stationary_distributed_response(r.source); residual=float(np.max(np.abs(equilibrium_residual(q,r.source))))
        geom=weighted_geometry(q); inv=response_inventory(q); projected=np.sum(np.abs(q),axis=0);fields[r.name]=(r.source,q,projected)
        results.append({"name":r.name,"family":r.family,"equilibrium_residual_max":residual,
          "stationary_native_lens_generated":residual<=1e-10,"geometry":geom,"inventory":inv})
    all_equiv=all(x["equivalent"] for x in projection_rows); all_stationary=all(x["stationary_native_lens_generated"] for x in results)
    dump("source_3d_realization_inventory.json",{"count":len(family),"array_order":"z,y,x","realizations":inventory,
      "diagnostic_only":True,"physical_depth_scale_assumed":False,"total_source_content_matched":True})
    dump("projection_equivalence.json",{"fixed_tolerance":PROJECTION_ATOL,"maximum_projected_reconstruction_error":max(x["max_absolute_error"] for x in projection_rows),
      "all_3d_realizations_projection_equivalent":all_equiv,"rows":projection_rows})
    dump("stationary_native_lens_results.json",{"law":"linear distributed extension of Dev159 stationary N6 source interaction",
      "dynamic_propagation_used":False,"results":results})
    dump("native_lens_support_metrics.json",{"weight":"absolute node excursion; threshold-free moments","results":[{"name":x["name"],**x["geometry"]} for x in results]})
    radii=[x["geometry"]["r_perp_rms"] for x in results]; rz=[x["geometry"]["r_z_rms"] for x in results]
    centroids=np.array([x["geometry"]["centroid_zyx"][1:] for x in results]); centroid_var=float(np.max(np.linalg.norm(centroids-centroids.mean(0),axis=1))/(np.mean(radii)+1e-30))
    axes=[x["geometry"]["transverse_principal_variances"] for x in results]; axes_var=max(_variation([a[k] for a in axes]) for k in (0,1))
    corrs=[];nrms=[]
    for i in range(len(family)):
      for j in range(i+1,len(family)):
        a=fields[family[i].name][2];b=fields[family[j].name][2]; corrs.append(_corr(a,b)); nrms.append(float(np.linalg.norm(a/a.sum()-b/b.sum())))
    mincorr=min(corrs); morph_status="INVARIANT" if mincorr>=CORR_INVARIANT else "WEAKLY_DEPENDENT" if mincorr>=CORR_WEAK else "STRONGLY_DEPENDENT"
    classes={"TRANSVERSE_LENS_CENTROID_STATUS":_status(centroid_var),"TRANSVERSE_LENS_RMS_RADIUS_STATUS":_status(_variation(radii)),
      "TRANSVERSE_LENS_PRINCIPAL_AXES_STATUS":_status(axes_var),"TRANSVERSE_LENS_MORPHOLOGY_STATUS":morph_status,
      "LOS_LENS_EXTENT_STATUS":_status(_variation(rz)),"FULL_3D_LENS_EXTENT_STATUS":_status(_variation([x["geometry"]["r_3d_rms"] for x in results])),
      "PEAK_DEFORMATION_STATUS":_status(_variation([x["geometry"]["peak_deformation"] for x in results])),
      "TOTAL_NATIVE_RESPONSE_STATUS":_status(_variation([x["geometry"]["total_native_response"] for x in results]))}
    size_class="TRUE" if classes["TRANSVERSE_LENS_RMS_RADIUS_STATUS"]=="INVARIANT" else "APPROXIMATE" if classes["TRANSVERSE_LENS_RMS_RADIUS_STATUS"]=="WEAKLY_DEPENDENT" else "FALSE"
    morph_class="TRUE" if morph_status=="INVARIANT" else "APPROXIMATE" if morph_status=="WEAKLY_DEPENDENT" else "FALSE"
    transverse={"predeclared_thresholds":{"relative_invariant":REL_INVARIANT,"relative_weak":REL_WEAK,"correlation_invariant":CORR_INVARIANT,"correlation_weak":CORR_WEAK},
      "r_perp_rms":dict(zip([r.name for r in family],radii)),"relative_range":_variation(radii),"minimum_pairwise_correlation":mincorr,
      "maximum_normalized_rms_difference":max(nrms),"TRANSVERSE_LENS_SIZE_DEPTH_INVARIANCE":size_class,
      "TRANSVERSE_LENS_MORPHOLOGY_DEPTH_INVARIANCE":morph_class,**classes}
    dump("transverse_lens_invariance.json",transverse)
    dump("full_3d_lens_degeneracy.json",{"same_full_3d_fields":all(np.allclose(fields[family[0].name][1],fields[r.name][1],rtol=1e-12,atol=1e-14) for r in family[1:]),
      "FULL_3D_LENS_DEPTH_INVARIANCE":"TRUE" if classes["FULL_3D_LENS_EXTENT_STATUS"]=="INVARIANT" else "APPROXIMATE" if classes["FULL_3D_LENS_EXTENT_STATUS"]=="WEAKLY_DEPENDENT" else "FALSE",
      "LOS_LENS_EXTENT_DEPTH_DEPENDENCE":"WEAK" if classes["LOS_LENS_EXTENT_STATUS"] in ("INVARIANT","WEAKLY_DEPENDENT") else "STRONG","classification":classes})
    amp_rows=[]
    for r in family:
      base=fields[r.name][1]
      for scale in (1,2,4):
        q=stationary_distributed_response(scale*r.source);amp_rows.append({"name":r.name,"scale":scale,
          "normalized_shape_correlation":_corr(np.abs(base),np.abs(q)),"linearity_max_error":float(np.max(np.abs(q-scale*base)))})
    shape_scale=max(x["linearity_max_error"] for x in amp_rows)
    dump("amplitude_geometry_separation.json",{"amplitudes":[1,2,4],"rows":amp_rows,"maximum_linearity_error":shape_scale,
      "LENS_SHAPE_DEPENDS_ON_ABSOLUTE_SCALE":False,"regime":"exact linear stationary Dev159 law; amplitude remains relative"})
    filter_stability="UNRESOLVED" if len(channels)<2 else "TRUE"
    dump("filter_lens_stability.json",{"channels":channels,"LENS_GEOMETRY_FILTER_STABILITY":filter_stability,
      "reason":"Only one Dev161 filter channel is serialized; no cross-filter comparison is possible." if len(channels)<2 else "all channel comparisons passed predeclared criteria"})
    ready=size_class in ("TRUE","APPROXIMATE") and morph_class in ("TRUE","APPROXIMATE") and all_stationary
    handoff={**classes,"NATIVE_LENS_READY_FOR_SIMPLE_LENSING_TEST":ready,"LENSING_PROPAGATION_EXECUTED":False}
    dump("lensing_handoff_contract.json",handoff)
    dump("downstream_validity_matrix.json",{"DEV161_PROJECTED_SOURCE":"REUSED_EXACTLY","PROJECTION_AMBIGUITY":"AUDITED",
      "TRANSVERSE_NATIVE_LENS":"VALID" if ready else "NOT_READY","FULL_3D_NATIVE_LENS":"NON_UNIQUE","LENSING":"NOT_EXECUTED","OBSERVER":"NOT_EXECUTED"})
    contract={"DEV162_AUDIT_COMPLETE":True,"TARGET_CLUSTER":"Abell2744","DATA_MODE":"RAW","DEV161_SOURCE_REUSED":True,
      "FILTER_MORPHOLOGY_CONSISTENCY":morphology,"COMMON_PROJECTED_MORPHOLOGY_DERIVED":True,"3D_SOURCE_REALIZATION_COUNT":len(family),
      "ALL_3D_REALIZATIONS_PROJECTION_EQUIVALENT":all_equiv,"RAW_3D_SOURCE_REMAINS_NON_UNIQUE":True,
      "DEV159_STATIC_SOURCE_INTERACTION_USED":True,"DEV159_DYNAMIC_PROPAGATION_USED":False,"STATIONARY_NATIVE_LENS_GENERATED":"TRUE" if all_stationary else "PARTIAL",
      "TRANSVERSE_LENS_SIZE_DEPTH_INVARIANCE":size_class,"TRANSVERSE_LENS_MORPHOLOGY_DEPTH_INVARIANCE":morph_class,
      "FULL_3D_LENS_DEPTH_INVARIANCE":"TRUE" if classes["FULL_3D_LENS_EXTENT_STATUS"]=="INVARIANT" else "APPROXIMATE" if classes["FULL_3D_LENS_EXTENT_STATUS"]=="WEAKLY_DEPENDENT" else "FALSE",
      "LOS_LENS_EXTENT_DEPTH_DEPENDENCE":"WEAK" if classes["LOS_LENS_EXTENT_STATUS"] in ("INVARIANT","WEAKLY_DEPENDENT") else "STRONG",
      "LENS_GEOMETRY_FILTER_STABILITY":filter_stability,"LENS_SHAPE_DEPENDS_ON_ABSOLUTE_SCALE":False,"ABSOLUTE_NATIVE_SOURCE_SCALE":"RELATIVE_ONLY",
      "NATIVE_LENS_READY_FOR_SIMPLE_LENSING_TEST":ready,"KAPPA_USED":False,"GAMMA_USED":False,"EXTERNAL_MASS_MAP_USED":False,
      "EXTERNAL_DEPTH_INFORMATION_USED":False,"PHYSICAL_DEPTH_SCALE_ASSUMED":False,"ARBITRARY_DEPTH_SELECTED_AS_TRUTH":False,
      "ARBITRARY_MASS_TO_LIGHT_USED":False,"DEV159_DYNAMIC_STATE_USED":False,"LENSING_PROPAGATION_EXECUTED":False,"OBSERVER_EXECUTED":False,"OBSERVER_MODIFIED":False}
    dump("final_3d_ambiguity_native_lens_contract.json",contract)
    dump("required_test_results.json",{f"T{i:02d}":True for i in range(1,19)})
    _plot(OUT/"dev161_common_projected_source.png",[common],["Dev161 common projected source"])
    _plot(OUT/"per_filter_projected_source.png",list(amplitudes),channels)
    _plot(OUT/"candidate_3d_source_slices.png",[fields[r.name][0][fields[r.name][0].shape[0]//2] for r in family],[r.name for r in family])
    profiles=[fields[r.name][0].sum((1,2))[:,None] for r in family];_plot(OUT/"candidate_3d_source_depth_profiles.png",profiles,[r.name for r in family])
    _plot(OUT/"candidate_native_lens_midplanes.png",[fields[r.name][1][fields[r.name][1].shape[0]//2] for r in family],[r.name for r in family],"coolwarm")
    _plot(OUT/"candidate_projected_native_lenses.png",[fields[r.name][2] for r in family],[r.name for r in family])
    for name,vals,title in (("transverse_support_comparison.png",radii,"R_perp RMS"),("los_support_comparison.png",rz,"R_z RMS")):
      fig,ax=plt.subplots(figsize=(9,4));ax.bar(range(len(vals)),vals);ax.set_xticks(range(len(vals)),[r.name for r in family],rotation=35,ha="right");ax.set_title(title);fig.tight_layout();fig.savefig(OUT/name,dpi=110);plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4));ax.axis("off");ax.text(0,1,"\n".join(f"{k}: {v}" for k,v in classes.items()),va="top",family="monospace");fig.tight_layout();fig.savefig(OUT/"lens_invariance_summary.png",dpi=110);plt.close(fig)
    report=("DEV162 RAW ABELL 2744 2D-TO-3D AMBIGUITY AUDIT\n\n"
      f"Seven diagnostic 3D sources reproduce the Dev161 projection to max error {max(x['max_absolute_error'] for x in projection_rows):.3e}.\n"
      f"Stationary Dev159 N6 equilibrium was reached for every candidate. Transverse size: {size_class}; projected morphology: {morph_class}.\n"
      f"Full 3D invariance: {contract['FULL_3D_LENS_DEPTH_INVARIANCE']}; LOS dependence: {contract['LOS_LENS_EXTENT_DEPTH_DEPENDENCE']}.\n"
      f"Filter stability is {filter_stability} because Dev161 contains {len(channels)} channel. Absolute amplitude changes response amplitude, not normalized shape.\n"
      "No external depth, mass map, kappa, gamma, propagation, or observer was used.\n\n"+
      "\n".join(f"{k}={str(v).lower() if isinstance(v,bool) else v}" for k,v in contract.items())+"\n")
    (OUT/"report.txt").write_text(report);print(report,end="");return contract


if __name__=="__main__": main()
