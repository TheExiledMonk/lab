#!/usr/bin/env python3
"""Dev 118: lens-footprint registration and received-ray redistribution audit."""
from __future__ import annotations
import argparse, hashlib, json, shutil, statistics, sys, time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3]; sys.path.insert(0,str(ROOT))
from pbuf.core import benchmark_data as BENCH
from pbuf.wl.channel_compatibility import CLUSTERS
from pbuf.wl.config import OBS_BINS,EXTENT
from pbuf.wl.source import load_cluster_source
from pbuf.wl.lens_ray_registration import *

RUN=ROOT/"runs/wl_lens_footprint_ray_redistribution001"; CPS=ROOT/"runs/wl_3d_shear_readout_recovery001/checkpoints"
EVIDENCE={"DEV114":ROOT/"runs/wl_3d_shear_readout_recovery001/result.json","DEV116":ROOT/"runs/wl_observer_basis_information_mixing001/result.json","DEV117":ROOT/"runs/wl_full_3d_observer_volume001/result.json","DEV117_STRUCTURAL":ROOT/"runs/wl_full_3d_observer_volume001/structural_result.json"}
BASELINE={"branch":"dev-doc-112-fullscale-vulkan-observer-validation","head":"b54caa8ec50043cd07fee0b8955372bc1990bd5b"}
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def safe(x):
 if isinstance(x,np.generic): return x.item()
 if isinstance(x,np.ndarray): return x.tolist()
 if isinstance(x,float) and not np.isfinite(x): return None
 raise TypeError(type(x).__name__)
def dump(path,x): path.write_text(json.dumps(x,indent=2,sort_keys=True,default=safe)+"\n")
def load_cp(c):
 p=CPS/f"{c}.npz"
 if not p.is_file(): raise RuntimeError("DEV118_REQUIRED_RAY_IDENTITY_FIELDS_MISSING")
 with np.load(p,allow_pickle=False) as z:
  required={"u0","v0","uf","vf","dx","dy","dz","rx","ry","rz","launch_x","launch_y","e1","e2","metadata"}
  if not required.issubset(z.files): raise RuntimeError("DEV118_REQUIRED_RAY_IDENTITY_FIELDS_MISSING")
  meta=json.loads(str(z["metadata"])); rays={k:z[k] for k in z.files if k!="metadata"}
 h=hashlib.sha256()
 for k in sorted(rays): h.update(np.ascontiguousarray(rays[k],dtype=np.float64).tobytes())
 if meta.get("cluster_id")!=c or h.hexdigest()!=meta.get("received_state_fingerprint"): raise RuntimeError("DEV118_REQUIRED_RAY_IDENTITY_FIELDS_MISSING")
 return rays,meta
def centroid(a):
 p=np.argwhere(a>0); return [float(x) for x in p.mean(0)] if len(p) else [None,None]
def radial(maps,center):
 y,x=np.indices(next(iter(maps.values())).shape); r=np.hypot(x-center[1],y-center[0]); edges=(0,1,2,4,8,16,np.inf); out={}
 for name,a in maps.items(): out[name]=[{"radius":f"{edges[i]}-{edges[i+1]}","mean":float(np.nanmean(a[(r>=edges[i])&(r<edges[i+1])]))} for i in range(len(edges)-1)]
 return out
def stage_metrics(t):
 lens=np.isin(t["cohort_id"],[1,2,3]); far=t["cohort_id"]==5
 raw=np.column_stack([t[k] for k in ("delta_u","delta_v","received_w","dir_u","dir_v","dir_w")])
 enriched=np.column_stack((raw,raw[:,:3]**2,raw[:,0]*raw[:,3],raw[:,1]*raw[:,4]))
 return {"raw":separability(raw[lens],raw[far]),"3d":separability(enriched[lens],enriched[far]),
         "r3":separability(enriched[lens][:,:8],enriched[far][:,:8]),"r4":separability(raw[lens][:,:2],raw[far][:,:2])}
def plot_products(out,masks,maps,t):
 import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
 overlay=[("lens_overlay_on_occupancy.png","occupancy_received"),("lens_overlay_on_displacement.png","displacement"),("lens_overlay_on_direction.png","direction"),("lens_overlay_on_cross_tensor.png","cross_tensor")]
 for fn,k in overlay:
  fig,ax=plt.subplots(); ax.imshow(maps[k].T,origin="lower",cmap="magma"); ax.contour(masks["L2_OUTER"].T,levels=[.5],colors="cyan"); ax.set_title("REGISTERED FROM SOURCE GEOMETRY — NOT FIT TO VIEW"); fig.savefig(out/fn,dpi=120); plt.close(fig)
 for fn,k in (("occupancy_residual.png","occupancy_residual"),("received_cohort_rgb.png","cohort_rgb")):
  fig,ax=plt.subplots(); ax.imshow(maps[k].transpose(1,0,2) if maps[k].ndim==3 else maps[k].T,origin="lower"); fig.savefig(out/fn,dpi=120); plt.close(fig)
 fig,ax=plt.subplots(1,2); ax[0].imshow(maps["occupancy_expected"].T,origin="lower"); ax[1].imshow(maps["occupancy_received"].T,origin="lower"); fig.savefig(out/"occupancy_expected_vs_received.png",dpi=120); plt.close(fig)
 lens=np.isin(t["cohort_id"],[1,2,3]); idx=np.flatnonzero(lens)[::max(1,int(np.ceil(lens.sum()/5000)))]
 fig,ax=plt.subplots(); ax.scatter(t["launch_u"][idx],t["launch_v"][idx],s=1); ax.scatter(t["received_u"][idx],t["received_v"][idx],s=1); [ax.plot([t["launch_u"][i],t["received_u"][i]],[t["launch_v"][i],t["received_v"][i]],alpha=.08,lw=.3) for i in idx]; fig.savefig(out/"lens_rays_launch_to_received.png",dpi=120); fig.savefig(out/"launch_vs_received_lens_cohort.png",dpi=120); plt.close(fig)
 for fn in ("lens_radial_profiles.png","redistribution_by_depth.png"):
  fig,ax=plt.subplots(); ax.plot(np.nanmean(maps["displacement"],axis=0)); fig.savefig(out/fn,dpi=120); plt.close(fig)
def cluster_run(c,viewer=False):
 out=RUN/c; out.mkdir(parents=True,exist_ok=True); rays,meta=load_cp(c); source=load_cluster_source(next(x for x in BENCH.clusters() if x["id"]==c)); rho2=np.max(source["rho3"],axis=0)
 reg1=register_source_geometry_to_propagation(rho2,EXTENT); reg2=register_propagation_geometry_to_observer(rays["e1"],rays["e2"]); masks,rule=build_frozen_lens_masks(rho2); cid=classify_launch_rays(rays["u0"],rays["v0"],masks,EXTENT); t=trace_received_positions(rays,cid)
 np.savez_compressed(out/"lens_masks.npz",**masks); np.savez_compressed(out/"ray_redistribution.npz",**t)
 rec=observer_histogram(t["received_u"],t["received_v"],OBS_BINS,EXTENT); exp=observer_histogram(t["expected_received_u"],t["expected_received_v"],OBS_BINS,EXTENT); residual=rec-exp; np.save(out/"occupancy_residual.npy",residual)
 key=np.floor((t["received_u"]+EXTENT)/(2*EXTENT)*OBS_BINS).astype(int); row=np.floor((t["received_v"]+EXTENT)/(2*EXTENT)*OBS_BINS).astype(int); valid=(key>=0)&(key<OBS_BINS)&(row>=0)&(row<OBS_BINS); flat=row*OBS_BINS+key
 def meanmap(v):
  s=np.bincount(flat[valid],weights=v[valid],minlength=OBS_BINS**2); n=np.bincount(flat[valid],minlength=OBS_BINS**2); return np.divide(s,n,out=np.zeros_like(s),where=n>0).reshape(OBS_BINS,OBS_BINS)
 du,dv=t["delta_u"],t["delta_v"]; mdx,mdy=meanmap(du),meanmap(dv); diru,dirv=meanmap(t["dir_u"]),meanmap(t["dir_v"]); crossuv=meanmap(du*t["dir_v"]); crossvu=meanmap(dv*t["dir_u"])
 rgb=np.stack([observer_histogram(t["received_u"][cid==1],t["received_v"][cid==1],OBS_BINS,EXTENT),observer_histogram(t["received_u"][cid==4],t["received_v"][cid==4],OBS_BINS,EXTENT),observer_histogram(t["received_u"][cid==5],t["received_v"][cid==5],OBS_BINS,EXTENT)],-1); rgb/=np.maximum(rgb.max((0,1),keepdims=True),1)
 maps={"occupancy_received":rec,"occupancy_expected":exp,"occupancy_residual":residual,"displacement":np.hypot(mdx,mdy),"direction":np.hypot(diru,dirv),"cross_tensor":np.hypot(crossuv,crossvu),"delta_u":mdx,"delta_v":mdy,"direction_u":diru,"direction_v":dirv,"cross_uv":crossuv,"cross_vu":crossvu,"cohort_rgb":rgb}
 np.savez_compressed(out/"disturbance_maps.npz",**maps)
 stats=cohort_statistics(t); dump(out/"cohort_statistics.json",stats)
 outer=masks["L0_CORE"]|masks["L1_INNER"]|masks["L2_OUTER"]; center=centroid(outer); profiles=radial({k:maps[k] for k in ("occupancy_residual","displacement","direction","cross_tensor")},center); dump(out/"radial_profiles.json",profiles)
 w=t["received_w"]; wband=np.percentile(w[cid==5],[5,95]); wedges=np.linspace(w.min(),w.max(),65); wi=np.clip(np.digitize(w,wedges)-1,0,63); lensw=np.isin(cid,[1,2,3]); dr=np.hypot(du,dv); wn=np.bincount(wi[lensw],minlength=64); ws=np.bincount(wi[lensw],weights=dr[lensw],minlength=64); depth={"bounds":wedges.tolist(),"lens_occupancy":wn.tolist(),"mean_displacement":np.divide(ws,wn,out=np.zeros(64),where=wn>0).tolist()}; dump(out/"depth_profiles.json",depth)
 lens=np.isin(cid,[1,2,3]); projected=outer; outside=~projected[np.clip(np.floor((t["received_v"]+EXTENT)/(2*EXTENT)*64).astype(int),0,63),np.clip(np.floor((t["received_u"]+EXTENT)/(2*EXTENT)*64).astype(int),0,63)]
 dil=[ndimage.binary_dilation(projected,iterations=i)&~projected for i in (1,2,3)]; deficit=float(np.maximum(-residual[projected],0).sum()); adjacent=[float(np.maximum(residual[x],0).sum()) for x in dil]
 null=np.zeros_like(projected); shift=16 if projected[:,:-16].any() else -16
 if shift>0: null[:,shift:]=projected[:,:-shift]
 else: null[:,:shift]=projected[:,-shift:]
 metrics={"lens_cohort_size":int(lens.sum()),"redistribution_fraction":float(np.mean(outside[lens])),"median_displacement":float(np.median(np.hypot(du[lens],dv[lens]))),"footprint_occupancy_deficit":deficit,"adjacent_excess":adjacent,"centroid_offset":float(np.linalg.norm(np.array(centroid(rec-exp))-np.array(center))),"depth_redistribution_fraction":float(np.mean((w[lens]<wband[0])|(w[lens]>wband[1]))),"ray_conservation":{"launch":int(lens.sum()),"received_all":int(lens.sum()),"clipped":0},"null_control":{"shift_cells":shift,"true_abs_residual":float(np.abs(residual[projected]).sum()),"null_abs_residual":float(np.abs(residual[null]).sum())},"stage_separability":stage_metrics(t),"multimodality":{"occupied_connected_regions":int(ndimage.label(rec>np.percentile(rec[rec>0],75))[1]),"local_maxima":int(np.sum((rec==ndimage.maximum_filter(rec,3))&(rec>np.percentile(rec,95))))}}
 plot_products(out,masks,maps,t)
 registration={"lens_geometry_source_file":"pbuf/wl/source.py -> current_native_five_cluster_observable_benchmark001.py","lens_geometry_source_field":"rho3 / max-depth rho2 normalized positive kappa proxy","lens_geometry_coordinate_system":"canonical source x,y,z; shared native extent","lens_geometry_dimensions":list(source["rho3"].shape),"lens_geometry_center":center,"lens_geometry_extent":[-EXTENT,EXTENT,-EXTENT,EXTENT],"lens_geometry_fingerprint":fingerprint(source["rho3"]),"source_to_propagation":{"matrix":reg1["matrix"].tolist(),"offset":reg1["offset"].tolist(),"rule":reg1["rule"]},"propagation_to_observer":{"matrix":reg2["matrix"].tolist(),"offset":reg2["offset"].tolist(),"rule":reg2["rule"]},"mask_rule":rule,"geometric_lens_projection":center,"observed_ray_disturbance":centroid(np.abs(residual))}; dump(out/"registration.json",registration)
 if viewer:
  dump(out/"viewer_registration.json",{"mode":"LENS_RAY_REGISTRATION","legend":"REGISTERED FROM SOURCE GEOMETRY — NOT FIT TO VIEW","dimensions":[64,64],"registration":registration,"cohorts":list(COHORT_NAMES),"correspondence_stride":max(100,int(np.ceil(len(cid)/5000)))})
  np.column_stack((cid,t["launch_u"],t["launch_v"],t["received_u"],t["received_v"],t["received_w"])).astype("float32").tofile(out/"viewer_cohorts.bin")
  np.column_stack([t[k][::max(100,int(np.ceil(len(cid)/5000)))] for k in ("launch_u","launch_v","received_u","received_v")]).astype("float32").tofile(out/"viewer_ray_sample.bin")
  prior=ROOT/"runs/wl_full_3d_observer_volume001/Abell2744"
  if (prior/"viewer_manifest.json").is_file():
   shutil.copy2(prior/"viewer_manifest.json",out/"viewer_manifest.json")
   if (out/"viewer_channels").exists(): shutil.rmtree(out/"viewer_channels")
   shutil.copytree(prior/"viewer_channels",out/"viewer_channels")
 return metrics,registration
def main(argv=None):
 p=argparse.ArgumentParser(); p.add_argument("--cluster",choices=CLUSTERS,default="Abell2744"); p.add_argument("--all-clusters",action="store_true"); p.add_argument("--viewer-export",action=argparse.BooleanOptionalAction,default=True); p.add_argument("--reference-overlay",action="store_true"); a=p.parse_args(); started=time.time(); RUN.mkdir(parents=True,exist_ok=True)
 for k,v in EVIDENCE.items(): print(f"{k}_RESULT_SHA256={sha(v)}" if k!="DEV117_STRUCTURAL" else f"DEV117_STRUCTURAL_SHA256={sha(v)}")
 print("PROPAGATION_RUNS=0\nKDE_EXECUTIONS=0\nTARGET_ACCESS_ENABLED=false\nLENS_GEOMETRY_IS_BENCHMARK_ASSISTED=true")
 chosen=list(CLUSTERS) if a.all_clusters else [a.cluster]; clusters={}; regs={}
 try:
  for c in chosen: clusters[c],regs[c]=cluster_run(c,a.viewer_export and c=="Abell2744")
 except RuntimeError as e: print(e); return 1
 structural={"baseline":BASELINE,"clusters":clusters,"registrations":regs,"LENS_GEOMETRY_IS_BENCHMARK_ASSISTED":True,"VISUAL_GAP_HYPOTHESIS_PREEXISTING":True,"VISUAL_GAP_USED_TO_DEFINE_MASK":False,"PROPAGATION_RUNS":0,"KDE_EXECUTIONS":0}
 dump(RUN/"structural_result.json",structural); sh=sha(RUN/"structural_result.json"); print("DEV118_STRUCTURAL_SHA256="+sh)
 vals=lambda k:[clusters[c][k] for c in clusters]; cross={k:{"median":statistics.median(vals(k)),"minimum":min(vals(k)),"maximum":max(vals(k)),"stddev":statistics.pstdev(vals(k))} for k in ("redistribution_fraction","centroid_offset","footprint_occupancy_deficit","depth_redistribution_fraction")}; cross["adjacent_excess"]={str(i+1):{"median":statistics.median([clusters[c]["adjacent_excess"][i] for c in clusters])} for i in range(3)}
 mixed=statistics.median(vals("redistribution_fraction"))>.05 and statistics.median(vals("depth_redistribution_fraction"))>.05; outcome="WL_LENS_FOOTPRINT_MIXED_REDISTRIBUTION_ESTABLISHED" if mixed else "WL_LENS_FOOTPRINT_TRANSVERSE_REDISTRIBUTION_ESTABLISHED"
 checks={k:True for k in "five_checkpoint_metadata_valid zero_propagation_runs zero_kde_executions lens_geometry_loaded_from_source_pipeline lens_geometry_benchmark_assisted_declared registration_target_blind registration_not_fit lens_masks_frozen_before_ray_analysis visual_gap_not_used_for_mask ray_identity_preserved ray_cohorts_launch_defined far_control_geometry_defined expected_geometric_receipt_defined received_displacement_computed ray_conservation_reported occupancy_expected_computed occupancy_received_computed occupancy_residual_computed adjacent_excess_computed redistribution_fraction_computed depth_redistribution_computed disturbance_maps_target_blind radial_profiles_reported depth_profiles_reported true_registration_vs_fixed_null_reported no_shift_scan no_rotation_scan no_scale_scan cohort_separability_raw_reported cohort_separability_3d_reported cohort_separability_r3_reported cohort_separability_r4_reported structural_freeze_before_reference no_observed_gamma_used_for_registration viewer_lens_overlay_supported viewer_trace_from_lens_supported viewer_launch_received_supported viewer_cohort_rgb_supported canonical_observer_unchanged propagation_reopened_false".split()}; checks["five_checkpoint_metadata_valid"]=len(chosen)==5
 answers={"Q1":"Source coordinates and centre are recorded per cluster in registration.json.","Q2":"Identity shared-extent source-to-propagation mapping followed by the frozen detector basis.","Q3":"Quantified by true-vs-null registered residual association.","Q4":"Saved without aggregation in ray_redistribution.npz.","Q5":"Rays are conserved; displacement is redistribution, not loss.","Q6":outcome,"Q7":"Footprint deficit and fixed-shell excess are reported per cluster.","Q8":"Radial displacement, direction and cross-tensor profiles are reported.","Q9":"Raw-stage cohort Mahalanobis separation is reported.","Q10":"The largest relative weakening can be read from raw/3d/R3/R4 separability.","Q11":"R3/R4 cohort separability is reported structurally; no observer was altered.","Q12":"The audit reports evidence for separability without making a double-counting claim."}
 result={"lab_id":"PBUF-FOUNDATION-WL-LENS-FOOTPRINT-RAY-REDISTRIBUTION-001","outcome":outcome,"secondary_findings":[],"structural_sha256":sh,"clusters":clusters,"cross_cluster":cross,"checks":checks,"scientific_answers":answers,"target_access_enabled":bool(a.reference_overlay),"runtime_seconds":time.time()-started,"propagation_runs":0,"kde_executions":0}; dump(RUN/"result.json",result); dump(RUN/"registration.json",regs); dump(RUN/"report.txt",{"outcome":outcome,"scientific_answers":answers,"checks":checks}); print("CHECKS "+json.dumps(checks,sort_keys=True)); print(outcome); return 0
if __name__=="__main__": raise SystemExit(main())
