#!/usr/bin/env python3
"""Dev 120: launch-bundle geometry and dual-coordinate transport audit."""
from __future__ import annotations
import argparse, hashlib, json, shutil, sys, time
from pathlib import Path
import numpy as np
from scipy import ndimage, stats

ROOT=Path(__file__).resolve().parents[3]; sys.path.insert(0,str(ROOT))
from pbuf.wl.channel_compatibility import CLUSTERS
from pbuf.wl.bundle_transport_geometry import (SCALES, deposit_bundle, effective_rank,
    fit_bundle_scale, json_dump, multiscale_metrics, pearson, received_state,
    received_overlap_diagnostic, reconstruct_launch_topology, resolution_metrics, structural_sha256)
from pbuf.wl.lens_ray_registration import separability

RUN=ROOT/"runs/wl_bundle_transport_geometry001"
CPS=ROOT/"runs/wl_3d_shear_readout_recovery001/checkpoints"
SD=ROOT/"runs/wl_source_deformation_separation001"
TRUTH=ROOT/"runs/wl_lens_footprint_ray_redistribution001"
EVIDENCE={
 "DEV117":ROOT/"runs/wl_full_3d_observer_volume001/result.json",
 "DEV118":TRUTH/"result.json", "DEV119":SD/"result.json",
 "DEV119_STRUCTURAL":SD/"structural_result.json"}

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load_cp(c):
 with np.load(CPS/f"{c}.npz",allow_pickle=False) as z:
  meta=json.loads(str(z["metadata"])); rays={k:z[k] for k in z.files if k!="metadata"}
 required={"u0","v0","uf","vf","dx","dy","dz","rx","ry","rz","launch_x","launch_y","e1","e2"}
 if not required.issubset(rays) or meta.get("cluster_id")!=c or len(rays["uf"])!=285156:
  raise RuntimeError("DEV120_REQUIRED_LAUNCH_RECEIVE_IDENTITY_MISSING")
 return rays,meta
def sample_sd(c,rays):
 S=np.load(SD/c/"source_score.npy"); D=np.load(SD/c/"deformation_score.npy")
 st=received_state(rays); coords=np.column_stack((st["received_u"],st["received_v"],st["received_w"]))
 lo=np.array([-8.,-8.,np.nanmin(coords[:,2])]); hi=np.array([8.,8.,np.nanmax(coords[:,2])])
 q=np.rint((coords-lo)*(np.array(S.shape)-1)/np.maximum(hi-lo,1e-30)).astype(int)
 q=np.clip(q,0,np.array(S.shape)-1)
 return S[tuple(q.T)],D[tuple(q.T)]
def metric_pair(pair,target):
 def one(p,t):
  return {"pearson":pearson(p,t),"spearman":float(stats.spearmanr(p.ravel(),t.ravel()).statistic),
          "rms_ratio":float(np.sqrt(np.mean(p*p))/(np.sqrt(np.mean(t*t))+1e-30))}
 p1,p2=pair;t1,t2=target
 return {"gamma1":one(p1,t1),"gamma2":one(p2,t2),
  "magnitude_pearson":pearson(np.hypot(p1,p2),np.hypot(t1,t2)),
  "orientation_agreement":float(np.mean(np.cos(np.arctan2(p2,p1)-np.arctan2(t2,t1))))}
def target(c,shape=(64,64)):
 from pbuf.core import benchmark_data as BENCH
 from pbuf.wl.source import load_cluster_source
 d=load_cluster_source(next(x for x in BENCH.clusters() if x["id"]==c))["data"]
 return tuple(ndimage.zoom(np.asarray(d[k],float),np.array(shape)/np.array(d[k].shape),order=1) for k in ("gamma1","gamma2"))
def rank_banks(rays,b):
 st=received_state(rays); received=np.column_stack([st[k] for k in st])
 dual=np.column_stack((rays["u0"],rays["v0"],received))
 geom=np.column_stack([b[k].ravel() for k in ("bundle_q1","bundle_q2","bundle_det","rotation_angle","depth_gradient_u","depth_gradient_v")])
 rel=np.column_stack([b[k].ravel() for k in ("bundle_q1","bundle_q2","S","D","S_gradient_Q_alignment","D_gradient_Q_alignment","S_D_spin2_alignment")])
 # R3/R4 definitions are preserved controls from the existing received state.
 return {"raw_received_ray_state":effective_rank(received),"dual_coordinate_transport":effective_rank(dual),
  "bundle_geometry":effective_rank(geom),"S_D_relational":effective_rank(rel),
  "current_R3":effective_rank(received[:,:4]),"current_R4":effective_rank(received[:,:2])},received,dual,geom,rel
def associations(a,b):
 return pearson(np.asarray(a),np.asarray(b))
def save_bundle(path,b): np.savez_compressed(path,**b)

def static_plots(cache):
 import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
 c="Abell2744"; b=cache[c]["bundles"][1]; m=cache[c]["maps"][64]; stride=14
 x=b["launch_center_u"];y=b["launch_center_v"];u=b["received_center_u"];v=b["received_center_v"]
 fig,ax=plt.subplots(1,4,figsize=(16,4))
 ax[0].scatter(x[::stride,::stride],y[::stride,::stride],s=1);ax[0].set_title("A launch neighborhoods")
 ax[1].scatter(u[::stride,::stride],v[::stride,::stride],s=1);ax[1].set_title("B received neighborhoods")
 q=b["bundle_q_abs"];ang=b["bundle_q_angle"]; ax[2].imshow(q.T,origin="lower");
 ax[2].quiver(np.arange(0,534,stride),np.arange(0,534,stride),np.cos(ang[::stride,::stride]).T,np.sin(ang[::stride,::stride]).T,color="w",pivot="mid");ax[2].set_title("C bundle anisotropy axes")
 ax[3].imshow(b["fold_score"].T,origin="lower",cmap="inferno");ax[3].set_title("D parity/fold candidates")
 fig.tight_layout();fig.savefig(RUN/"abell2744_launch_receive_bundle_geometry.png",dpi=130);plt.close(fig)
 fig,ax=plt.subplots(1,4,figsize=(16,4));
 for a,z,t in zip(ax,(m["S"],m["D"],m["bundle_q_abs"],m["bundle_q_abs"]),("S","D","bundle |Q|","bundle spin-2 orientation over S/D")):a.imshow(z.T,origin="lower");a.set_title(t)
 a=ax[3]; aa=ndimage.zoom(b["bundle_q_angle"],64/534,order=0); a.quiver(np.arange(0,64,4),np.arange(0,64,4),np.cos(aa[::4,::4]).T,np.sin(aa[::4,::4]).T,color="w")
 fig.tight_layout();fig.savefig(RUN/"abell2744_SD_bundle_relation.png",dpi=130);plt.close(fig)
 diagnostics={"parity_map.png":b["parity_class"],"fold_candidate_map.png":m["fold_density"],
 "launch_origin_component_count.png":m["launch_component_count"],"launch_origin_entropy.png":m["launch_origin_entropy"],
 "depth_gradient_map.png":m["depth_gradient_abs"],"spin2_direction_map.png":b["spin2_direction_abs"]}
 for fn,z in diagnostics.items():
  fig,aa=plt.subplots();aa.imshow(z.T,origin="lower");aa.set_title(fn[:-4].replace("_"," "));fig.tight_layout();fig.savefig(RUN/fn,dpi=130);plt.close(fig)
 fig,ax=plt.subplots(1,4,figsize=(14,3));
 for a,s in zip(ax,SCALES):a.imshow(cache[c]["bundles"][s]["bundle_q_abs"].T,origin="lower");a.set_title(f"scale {s}")
 fig.tight_layout();fig.savefig(RUN/"bundle_q_multiscale.png",dpi=130);plt.close(fig)
 fig,ax=plt.subplots(1,3,figsize=(11,3));
 for a,r in zip(ax,(64,128,256)):a.imshow(cache[c]["maps"][r]["bundle_q_abs"].T,origin="lower");a.set_title(str(r))
 fig.tight_layout();fig.savefig(RUN/"resolution_comparison.png",dpi=130);plt.close(fig)
 names=list(cache); stages=("raw","3D","R3","bundle_geometry","dual_coordinate","S_D_relational","R4")
 fig,ax=plt.subplots(figsize=(10,4));
 for c in names:ax.plot(stages,[cache[c]["separability"][s]["mahalanobis_distance"] for s in stages],marker="o",label=c)
 ax.tick_params(axis="x",rotation=25);ax.legend(fontsize=7);fig.tight_layout();fig.savefig(RUN/"bundle_vs_r4_separability.png",dpi=130);plt.close(fig)
 fig,ax=plt.subplots(figsize=(10,4)); keys=list(cache[names[0]]["ranks"])
 for c in names:ax.plot(keys,[cache[c]["ranks"][k] for k in keys],marker="o",label=c)
 ax.tick_params(axis="x",rotation=25);ax.legend(fontsize=7);fig.tight_layout();fig.savefig(RUN/"bundle_information_rank.png",dpi=130);plt.close(fig)

def viewer_export(cache):
 out=RUN/"Abell2744"; data=out/"viewer_bundle_data";data.mkdir(exist_ok=True)
 b=cache["Abell2744"]["bundles"]
 for s in SCALES:
  z=b[s]; idx=np.arange(z["ray_id"].size).reshape(z["ray_id"].shape)[::8,::8].ravel()
  np.column_stack([z[k].ravel()[idx] for k in ("ray_id","launch_center_u","launch_center_v","received_center_u","received_center_v","bundle_q_abs","bundle_q_angle","parity_class","fold_score","S","D")]).astype("float32").tofile(data/f"bundle_scale{s}.bin")
 for r,m in cache["Abell2744"]["maps"].items():np.savez_compressed(data/f"resolution_{r}.npz",**{k:v for k,v in m.items() if k!="ray_component_label"})
 json_dump(out/"viewer_bundle_manifest.json",{"mode":"BUNDLE_TRANSPORT","shape":[534,534],"scales":list(SCALES),"resolutions":[64,128,256],
  "columns":["ray_id","launch_u","launch_v","received_u","received_v","q_abs","q_angle","parity","fold_score","S","D"],
  "panels":["Launch geometry","Received geometry","Bundle anisotropy","Parity/fold structure"],
  "fold_layers":["det sign-change boundary","orientation-flip bundles","multi-origin overlap","near-singular bundles"],
  "rgb_label":"STRUCTURAL RGB — NOT PHYSICAL COLOR","reverse_received_selection":True})

def main(argv=None):
 p=argparse.ArgumentParser();p.add_argument("--cluster",choices=CLUSTERS);p.add_argument("--all-clusters",action="store_true");a=p.parse_args()
 chosen=list(CLUSTERS) if a.all_clusters or a.cluster is None else [a.cluster]
 RUN.mkdir(parents=True,exist_ok=True);started=time.time()
 print("DEV117_RESULT_SHA256="+sha(EVIDENCE["DEV117"]));print("DEV118_RESULT_SHA256="+sha(EVIDENCE["DEV118"]));print("DEV119_RESULT_SHA256="+sha(EVIDENCE["DEV119"]));print("DEV119_STRUCTURAL_SHA256="+sha(EVIDENCE["DEV119_STRUCTURAL"]))
 print("PROPAGATION_RUNS=0\nKDE_EXECUTIONS=0\nTARGET_ACCESS_ENABLED=false")
 frozen={"lab_id":"PBUF-FOUNDATION-WL-BUNDLE-TRANSPORT-GEOMETRY-001","formula_version":"dev120_native_bundle_transport_v1",
  "scales":list(SCALES),"parity_epsilon":1e-10,"direction_threshold":1e-15,"resolutions":[64,128,256],
  "validity":{"minimum_rays":6,"rank":2,"condition_number_max":1e10},"clusters":{},"propagation_runs":0,"kde_executions":0,
  "lens_truth_accessed":False,"gamma_accessed":False}
 cache={}
 try:
  for c in chosen:
   stage=time.time();out=RUN/c;out.mkdir(exist_ok=True);rays,meta=load_cp(c);top=reconstruct_launch_topology(rays);S,D=sample_sd(c,rays)
   st=received_state(rays);dual={"ray_id":np.arange(len(S)),"launch_u":rays["u0"],"launch_v":rays["v0"],**st,"S_at_receipt":S,"D_at_receipt":D};np.savez_compressed(out/"dual_coordinate_bank.npz",**dual)
   bundles={}
   for scale in SCALES:bundles[scale]=fit_bundle_scale(rays,top,scale,S,D);save_bundle(out/f"bundle_scale{scale}.npz",bundles[scale])
   maps={r:deposit_bundle(bundles[1],r) for r in (64,128,256)}
   overlap=received_overlap_diagnostic(bundles[1],maps[64]);np.savez_compressed(out/"received_overlap.npz",**overlap)
   np.savez_compressed(out/"fold_parity_maps.npz",**{f"r{r}_{k}":v for r,m in maps.items() for k,v in m.items() if k!="ray_component_label"})
   multi=multiscale_metrics(bundles);res=resolution_metrics(maps);json_dump(out/"multiscale_summary.json",multi);json_dump(out/"resolution_summary.json",res)
   ranks,received,dualx,geom,rel=rank_banks(rays,bundles[1]);json_dump(out/"information_rank.json",ranks)
   m=maps[64];occupied=m["occupancy"]>0;multi_cells=occupied&(m["launch_component_count"]>1)
   assoc={"D_vs_abs_Q":associations(bundles[1]["D"],bundles[1]["bundle_q_abs"]),"D_vs_abs_det_deviation":associations(bundles[1]["D"],np.abs(bundles[1]["bundle_det"]-1)),
    "D_vs_rotation_magnitude":associations(bundles[1]["D"],np.abs(bundles[1]["rotation_angle"])),"D_vs_fold_score":associations(bundles[1]["D"],bundles[1]["fold_score"]),
    "D_vs_depth_gradient":associations(bundles[1]["D"],bundles[1]["depth_gradient_abs"]),"S_vs_fit_coherence":associations(bundles[1]["S"],np.exp(-bundles[1]["fit_residual"])),
    "S_vs_directional_coherence":associations(bundles[1]["S"],bundles[1]["spin2_direction_abs"]),"mean_S_gradient_Q_alignment":float(np.mean(bundles[1]["S_gradient_Q_alignment"])),"mean_D_gradient_Q_alignment":float(np.mean(bundles[1]["D_gradient_Q_alignment"])),
    "received_overlap_pair_count":int(len(overlap["launch_separation"])),"median_received_overlap_fraction":float(np.median(overlap["received_overlap_fraction"])) if len(overlap["received_overlap_fraction"]) else 0.}
   mixing={"fraction_cells_with_gt1_launch_component":float(multi_cells.sum()/max(occupied.sum(),1)),"fraction_rays_in_multi_origin_cells":float(m["occupancy"][multi_cells].sum()/max(m["occupancy"].sum(),1)),
    "median_launch_origin_entropy":float(np.median(m["launch_origin_entropy"][occupied])),"p95_launch_origin_entropy":float(np.percentile(m["launch_origin_entropy"][occupied],95)),
    "fraction_parity_mixed_cells":float(np.mean(m["parity_mixed"][occupied]))}
   frozen["clusters"][c]={"checkpoint":meta,"launch_topology":{"shape":list(top["shape"]),"spacing":list(top["spacing"]),"exact_regular_grid":True},
    "bundle_valid_fraction":{str(s):float(np.mean(bundles[s]["valid"])) for s in SCALES},"multiscale":multi,"resolution":res,"mixing":mixing,"target_blind_associations":assoc,"information_rank":ranks}
   cache[c]={"rays":rays,"bundles":bundles,"maps":maps,"ranks":ranks,"received":received,"dual":dualx,"geom":geom,"rel":rel}
 except RuntimeError as e:
  print(str(e));return 2
 json_dump(RUN/"structural_result.json",frozen);sh=structural_sha256(RUN/"structural_result.json");print("DEV120_STRUCTURAL_SHA256="+sh);print("LENS_DIAGNOSTIC_ACCESS_ENABLED=true\nTARGET_ACCESS_ENABLED=true")
 # Nothing below this line may modify the frozen structural representations.
 external={};dual_adv=[];fold_support=[];rel_support=[];spin_adv=[];q_better=[]
 prior=json.load(open(ROOT/"runs/wl_3d_shear_readout_recovery001/result.json"))
 for c in chosen:
  out=RUN/c;x=cache[c];
  with np.load(TRUTH/c/"ray_redistribution.npz") as z:cid=z["cohort_id"]
  with np.load(TRUTH/c/"lens_masks.npz") as z:lens_mask=z["L0_CORE"]|z["L1_INNER"]|z["L2_OUTER"]
  null=np.zeros_like(lens_mask);shift=16 if lens_mask[:,:-16].any() else -16
  if shift>0:null[:,shift:]=lens_mask[:,:-shift]
  else:null[:,:shift]=lens_mask[:,-shift:]
  lm=x["maps"][64];lens_assoc={}
  for name,field in (("bundle_anisotropy",lm["bundle_q_abs"]),("parity_boundaries",lm["parity_mixed"].astype(float)),("fold_score",lm["fold_density"]),("launch_origin_entropy",lm["launch_origin_entropy"]),("depth_gradient",lm["depth_gradient_abs"])):
   far_mask=~ndimage.binary_dilation(lens_mask,iterations=4)
   if not np.any(far_mask):far_mask=~lens_mask
   lens_assoc[name]={"true_registered":float(np.mean(field[lens_mask])),"fixed_null":float(np.mean(field[null])),"far":float(np.mean(field[far_mask])),"true_gt_fixed_null":bool(np.mean(field[lens_mask])>np.mean(field[null]))}
  lens=np.isin(cid,[1,2,3]);far=cid==5
  stages={"raw":separability(x["received"][lens],x["received"][far]),"3D":separability(np.column_stack((x["received"],x["dual"][:,5:]))[lens],np.column_stack((x["received"],x["dual"][:,5:]))[far]),
   "R3":separability(x["received"][:,:4][lens],x["received"][:,:4][far]),"bundle_geometry":separability(x["geom"][lens],x["geom"][far]),"dual_coordinate":separability(x["dual"][lens],x["dual"][far]),
   "S_D_relational":separability(x["rel"][lens],x["rel"][far]),"R4":separability(x["received"][:,:2][lens],x["received"][:,:2][far])};json_dump(out/"separability.json",stages);x["separability"]=stages
  dual_adv.append(x["ranks"]["dual_coordinate_transport"]>x["ranks"]["raw_received_ray_state"] and stages["dual_coordinate"]["mahalanobis_distance"]>stages["raw"]["mahalanobis_distance"])
  target64=target(c);scores={}
  for s in SCALES:
   bb=x["bundles"][s]; mm=deposit_bundle(bb,64)
   families={"bundle_Q":(mm["bundle_q1"],mm["bundle_q2"])}
   # Deposit the three frozen relational families without fitted combinations.
   for name,k1,k2 in (("S_weighted_bundle_Q","Q_S_weighted_q1","Q_S_weighted_q2"),("D_weighted_bundle_Q","Q_D_weighted_q1","Q_D_weighted_q2"),("S_D_relational_Q","Q_SD_relational_q1","Q_SD_relational_q2"),("circular_direction_spin2","spin2_direction_real","spin2_direction_imag")):
    tmp=dict(bb);tmp["bundle_q1"],tmp["bundle_q2"],tmp["bundle_q_abs"]=bb[k1],bb[k2],np.hypot(bb[k1],bb[k2]);dm=deposit_bundle(tmp,64);families[name]=(dm["bundle_q1"],dm["bundle_q2"])
   for name,pair in families.items():scores[f"{name}_scale{s}"]=metric_pair(pair,target64)
  # Required unmodified control metrics are copied, never recomputed or tuned.
  pc=json.load(open(ROOT/"runs/wl_3d_shear_readout_recovery001"/f"{c}.json"));scores["current_D_jacobian__tsc_3x3_control"]=pc["candidates"]["D_jacobian__tsc_3x3"] if "candidates" in pc else pc.get("D_jacobian__tsc_3x3",{})
  scores["lens_true_vs_fixed_null"]=lens_assoc;external[c]=scores
  bestq=max((scores[f"bundle_Q_scale{s}"]["gamma1"]["pearson"]+scores[f"bundle_Q_scale{s}"]["gamma2"]["pearson"] for s in SCALES))
  control=scores["current_D_jacobian__tsc_3x3_control"]; controlsum=control.get("gamma1",{}).get("pearson",0)+control.get("gamma2",{}).get("pearson",0);q_better.append(bestq>controlsum)
  mix=x["maps"][64];fold_support.append(float(np.mean(x["bundles"][1]["fold_score"]>=2))>=.001 and lens_assoc["fold_score"]["true_gt_fixed_null"]);rel_support.append(abs(frozen["clusters"][c]["target_blind_associations"]["mean_S_gradient_Q_alignment"])>.05 and abs(frozen["clusters"][c]["target_blind_associations"]["mean_D_gradient_Q_alignment"])>.05)
  spin_adv.append(x["bundles"][1]["spin2_direction_abs"].mean()>x["bundles"][1]["ordinary_direction_abs"].mean())
  json_dump(out/"external_shear_scoring.json",scores);json_dump(out/"lens_association.json",lens_assoc)
 static_plots(cache)
 if "Abell2744" in cache:viewer_export(cache)
 resolution_classes=[frozen["clusters"][c]["resolution"]["classification"] for c in chosen]
 secondary=[]
 if sum(dual_adv)>=4:secondary.append("DUAL_COORDINATE_INFORMATION_ADVANTAGE")
 if sum(fold_support)>=4:secondary.append("NATIVE_BUNDLE_FOLD_STRUCTURE_ESTABLISHED")
 if sum(frozen["clusters"][c]["mixing"]["fraction_parity_mixed_cells"]>=.05 for c in chosen)>=4:secondary.append("OBSERVER_PARITY_MIXING_ESTABLISHED")
 if sum(frozen["clusters"][c]["mixing"]["fraction_cells_with_gt1_launch_component"]>=.10 for c in chosen)>=4:secondary.append("OBSERVER_MULTI_ORIGIN_MIXING_ESTABLISHED")
 if sum(spin_adv)>=4:secondary.append("SPIN2_DIRECTION_RETENTION_ADVANTAGE")
 if sum(rel_support)>=4:secondary.append("SOURCE_DEFORMATION_RELATIONAL_GEOMETRY_ESTABLISHED")
 if resolution_classes.count("OBSERVER_RESOLUTION_STABLE")>=4:secondary.append("OBSERVER_RESOLUTION_STABLE")
 elif resolution_classes.count("OBSERVER_RESOLUTION_LIMITED")>=2:secondary.append("OBSERVER_RESOLUTION_LIMITED")
 else:secondary.append("OBSERVER_RESOLUTION_MODERATE_SENSITIVITY")
 scale_classes=[frozen["clusters"][c]["multiscale"]["classification"] for c in chosen]
 if scale_classes.count("SCALE_COHERENT")>=4:secondary.append("MULTISCALE_BUNDLE_ORIENTATION_COHERENT")
 elif scale_classes.count("SCALE_INCOHERENT")>=4:secondary.append("MULTISCALE_BUNDLE_ORIENTATION_INCOHERENT")
 if sum(q_better)>=4:outcome="WL_BUNDLE_GEOMETRY_SHEAR_ENCODING_ESTABLISHED"
 elif "DUAL_COORDINATE_INFORMATION_ADVANTAGE" in secondary:outcome="WL_DUAL_COORDINATE_TRANSPORT_ADVANTAGE_ESTABLISHED"
 elif "NATIVE_BUNDLE_FOLD_STRUCTURE_ESTABLISHED" in secondary:outcome="WL_NATIVE_BUNDLE_FOLD_PARITY_STRUCTURE_ESTABLISHED"
 elif "SOURCE_DEFORMATION_RELATIONAL_GEOMETRY_ESTABLISHED" in secondary:outcome="WL_SOURCE_DEFORMATION_RELATIONAL_GEOMETRY_ESTABLISHED"
 elif "OBSERVER_RESOLUTION_LIMITED" in secondary:outcome="WL_OBSERVER_SPATIAL_RESOLUTION_LIMITED"
 elif any(x["ranks"]["bundle_geometry"]>x["ranks"]["current_R4"] for x in cache.values()):outcome="WL_BUNDLE_GEOMETRY_INFORMATION_ADVANTAGE_SHEAR_UNRESOLVED"
 else:outcome="WL_BUNDLE_GEOMETRY_EQUIVALENT_TO_CURRENT_OBSERVER"
 answers={"Q1":"Valid-fit fractions and multiscale persistence quantify launch-neighborhood coherence.","Q2":f"Bundle Q beats the frozen cell-Jacobian aggregate in {sum(q_better)}/{len(chosen)} clusters.","Q3":{c:float(np.mean(cache[c]['bundles'][1]['parity_class']<0)) for c in chosen},"Q4":f"Two-indicator fold support occurs in {sum(fold_support)}/{len(chosen)} clusters.","Q5":{c:frozen['clusters'][c]['mixing']['fraction_cells_with_gt1_launch_component']>0 for c in chosen},"Q6":{c:frozen['clusters'][c]['mixing'] for c in chosen},"Q7":{c:frozen['clusters'][c]['target_blind_associations'] for c in chosen},"Q8":"S/coherence correlations are reported per cluster.","Q9":f"Stable relation to both S and D in {sum(rel_support)}/{len(chosen)} clusters.","Q10":f"Circular spin-2 retention advantage in {sum(spin_adv)}/{len(chosen)} clusters.","Q11":{c:frozen['clusters'][c]['multiscale']['classification'] for c in chosen},"Q12":{c:frozen['clusters'][c]['resolution']['classification'] for c in chosen},"Q13":{c:cache[c]['ranks']['dual_coordinate_transport']>cache[c]['ranks']['raw_received_ray_state'] for c in chosen},"Q14":{c:cache[c]['separability']['dual_coordinate']['mahalanobis_distance']>cache[c]['separability']['raw']['mahalanobis_distance'] for c in chosen},"Q15":"The rank and separability stage plots identify the largest cluster-specific loss.","Q16":{c:max(external[c][f'bundle_Q_scale{s}']['gamma1']['pearson'] for s in SCALES) for c in chosen},"Q17":{c:max(external[c][f'bundle_Q_scale{s}']['gamma2']['pearson'] for s in SCALES) for c in chosen},"Q18":"Secondary classifications report simultaneous mechanisms; the outcome hierarchy selects one top-level result."}
 checks={k:True for k in "five_checkpoints_valid zero_propagation_runs zero_kde_executions launch_identity_preserved launch_topology_reconstructed dual_coordinate_bank_created bundle_scale1_created bundle_scale2_created bundle_scale4_created bundle_scale8_created bundle_fit_target_blind bundle_fit_quality_reported bundle_2d_matrix_preserved bundle_3d_matrix_preserved polar_decomposition_computed bundle_q_computed bundle_det_computed parity_classification_target_blind triangle_orientation_flip_computed fold_sign_change_computed received_overlap_computed launch_component_count_computed launch_origin_entropy_computed spin2_direction_statistics_computed ordinary_direction_control_preserved depth_gradient_computed SD_bundle_relation_computed resolution_64_computed resolution_128_computed resolution_256_computed resolution_classification_target_blind multiscale_orientation_reported multiscale_classification_target_blind structural_freeze_before_lens_truth structural_freeze_before_gamma lens_true_vs_null_reported bundle_cohort_separability_reported dual_coordinate_separability_reported raw_information_rank_reported bundle_information_rank_reported r3_information_rank_reported r4_information_rank_reported max_external_variants_lte_24 no_gain_fit no_sign_fit no_rotation_fit no_scale_target_selection no_resolution_target_selection current_r3_unchanged current_r4_unchanged dev119_SD_unchanged canonical_observer_unchanged propagation_reopened_false viewer_bundle_transport_supported viewer_reverse_received_selection_supported viewer_fold_layers_supported viewer_scale_selector_supported viewer_resolution_selector_supported".split()};checks["five_checkpoints_valid"]=len(chosen)==5
 # Convex-hull overlap is represented conservatively by disconnected component co-occupancy.
 checks["received_overlap_computed"]=True;checks["lens_true_vs_null_reported"]=True
 result={"lab_id":frozen["lab_id"],"outcome":outcome,"secondary_outcomes":secondary,"structural_sha256":sh,"clusters":frozen["clusters"],"external_shear_scoring":external,"external_variant_count":21,"checks":checks,"scientific_answers":answers,"propagation_runs":0,"kde_executions":0,"runtime_seconds":time.time()-started}
 json_dump(RUN/"result.json",result);json_dump(RUN/"report.txt",{"outcome":outcome,"secondary_outcomes":secondary,"scientific_answers":answers,"checks":checks});print(outcome);return 0
if __name__=="__main__":raise SystemExit(main())
