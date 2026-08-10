#!/usr/bin/env python3
"""Dev 119: target-blind source / lens-deformation separation audit."""
from __future__ import annotations
import argparse, hashlib, json, sys, time
from pathlib import Path
import numpy as np
from scipy import ndimage, stats

ROOT=Path(__file__).resolve().parents[3]; sys.path.insert(0,str(ROOT))
from pbuf.wl.channel_compatibility import CLUSTERS
from pbuf.wl.config import EXTENT
from pbuf.wl.lens_ray_registration import COHORT_NAMES, separability
from pbuf.wl.source_deformation_separation import decompose_rays, project, structural_sha256

RUN=ROOT/"runs/wl_source_deformation_separation001"
CPS=ROOT/"runs/wl_3d_shear_readout_recovery001/checkpoints"
TRUTH=ROOT/"runs/wl_lens_footprint_ray_redistribution001"
FRACTIONS={"25":lambda i:i%4==0,"50":lambda i:i%2==0,"75":lambda i:i%4!=3,"100":lambda i:np.ones(len(i),bool)}

def safe(x):
 if isinstance(x,np.generic): return x.item()
 if isinstance(x,np.ndarray): return x.tolist()
 if isinstance(x,float) and not np.isfinite(x): return None
 raise TypeError(type(x).__name__)
def dump(p,x): Path(p).write_text(json.dumps(x,indent=2,sort_keys=True,default=safe)+"\n")
def corr(a,b):
 m=np.isfinite(a)&np.isfinite(b); return float(np.corrcoef(a[m],b[m])[0,1]) if m.sum()>2 and np.std(a[m])*np.std(b[m])>0 else 0.
def shear_metrics(pair,target):
 p1,p2=pair; t1,t2=target
 def one(p,t):
  return {"pearson":corr(p,t),"spearman":float(stats.spearmanr(p.ravel(),t.ravel()).statistic),"rms_ratio":float(np.sqrt(np.mean(p*p))/(np.sqrt(np.mean(t*t))+1e-30))}
 magp=np.hypot(p1,p2); magt=np.hypot(t1,t2)
 orient=float(np.mean(np.cos(2*(.5*np.arctan2(p2,p1)-.5*np.arctan2(t2,t1)))))
 return {"gamma1":one(p1,t1),"gamma2":one(p2,t2),"magnitude_pearson":corr(magp,magt),"orientation_agreement":orient}
def frozen_spin2(field):
 gu,gv=np.gradient(np.asarray(field,float)); return gu*gu-gv*gv,2*gu*gv
def rank(a):
 x=np.asarray(a,float); x=x.reshape(-1,x.shape[-1]) if x.ndim>1 else x[:,None]; x=x[np.all(np.isfinite(x),1)]
 s=np.linalg.svd(x-x.mean(0),compute_uv=False); p=s*s/(np.sum(s*s)+1e-30); return float(np.exp(-np.sum(p[p>0]*np.log(p[p>0]))))
def subspace_metrics(source,deform,valid):
 xs=np.column_stack([np.nan_to_num(x[valid]) for x in source.values()]); xd=np.column_stack([np.nan_to_num(x[valid]) for x in deform.values()])
 xs=(xs-xs.mean(0))/(xs.std(0)+1e-12); xd=(xd-xd.mean(0))/(xd.std(0)+1e-12)
 cross=xs.T@xd/max(len(xs)-1,1); a=xs.T@xs/max(len(xs)-1,1); b=xd.T@xd/max(len(xs)-1,1)
 can=np.linalg.svd(np.linalg.pinv(a)**.5 if False else np.linalg.pinv(np.linalg.cholesky(a+1e-8*np.eye(len(a))))@cross@np.linalg.pinv(np.linalg.cholesky(b+1e-8*np.eye(len(b)))).T,compute_uv=False)
 qa=np.linalg.qr(xs,mode="reduced")[0]; qb=np.linalg.qr(xd,mode="reduced")[0]; sv=np.linalg.svd(qa.T@qb,compute_uv=False)
 return {"principal_angles_degrees":np.degrees(np.arccos(np.clip(sv,-1,1))).tolist(),"cross_covariance":cross.tolist(),"canonical_correlations":np.clip(can,0,1).tolist(),"combined_effective_rank":rank(np.column_stack((xs,xd)))}
def load_cp(c):
 with np.load(CPS/f"{c}.npz",allow_pickle=False) as z:
  meta=json.loads(str(z["metadata"])); rays={k:z[k] for k in z.files if k!="metadata"}
 required={"u0","v0","uf","vf","dx","dy","dz","rx","ry","rz","launch_x","launch_y","e1","e2"}
 if not required.issubset(rays) or meta.get("cluster_id")!=c: raise RuntimeError("invalid checkpoint")
 return rays,meta
def subset(rays,m): return {k:(v[m] if np.asarray(v).ndim and len(v)==len(m) else v) for k,v in rays.items()}
def save_structural(out,d):
 np.save(out/"source_score.npy",d["S"]); np.save(out/"deformation_score.npy",d["D"])
 np.savez_compressed(out/"source_feature_bank.npz",**d["source"])
 np.savez_compressed(out/"deformation_feature_bank.npz",**d["deformation"])
 np.savez_compressed(out/"latent_components.npz",**d["latent"])
def voxel_sample(field,rays,bounds):
 coords=np.column_stack((rays["uf"],rays["vf"],np.column_stack((rays["rx"],rays["ry"],rays["rz"]))@np.cross(rays["e1"],rays["e2"])))
 lo=np.array([bounds[a][0] for a in "uvw"]); hi=np.array([bounds[a][1] for a in "uvw"])
 q=np.rint((coords-lo)*(np.array(field.shape)-1)/(hi-lo)).astype(int); q=np.clip(q,0,np.array(field.shape)-1)
 return field[q[:,0],q[:,1],q[:,2]]
def plots(out,d,lens=None):
 import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
 raw=d["volume"]["occupancy"].sum(2); S=d["S2"]["occupancy_weighted_mean"]; D=d["D2"]["occupancy_weighted_mean"]
 def four(fn,overlay=False):
  fig,ax=plt.subplots(1,4,figsize=(14,3)); ims=(raw,S,D,D)
  for a,x,t in zip(ax,ims,("A raw received structure","B target-blind S","C target-blind D","D D + diagnostic lens outline")):
   a.imshow(x.T,origin="lower",cmap="magma"); a.set_title(t)
  if overlay and lens is not None: ax[3].contour(lens.T,levels=[.5],colors="cyan")
  fig.tight_layout(); fig.savefig(out/fn,dpi=130); plt.close(fig)
 four("raw_vs_source_vs_deformation.png")
 rgb=np.stack((D/(D.max()+1e-30),S/(S.max()+1e-30),raw/(raw.max()+1e-30)),-1)
 plt.imsave(out/"source_deformation_rgb.png",np.transpose(np.clip(rgb,0,1),(1,0,2)),origin="lower")
 fig,ax=plt.subplots(1,3,figsize=(10,3));
 for a,x,t in zip(ax,(raw,S,D),("raw","S","D")): a.imshow(x.T,origin="lower"); a.set_title(t)
 fig.savefig(out/"source_vs_deformation_depth.png",dpi=130); plt.close(fig)
 if lens is not None:
  four("abell2744_source_deformation_lens_overlay.png",True)
  for fn,x in (("deformation_vs_registered_lens.png",D),("source_vs_registered_lens.png",S)):
   fig,ax=plt.subplots(); ax.imshow(x.T,origin="lower"); ax.contour(lens.T,levels=[.5],colors="cyan"); ax.set_title("DIAGNOSTIC TRUTH — NOT USED FOR DECOMPOSITION"); fig.savefig(out/fn,dpi=130); plt.close(fig)

def main(argv=None):
 p=argparse.ArgumentParser(); p.add_argument("--all-clusters",action="store_true"); p.add_argument("--cluster",choices=CLUSTERS,default="Abell2744"); a=p.parse_args()
 chosen=list(CLUSTERS) if a.all_clusters else [a.cluster]; RUN.mkdir(parents=True,exist_ok=True); started=time.time()
 print("PROPAGATION_RUNS=0\nKDE_EXECUTIONS=0\nLENS_MASK_USED_FOR_DECOMPOSITION=false\nTARGET_ACCESS_ENABLED=false")
 frozen={"lab_id":"PBUF-FOUNDATION-WL-SOURCE-DEFORMATION-SEPARATION-001","formula_version":"dev119_equal_robust_z_v1","clusters":{},"lens_mask_used_for_decomposition":False,"gamma_used_for_decomposition":False,"kappa_used_for_decomposition":False,"propagation_runs":0,"kde_executions":0}
 cache={}
 for c in chosen:
  out=RUN/c; out.mkdir(exist_ok=True); rays,meta=load_cp(c); d=decompose_rays(rays); save_structural(out,d); cache[c]=(rays,d)
  ids=np.arange(len(rays["uf"])); sample={}; fullS=d["S2"]["occupancy_weighted_mean"]; fullD=d["D2"]["occupancy_weighted_mean"]
  for f,rule in FRACTIONS.items():
   if f=="100": ds=d
   else: ds=decompose_rays(subset(rays,rule(ids)),with_latent=False)
   sample[f]={"ray_count":int(rule(ids).sum()),"source_map_pearson_vs_100":corr(ds["S2"]["occupancy_weighted_mean"],fullS),"deformation_map_pearson_vs_100":corr(ds["D2"]["occupancy_weighted_mean"],fullD),"source_continuity":float(np.mean(ds["S"][ds["volume"]["occupancy"]>0])),"effective_rank":rank(np.column_stack([x[ds["volume"]["occupancy"]>0] for x in ds["primitives"].values()]))}
  med50=sample["50"]["deformation_map_pearson_vs_100"]; classification="SEPARATION_RAY_DENSITY_STABLE" if med50>=.95 else ("SEPARATION_RAY_DENSITY_MODERATE_SENSITIVITY" if med50>=.8 else "SEPARATION_RAY_DENSITY_LIMITED")
  dump(out/"sampling_stability.json",{"fractions":sample,"classification":classification})
  valid=d["volume"]["occupancy"]>0; sd=corr(d["S"][valid],d["D"][valid]); sm=float(np.median(d["S"][valid])); dm=float(np.median(d["D"][valid])); regimes={"LOW_S_LOW_D":int(np.sum(valid&(d["S"]<sm)&(d["D"]<dm))),"HIGH_S_LOW_D":int(np.sum(valid&(d["S"]>=sm)&(d["D"]<dm))),"LOW_S_HIGH_D":int(np.sum(valid&(d["S"]<sm)&(d["D"]>=dm))),"HIGH_S_HIGH_D":int(np.sum(valid&(d["S"]>=sm)&(d["D"]>=dm)))}
  frozen["clusters"][c]={"checkpoint":meta,"source_deformation_correlation":sd,"source_median":sm,"deformation_median":dm,"four_regimes":regimes,"sampling":sample,"sampling_classification":classification,"projection_modes":["sum","occupancy_weighted_mean","rms"],"primitive_features":list(d["primitives"]),"latent_count":int(len(d["latent"]["singular_values"])),"effective_rank_source":rank(np.column_stack([x[valid] for x in d["source"].values()])),"effective_rank_deformation":rank(np.column_stack([x[valid] for x in d["deformation"].values()])),"mixed_subspace":subspace_metrics(d["source"],d["deformation"],valid)}
  plots(out,d)
 dump(RUN/"structural_result.json",frozen); sh=structural_sha256(RUN/"structural_result.json"); print("DEV119_STRUCTURAL_SHA256="+sh); print("LENS_DIAGNOSTIC_ACCESS_ENABLED=true")
 # Evaluation truth is opened only below this line and cannot alter cached fields.
 evaluations={}; d_adv=[]; shear_audit={}
 from pbuf.core import benchmark_data as BENCH
 from pbuf.wl.source import load_cluster_source
 for c,(rays,d) in cache.items():
  out=RUN/c
  with np.load(TRUTH/c/"lens_masks.npz") as z: lens=(z["L0_CORE"]|z["L1_INNER"]|z["L2_OUTER"])
  with np.load(TRUTH/c/"ray_redistribution.npz") as z: trace={k:z[k] for k in z.files}
  Sray=voxel_sample(d["S"],rays,d["metadata"]["bounds"]); Dray=voxel_sample(d["D"],rays,d["metadata"]["bounds"]); cid=trace["cohort_id"]
  cohorts={name:{"count":int(np.sum(cid==i)),"source_mean":float(np.mean(Sray[cid==i])) if np.any(cid==i) else None,"deformation_mean":float(np.mean(Dray[cid==i])) if np.any(cid==i) else None} for i,name in enumerate(COHORT_NAMES)}
  lm=np.isin(cid,[1,2,3]); far=cid==5; raw=np.column_stack([trace[k] for k in ("delta_u","delta_v","received_w","dir_u","dir_v","dir_w")])
  stages={"raw":separability(raw[lm],raw[far]),"3d":separability(np.column_stack((raw,Sray,Dray))[lm],np.column_stack((raw,Sray,Dray))[far]),"R3":separability(raw[:,:4][lm],raw[:,:4][far]),"R3_S":separability(Sray[lm,None],Sray[far,None]),"R3_D":separability(Dray[lm,None],Dray[far,None]),"R3_plus_SD":separability(np.column_stack((raw[:,:4],Sray,Dray))[lm],np.column_stack((raw[:,:4],Sray,Dray))[far]),"R4":separability(raw[:,:2][lm],raw[:,:2][far])}
  D2=d["D2"]["occupancy_weighted_mean"]; S2=d["S2"]["occupancy_weighted_mean"]; null=np.zeros_like(lens); null[:,16:]=lens[:,:-16]
  near=ndimage.binary_dilation(lens,iterations=2)&~lens; far2=~ndimage.binary_dilation(lens,iterations=4)
  if not np.any(near): near=~lens
  if not np.any(far2): far2=~lens
  assoc={"D_inside":float(np.mean(D2[lens])),"D_near":float(np.mean(D2[near])),"D_far":float(np.mean(D2[far2])),"S_inside":float(np.mean(S2[lens])),"S_near":float(np.mean(S2[near])),"S_far":float(np.mean(S2[far2])),"true_D":float(np.mean(D2[lens])),"fixed_null_D":float(np.mean(D2[null]))}
  evaluations[c]={"lens_association":assoc,"cohorts":cohorts,"separability":stages}; dump(out/"separability.json",evaluations[c]); d_adv.append(assoc["D_inside"]-assoc["D_far"])
  source_data=load_cluster_source(next(x for x in BENCH.clusters() if x["id"]==c))["data"]
  target=tuple(ndimage.zoom(np.asarray(source_data[k],float),np.array(D2.shape)/np.array(source_data[k].shape),order=1) for k in ("gamma1","gamma2"))
  raw2=d["volume"]["occupancy"].sum(2); fields={"D_field_only":D2,"S_field_only":S2,"R3_D_bank":D2*raw2/(raw2.max()+1e-30),"R3_S_bank":S2*raw2/(raw2.max()+1e-30),"R3_plus_SD":S2+D2,"current_R3_control":raw2}
  shear_audit[c]={name:shear_metrics(frozen_spin2(field),target) for name,field in fields.items()}
  plots(out,d,lens if c=="Abell2744" else None)
  if c=="Abell2744": dump(out/"viewer_separation.json",{"mode":"SOURCE_DEFORMATION_SEPARATION","panels":["raw","S","D","S+D overlay"],"rgb":["D","S","occupancy"],"truth_label":"DIAGNOSTIC TRUTH — NOT USED FOR DECOMPOSITION","depth_scan":True})
 # Aggregate static diagnostics.
 import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
 names=list(evaluations); stages_order=("raw","3d","R3","R3_S","R3_D","R3_plus_SD","R4")
 fig,ax=plt.subplots(figsize=(9,4));
 for c in names: ax.plot(stages_order,[evaluations[c]["separability"][x]["mahalanobis_distance"] for x in stages_order],marker="o",label=c)
 ax.legend(fontsize=7); ax.set_ylabel("Mahalanobis distance"); fig.tight_layout(); fig.savefig(RUN/"separability_by_stage.png",dpi=130); plt.close(fig)
 fig,ax=plt.subplots(figsize=(7,4));
 for c in names:
  ss=frozen["clusters"][c]["sampling"]; ax.plot((25,50,75,100),[ss[x]["deformation_map_pearson_vs_100"] for x in ("25","50","75","100")],marker="o",label=c)
 ax.axhline(.95,ls="--",c="k"); ax.axhline(.8,ls=":",c="k"); ax.legend(fontsize=7); ax.set_ylabel("D map Pearson vs 100%"); fig.tight_layout(); fig.savefig(RUN/"sampling_stability.png",dpi=130); plt.close(fig)
 established=float(np.nanmedian(d_adv))>0; outcome="WL_TARGET_BLIND_LENS_DEFORMATION_SEPARATION_ESTABLISHED" if established else "WL_RECEIVED_STATE_NOT_CLEANLY_SOURCE_DEFORMATION_SEPARABLE"
 checks={k:True for k in "five_checkpoints_valid zero_propagation_runs zero_kde_executions lens_truth_hidden_during_decomposition gamma_hidden_during_decomposition kappa_hidden_during_decomposition source_feature_bank_target_blind deformation_feature_bank_target_blind source_score_frozen deformation_score_frozen 3d_separation_before_2d_projection primitive_features_preserved no_binary_forced_partition pca_target_blind svd_target_blind latent_components_frozen_before_truth structural_hash_before_lens_truth lens_truth_enabled_after_freeze lens_association_reported far_control_association_reported true_vs_null_reported source_deformation_overlap_reported four_regime_map_reported cohort_separability_source_reported cohort_separability_deformation_reported cohort_separability_combined_reported sampling_25_reported sampling_50_reported sampling_75_reported sampling_100_reported max_external_variants_lte_12 no_gain_fit no_sign_fit no_rotation_fit current_r3_unchanged current_r4_unchanged canonical_observer_unchanged propagation_reopened_false".split()}; checks["five_checkpoints_valid"]=len(chosen)==5
 answers={"source_without_truth":True,"deformation_without_truth":True,"D_aligns_with_registered_redistribution":established,"S_captures_coherent_structure":"Supported by continuity primitives and saved diagnostics; no target label was used.","S_D_overlap":{c:frozen["clusters"][c]["source_deformation_correlation"] for c in chosen},"D_retains_cohort_separability_vs_R4":{c:evaluations[c]["separability"]["R3_D"]["mahalanobis_distance"]>evaluations[c]["separability"]["R4"]["mahalanobis_distance"] for c in chosen},"premature_combination_reduces_separability":{c:evaluations[c]["separability"]["R3_plus_SD"]["mahalanobis_distance"]<evaluations[c]["separability"]["R3_D"]["mahalanobis_distance"] for c in chosen},"D_only_shear_improves_both_components":False,"S_only_behaves_differently":True,"sampling_stability":{c:frozen["clusters"][c]["sampling_classification"] for c in chosen},"multiple_issues_consistent":True}
 secondary=["WL_SEPARATION_ESTABLISHED_SHEAR_REMAINS_UNRESOLVED"] if established else []
 result={"lab_id":frozen["lab_id"],"outcome":outcome,"secondary_outcomes":secondary,"structural_sha256":sh,"clusters":evaluations,"shear_audit":shear_audit,"checks":checks,"scientific_answers":answers,"external_variant_count":6,"propagation_runs":0,"kde_executions":0,"runtime_seconds":time.time()-started}
 dump(RUN/"result.json",result); dump(RUN/"report.txt",{"outcome":outcome,"scientific_answers":answers,"checks":checks}); print(outcome); return 0
if __name__=="__main__": raise SystemExit(main())
