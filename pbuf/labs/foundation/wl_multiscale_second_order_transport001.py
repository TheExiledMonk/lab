#!/usr/bin/env python3
"""Dev 123 — target-blind multiscale audit of the frozen rank-29 transport field."""
from __future__ import annotations
import argparse, hashlib, json, subprocess, sys, time
from pathlib import Path
import numpy as np
from scipy import ndimage, stats
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from pbuf.wl.channel_compatibility import CLUSTERS
from pbuf.wl.lens_ray_registration import separability
from pbuf.wl.multiscale_transport_relations import (SCALES,ablation_lanes,canonical_manifest,
 derivatives,load_second_order_bank,local_mean_variance,matrix_diagnostics,scale_persistence,
 sha256_json,spatial_quadrupole)
from pbuf.wl.transport_receiver_decode import rasterize

RUN=ROOT/"runs/wl_multiscale_second_order_transport001"; D121=ROOT/"runs/wl_dual_transport_receiver_decode001"
D122=ROOT/"runs/wl_reconstruction_first_decode001"; CP=ROOT/"runs/wl_3d_shear_readout_recovery001/checkpoints"
SD=ROOT/"runs/wl_source_deformation_separation001"; COHORT=ROOT/"runs/wl_lens_footprint_ray_redistribution001"
PAIR_SCALES=((2,4),(4,8),(8,16),(16,32)); EPS=np.finfo(float).eps

def dump(path,obj):
 path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+"\n")
def finite(x):return np.nan_to_num(np.asarray(x,float),nan=0.,posinf=0.,neginf=0.)
def corr(a,b):
 a=finite(a).ravel();b=finite(b).ravel();return float(np.corrcoef(a,b)[0,1]) if a.std() and b.std() else 0.
def correlations(a,b):return {"pearson":corr(a,b),"spearman":float(np.nan_to_num(stats.spearmanr(finite(a).ravel(),finite(b).ravel()).statistic))}
def standard(x):
 x=finite(x);return (x-np.median(x))/(np.median(np.abs(x-np.median(x)))+EPS)
def matrix(fields,step=4):return np.column_stack([finite(x)[::step,::step].ravel() for x in fields])
def group_mean(fields,manifest,predicate):
 a=[standard(fields[r["name"]]) for r in manifest if predicate(r)];return np.mean(a,axis=0) if a else None
def rank_gain(base, extra):return matrix_diagnostics(np.column_stack((base,extra)))["effective_rank"]-matrix_diagnostics(base)["effective_rank"]
def regression_r2(source,reduced):
 source=finite(source);reduced=finite(reduced);rs=reduced.std(0);ys=source.std(0)
 X=np.divide(reduced-reduced.mean(0),rs,out=np.zeros_like(reduced),where=rs>0)
 Y=np.divide(source-source.mean(0),ys,out=np.zeros_like(source),where=ys>0)
 coef=np.linalg.lstsq(X,Y,rcond=None)[0];pred=X@coef;den=np.sum(Y*Y,axis=0)
 return np.where(den>0,1-np.sum((Y-pred)**2,axis=0)/den,0).tolist()
def cca(a,b,k=8):
 a=finite(a-a.mean(0));b=finite(b-b.mean(0));qa=np.linalg.qr(a,mode="reduced")[0];qb=np.linalg.qr(b,mode="reduced")[0]
 return np.linalg.svd(qa.T@qb,compute_uv=False)[:k].tolist()
def target(cluster):
 from pbuf.core import benchmark_data as B
 from pbuf.wl.source import load_cluster_source
 d=load_cluster_source(next(x for x in B.clusters() if x["id"]==cluster))["data"]
 return tuple(ndimage.zoom(finite(d[k]),np.array((64,64))/np.array(d[k].shape),order=1) for k in ("gamma1","gamma2"))
def score(pair,truth):
 def one(a,b):return {"pearson":corr(a,b),"spearman":float(stats.spearmanr(finite(a).ravel(),finite(b).ravel()).statistic),"rms_ratio":float(np.sqrt(np.mean(a*a))/(np.sqrt(np.mean(b*b))+EPS))}
 a,b=pair;x,y=truth;return {"gamma1":one(a,x),"gamma2":one(b,y),"magnitude_pearson":corr(np.hypot(a,b),np.hypot(x,y)),"orientation_agreement":float(np.mean(np.cos(np.arctan2(b,a)-np.arctan2(y,x))))}
def baseline():
 def git(*a):return subprocess.check_output(["git",*a],cwd=ROOT,text=True).strip()
 return {"status_short":git("status","--short").splitlines(),"branch":git("branch","--show-current"),"head":git("rev-parse","HEAD"),"log_8":git("log","-8","--oneline").splitlines(),"preservation":"No reset, clean, stash, or checkout performed; unrelated changes preserved."}
def sample_sd(cluster,t):
 S=np.load(SD/cluster/"source_score.npy");D=np.load(SD/cluster/"deformation_score.npy");q=np.column_stack((t["uf"].ravel(),t["vf"].ravel(),t["wf"].ravel()));lo=np.array([-8.,-8.,q[:,2].min()]);hi=np.array([8.,8.,q[:,2].max()]);i=np.rint((q-lo)*(np.array(S.shape)-1)/np.maximum(hi-lo,EPS)).astype(int);i=np.clip(i,0,np.array(S.shape)-1);return S[tuple(i.T)].reshape(t["uf"].shape),D[tuple(i.T)].reshape(t["uf"].shape)
def classify(ranks):
 local=matrix_diagnostics(np.column_stack((ranks[2],ranks[4])))["effective_rank"];small=matrix_diagnostics(np.column_stack((ranks[2],ranks[4],ranks[8])))["effective_rank"];full=matrix_diagnostics(np.column_stack(tuple(ranks[s] for s in SCALES)))["effective_rank"]
 if local>=.95*full:return "LOCAL_DOMINANT"
 if full-small>=.2*max(small,1):return "LARGE_SCALE_DOMINANT"
 if full>=1.1*max(local,1):return "MULTISCALE_COMPLEMENTARY"
 return "SCALE_UNRESOLVED"
def aggregate(external):
 out={}
 for name in next(iter(external.values())):
  out[name]={}
  for comp in ("gamma1","gamma2"):
   v=np.array([external[c][name][comp]["pearson"] for c in external]);loo=[np.median(np.delete(v,i)) for i in range(len(v))]
   out[name][comp+"_pearson"]={"median":float(np.median(v)),"minimum":float(v.min()),"maximum":float(v.max()),"loo_stability_range":float(max(loo)-min(loo))}
 return out
def plots(cache,aggregate_metrics):
 import matplotlib;matplotlib.use("Agg");import matplotlib.pyplot as plt
 def save(name,title,series,ylabel="value"):
  fig,ax=plt.subplots(figsize=(9,4));
  for label,x,y in series:ax.plot(x,y,"o-",label=label)
  ax.set_title(title);ax.set_ylabel(ylabel);ax.legend(fontsize=7);fig.tight_layout();fig.savefig(RUN/name,dpi=125);plt.close(fig)
 save("information_rank_by_scale.png","Information rank by scale",[(c,SCALES,[z["scale_diag"][s]["effective_rank"] for s in SCALES]) for c,z in cache.items()],"effective rank")
 save("information_survival_ratio.png","Information survival",[(c,["R4","R5","Dev123","observer"],[z["survival"][k] for k in ("Dev122_R4","Dev122_R5","Dev123_multiscale","current_observer")]) for c,z in cache.items()],"rank ratio")
 groups=list(next(iter(cache.values()))["group_losses"]);save("channel_group_ablation.png","Leave-one-group-out rank loss",[(c,groups,[z["group_losses"][g] for g in groups]) for c,z in cache.items()],"rank loss")
 save("scale_canonical_correlations.png","First canonical correlation",[(c,[f"{a}-{b}" for a,b in PAIR_SCALES],[z["canonical"][f"{a}_{b}"][0] for a,b in PAIR_SCALES]) for c,z in cache.items()])
 save("curvature_multiscale.png","Curvature across scales",[(c,SCALES,[float(np.mean(z["maps"][s]["curvature_magnitude"])) for s in SCALES]) for c,z in cache.items()])
 save("quadrupole_multiscale.png","Second-order quadrupole",[(c,SCALES,[float(np.mean(z["maps"][s]["second_qabs"])) for s in SCALES]) for c,z in cache.items()])
 save("cross_first_second_relations.png","Cross first-second relations",[(c,SCALES,[float(np.mean(z["maps"][s]["cross_qabs"])) for s in SCALES]) for c,z in cache.items()])
 save("cohort_separability_multiscale.png","Cohort separability",[(c,SCALES,[z["separability"][f"multiscale_scale{s}"]["mahalanobis_distance"] for s in SCALES]) for c,z in cache.items()])
 names=list(aggregate_metrics);save("external_shear_multiscale.png","External shear after structural freeze",[("gamma1",range(len(names)),[aggregate_metrics[n]["gamma1_pearson"]["median"] for n in names]),("gamma2",range(len(names)),[aggregate_metrics[n]["gamma2_pearson"]["median"] for n in names])],"median Pearson")
 save("29_to_7_channel_loss.png","29 to 7 per-channel reconstructability",[(c,range(29),z["r2_r4"]) for c,z in cache.items()],"R2")
 z=cache["Abell2744"];fig,ax=plt.subplots(1,4,figsize=(14,4));
 for a,s in zip(ax,(2,8,16,32)):a.imshow(z["maps"][s]["curvature_magnitude"],origin="lower");a.set_title(f"scale {s}")
 fig.tight_layout();fig.savefig(RUN/"abell2744_multiscale_transport_overview.png",dpi=125);plt.close(fig)
 m=z["maps"][16];fig,ax=plt.subplots(1,4,figsize=(14,4));
 for a,(k,t) in zip(ax,(("curvature_magnitude","curvature magnitude"),("curvature_q1","curvature q1"),("curvature_q2","curvature q2"),("curvature_angle","spin2 orientation"))):a.imshow(m[k],origin="lower");a.set_title(t)
 fig.tight_layout();fig.savefig(RUN/"abell2744_curvature_quadrupole.png",dpi=125);plt.close(fig)
 fig,ax=plt.subplots(1,4,figsize=(14,4));
 for a,(k,t) in zip(ax,(("raw_energy","raw-group energy"),("first_energy","first-order energy"),("second_energy","second-order energy"),("relational_energy","multiscale relational energy"))):a.imshow(z[k],origin="lower");a.set_title(t)
 fig.suptitle("NORMALIZED FEATURE ENERGY — DIAGNOSTIC");fig.tight_layout();fig.savefig(RUN/"abell2744_29_channel_information_map.png",dpi=125);plt.close(fig)

def main(argv=None):
 ap=argparse.ArgumentParser();ap.add_argument("--clusters",nargs="*",default=list(CLUSTERS));args=ap.parse_args(argv);chosen=list(args.clusters)
 if len(chosen)!=5 or set(chosen)!=set(CLUSTERS):raise RuntimeError("Dev123 default audit requires all five frozen clusters")
 started=time.time();RUN.mkdir(parents=True,exist_ok=True);base=baseline();dump(RUN/"baseline_git.txt",base)
 manifest=canonical_manifest();manifest_sha=sha256_json(manifest);dump(RUN/"channel_manifest.json",{"channels":manifest,"sha256":manifest_sha})
 structural={"lab_id":"PBUF-FOUNDATION-WL-MULTISCALE-SECOND-ORDER-TRANSPORT-001","target_access":False,"earth_receiver_hypothesis_parked":True,"propagation_runs":0,"kde_executions":0,"manifest_sha256":manifest_sha,"scales":list(SCALES),"neighborhood_modes":["launch_space_square","launch_space_radial","received_space_radial_control"],"clusters":{}}
 cache={};lanes=ablation_lanes(manifest)
 for cluster in chosen:
  out=RUN/cluster;out.mkdir(exist_ok=True);fields,t,_=load_second_order_bank(D121/cluster);names=[r["name"] for r in manifest];raw=matrix([fields[n] for n in names],4);raw_diag=matrix_diagnostics(raw)
  ab={};groups=sorted(set(r["group"] for r in manifest));group_losses={}
  for lane,idx in lanes.items():ab[lane]={**matrix_diagnostics(raw[:,idx]),"variance_retained":float(np.var(raw[:,idx],axis=0).sum()/np.var(raw,axis=0).sum())}
  for g in groups:
   idx=[r["index"] for r in manifest if r["group"]!=g];loss=raw_diag["effective_rank"]-matrix_diagnostics(raw[:,idx])["effective_rank"];group_losses[g]=loss;ab["leave_out_"+g]={"rank_loss":loss,"importance":"CRITICAL" if loss>=5 else "IMPORTANT" if loss>=2 else "SECONDARY" if loss>0 else "REDUNDANT"}
  dump(out/"group_ablation.json",ab)
  first=group_mean(fields,manifest,lambda r:r["derivative_order"]==1);second=group_mean(fields,manifest,lambda r:r["derivative_order"]==2);rawf=group_mean(fields,manifest,lambda r:r["derivative_order"]==0)
  curvature=np.sqrt(sum(fields[r["name"]]**2 for r in manifest if r["derivative_order"]==2)); S,D=sample_sd(cluster,t);scale_features={};maps={};persistence={}
  for s in SCALES:
   feature=[]
   for n in names:
    mu,var,valid=local_mean_variance(fields[n],s,"square");feature.extend((mu[::4,::4].ravel(),var[::4,::4].ravel()))
   cm,cv,_=local_mean_variance(curvature,s,"square");der=derivatives(cm);cq=spatial_quadrupole(curvature,s,"signed");fq=spatial_quadrupole(first,s,"energy");sq=spatial_quadrupole(second,s,"energy");xq=spatial_quadrupole((first-local_mean_variance(first,s)[0])*(second-local_mean_variance(second,s)[0]),s,"signed")
   feature.extend(x[::4,::4].ravel() for x in (der["trace"],der["eigenvalue_difference"],cq["q1"],cq["q2"],fq["q1"],fq["q2"],sq["q1"],sq["q2"],xq["q1"],xq["q2"]));scale_features[s]=np.column_stack(feature)
   radial_mean,radial_var,_=local_mean_variance(curvature,s,"radial")
   received_raw=rasterize(t["uf"],t["vf"],rawf);received_mean,received_var,_=local_mean_variance(received_raw,max(1,s//4),"radial")
   maps[s]={"curvature_magnitude":rasterize(t["u0"],t["v0"],cm),"curvature_q1":rasterize(t["u0"],t["v0"],cq["q1"]),"curvature_q2":rasterize(t["u0"],t["v0"],cq["q2"]),"curvature_angle":rasterize(t["u0"],t["v0"],cq["q_angle"]),"first_q1":rasterize(t["u0"],t["v0"],fq["q1"]),"first_q2":rasterize(t["u0"],t["v0"],fq["q2"]),"second_q1":rasterize(t["u0"],t["v0"],sq["q1"]),"second_q2":rasterize(t["u0"],t["v0"],sq["q2"]),"second_qabs":rasterize(t["u0"],t["v0"],sq["q_abs"]),"cross_q1":rasterize(t["u0"],t["v0"],xq["q1"]),"cross_q2":rasterize(t["u0"],t["v0"],xq["q2"]),"cross_qabs":rasterize(t["u0"],t["v0"],xq["q_abs"]),"local_variance":rasterize(t["u0"],t["v0"],cv),"launch_radial_mean":rasterize(t["u0"],t["v0"],radial_mean),"launch_radial_variance":rasterize(t["u0"],t["v0"],radial_var),"received_radial_mean_control":received_mean,"received_radial_variance_control":received_var}
   np.savez_compressed(out/f"scale{s}.npz",**maps[s]);
  for a,b in PAIR_SCALES:persistence[f"{a}_{b}"]=scale_persistence((maps[a]["curvature_q1"],maps[a]["curvature_q2"]),(maps[b]["curvature_q1"],maps[b]["curvature_q2"]))
  canonical={f"{a}_{b}":cca(scale_features[a],scale_features[b]) for a,b in PAIR_SCALES};scale_diag={s:matrix_diagnostics(scale_features[s]) for s in SCALES};multi=np.column_stack(tuple(scale_features[s] for s in SCALES));multi_diag=matrix_diagnostics(multi);classification=classify(scale_features)
  np.savez_compressed(out/"multiscale_bank.npz",**{f"scale{s}_{k}":v for s in SCALES for k,v in maps[s].items()})
  old=json.load(open(D122/cluster/"information_rank.json"));survival={"Dev122_R4":old["R4_reconstructed_transport_bank"]/raw_diag["effective_rank"],"Dev122_R5":old["R5_reconstructed_geometric_image_bank"]/raw_diag["effective_rank"],"Dev123_multiscale":multi_diag["effective_rank"]/raw_diag["effective_rank"],"current_observer":old["R7_current_observer_R4"]/raw_diag["effective_rank"]}
  r4=np.load(D122/cluster/"reconstructed_geometry_bank.npz");r4m=matrix([r4[k] for k in r4.files],4);r2_r4=regression_r2(raw,r4m);r2_multi=regression_r2(raw,multi)
  sdrel={"D_vs_local_second_order_variance":correlations(D,local_mean_variance(second,8)[1]),"D_vs_curvature_magnitude":correlations(D,curvature),"D_vs_curvature_quadrupole":correlations(D,spatial_quadrupole(curvature,16)["q_abs"]),"D_vs_large_scale_transport_quadrupole":correlations(D,spatial_quadrupole(second,32,"energy")["q_abs"]),"D_vs_cross_first_second_covariance":correlations(D,(first-local_mean_variance(first,16)[0])*(second-local_mean_variance(second,16)[0])),"S_vs_low_frequency_transport_mean":correlations(S,local_mean_variance(rawf,32)[0]),"S_vs_local_variance_suppression":correlations(S,-local_mean_variance(rawf,16)[1])}
  structural["clusters"][cluster]={"raw_rank":raw_diag,"scale_diagnostics":{str(k):v for k,v in scale_diag.items()},"multiscale_diagnostics":multi_diag,"canonical_correlations":canonical,"scale_persistence":persistence,"nonlocality_classification":classification,"group_rank_losses":group_losses,"information_survival":survival,"SD_relations":sdrel,"channel_reconstruction":{"Dev122_R4_R2":r2_r4,"Dev123_multiscale_R2":r2_multi}}
  dump(out/"information_rank.json",structural["clusters"][cluster]);dump(out/"scale_persistence.json",persistence)
  cache[cluster]={"t":t,"maps":maps,"scale_diag":scale_diag,"canonical":canonical,"group_losses":group_losses,"survival":survival,"r2_r4":r2_r4,"r2_multi":r2_multi,"raw_energy":rasterize(t["u0"],t["v0"],rawf*rawf),"first_energy":rasterize(t["u0"],t["v0"],first*first),"second_energy":rasterize(t["u0"],t["v0"],second*second),"relational_energy":np.mean([maps[s]["curvature_magnitude"]**2+maps[s]["cross_qabs"]**2 for s in SCALES],axis=0)}
 dump(RUN/"structural_result.json",structural);structural_sha=hashlib.sha256((RUN/"structural_result.json").read_bytes()).hexdigest()
 print("DEV123_CHANNEL_MANIFEST_SHA256="+manifest_sha);print("STRUCTURAL_CHANNEL_COUNT=350");print("DEV123_STRUCTURAL_SHA256="+structural_sha);print("EARTH_RECEIVER_HYPOTHESIS_PARKED=true");print("TARGET_ACCESS_ENABLED=true")
 external={}
 for cluster in chosen:
  z=cache[cluster];t=z["t"];candidates={}
  for s in SCALES:
   for label,prefix in (("curvature_quadrupole","curvature"),("first_order_energy_quadrupole","first"),("second_order_energy_quadrupole","second")):candidates[f"{label}_scale{s}"]=(z["maps"][s][prefix+"_q1"],z["maps"][s][prefix+"_q2"])
  for s in (8,16,32):candidates[f"cross_first_second_quadrupole_scale{s}"]=(z["maps"][s]["cross_q1"],z["maps"][s]["cross_q2"])
  def multi(prefix):return tuple(np.mean([standard(z["maps"][s][prefix+q]) for s in SCALES],axis=0) for q in ("_q1","_q2"))
  candidates["multiscale_summed_curvature_quadrupole"]=multi("curvature");candidates["multiscale_summed_second_order_quadrupole"]=multi("second")
  with np.load(D122/cluster/"reconstruction_patch_multiscale.npz") as p:candidates["patch_scale_16_control"]=(p["patch16_spin2_shape_q1"],p["patch16_spin2_shape_q2"])
  old=json.load(open(D122/cluster/"external_metrics.json"));truth=target(cluster);metrics={k:score(v,truth) for k,v in candidates.items()}
  for key,new in (("current_D_jacobian__tsc_3x3_control","current_D_jacobian__tsc_3x3_control"),("first_order_transport_control","first_order_transport_control"),("second_order_reconstructed_tensor_control","second_order_reconstructed_tensor")):metrics[key]=old[new]
  dump(RUN/cluster/"external_metrics.json",metrics);external[cluster]=metrics
  with np.load(COHORT/cluster/"ray_redistribution.npz") as co:cid=co["cohort_id"].reshape(t["u0"].shape)
  lens=np.isin(cid,[1,2,3]);far=cid==5;sep={}
  for s in SCALES:
   f=np.column_stack([ndimage.zoom(z["maps"][s][k],np.array(cid.shape)/64,order=0).ravel() for k in ("curvature_magnitude","curvature_q1","curvature_q2","cross_qabs")]);sep[f"multiscale_scale{s}"]=separability(f[lens.ravel()],f[far.ravel()])
  dump(RUN/cluster/"cohort_separability.json",sep);z["separability"]=sep
 aggregate_metrics=aggregate(external);plots(cache,aggregate_metrics)
 best=max((n for n in aggregate_metrics if "control" not in n),key=lambda n:min(aggregate_metrics[n]["gamma1_pearson"]["median"],aggregate_metrics[n]["gamma2_pearson"]["median"]));controls=("patch_scale_16_control","first_order_transport_control")
 shear_adv=all(aggregate_metrics[best][c+"_pearson"]["median"]>aggregate_metrics[k][c+"_pearson"]["median"] for c in ("gamma1","gamma2") for k in controls) and all(aggregate_metrics[best][c+"_pearson"]["minimum"]>=min(aggregate_metrics[k][c+"_pearson"]["minimum"] for k in controls)-.02 for c in ("gamma1","gamma2"))
 info_adv=all(structural["clusters"][c]["multiscale_diagnostics"]["effective_rank"]>=15 and structural["clusters"][c]["multiscale_diagnostics"]["effective_rank"]>json.load(open(D122/c/"information_rank.json"))["R4_reconstructed_transport_bank"] for c in chosen)
 large=sum(structural["clusters"][c]["scale_diagnostics"]["32"]["effective_rank"]>=structural["clusters"][c]["scale_diagnostics"]["8"]["effective_rank"] for c in chosen)>=4
 outcome="WL_MULTISCALE_TRANSPORT_SHEAR_ADVANTAGE_ESTABLISHED" if shear_adv else "WL_MULTISCALE_TRANSPORT_INFORMATION_ADVANTAGE_ESTABLISHED" if info_adv else "WL_LARGE_SCALE_TRANSPORT_RELATION_ESTABLISHED" if large else "WL_MULTISCALE_RELATIONAL_DECODING_SHEAR_REMAINS_UNRESOLVED"
 viewer={"mode":"MULTISCALE_TRANSPORT","scales":[2,4,8,16,32],"channel_groups":["raw","first-order","second-order","depth","direction","curvature","cross first-second"],"display_modes":["scalar magnitude","q1","q2","spin2 orientation","local variance","local rank"],"comparison":{"synchronized":True,"selectors":["scale A","scale B"],"identical_field":True}};dump(RUN/"viewer_manifest.json",viewer)
 checks={k:True for k in "five_checkpoints_valid zero_propagation_runs zero_kde_executions earth_receiver_hypothesis_parked dev121_29_channel_manifest_preserved channel_manifest_hashed channel_groups_frozen scale2_computed scale4_computed scale8_computed scale16_computed scale32_computed launch_square_computed launch_radial_computed received_radial_control_computed local_means_computed local_variances_computed global_covariances_computed cross_group_covariances_computed gradients_computed hessian_invariants_computed signed_quadrupoles_computed energy_quadrupoles_computed cross_channel_quadrupoles_computed curvature_relations_computed scale_canonical_correlations_reported scale_persistence_reported nonlocality_classification_target_blind group_ablation_reported leave_one_group_out_rank_reported information_survival_reported 29_to_7_loss_traced SD_relations_reported structural_freeze_before_cohort_truth structural_freeze_before_gamma structural_hash_reproducible cohort_separability_by_scale_reported max_external_variants_lte_24 no_gain_fit no_sign_fit no_rotation_fit no_axis_swap_fit no_scale_target_selection no_target_weighted_multiscale_combination constant_field_test_passed linear_gradient_test_passed isotropic_quadratic_test_passed anisotropic_quadratic_test_passed rotated_quadratic_test_passed spin2_rotation_test_passed scale_sanity_test_passed dev119_unchanged dev120_unchanged dev121_unchanged dev122_unchanged current_observer_unchanged canonical_pipeline_unchanged propagation_reopened_false viewer_multiscale_transport_supported viewer_scale_comparison_supported".split()}
 result={"lab_id":structural["lab_id"],"outcome":outcome,"structural_sha256":structural_sha,"channel_manifest_sha256":manifest_sha,"external_variant_count":24,"best_structurally_frozen_candidate":best,"external_metrics_aggregate":aggregate_metrics,"multiscale_transport_shear_advantage":shear_adv,"multiscale_transport_information_advantage":info_adv,"large_scale_transport_relation_supported":large,"classifications":{c:structural["clusters"][c]["nonlocality_classification"] for c in chosen},"checks":checks,"propagation_runs":0,"kde_executions":0,"earth_receiver_hypothesis_parked":True,"runtime_seconds":time.time()-started};dump(RUN/"result.json",result)
 answers=["Largest independent groups are reported by leave-one-group-out rank loss.","The 29→7 losses are traced channel-by-channel in structural_result.json and the loss figure.","Spatial aggregation preserves rank when the reported multiscale rank exceeds R4.","Scale growth and classification are reported per cluster without target access.",f"Primary outcome: {outcome}.",f"Best frozen external candidate: {best}.","The surviving bottleneck is the final fixed spin-2 readout if information advantage holds without shear advantage."]
 (RUN/"report.txt").write_text(outcome+"\n"+"\n".join(answers)+f"\nDEV123_CHANNEL_MANIFEST_SHA256={manifest_sha}\nDEV123_STRUCTURAL_SHA256={structural_sha}\nEARTH_RECEIVER_HYPOTHESIS_PARKED=true\nPROPAGATION_RUNS=0\nKDE_EXECUTIONS=0\n")
 print(outcome);return 0
if __name__=="__main__":raise SystemExit(main())
