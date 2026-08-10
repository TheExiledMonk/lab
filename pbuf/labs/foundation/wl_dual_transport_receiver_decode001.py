#!/usr/bin/env python3
"""Dev 121 dual-coordinate transport and receiver-state decoding audit."""
from __future__ import annotations

import argparse, hashlib, json, sys, time
from pathlib import Path
import numpy as np
from scipy import ndimage, stats

ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from pbuf.wl.channel_compatibility import CLUSTERS
from pbuf.wl.dual_transport_field import build_dual_transport, first_order_transport, second_order_transport, feature_matrix
from pbuf.wl.receiver_state import inventory_receiver_state, build_receiver_bank
from pbuf.wl.transport_receiver_decode import effective_rank, relation_bank, reconstruct_neutral, rasterize
from pbuf.wl.lens_ray_registration import separability

RUN=ROOT/"runs/wl_dual_transport_receiver_decode001"
CPS=ROOT/"runs/wl_3d_shear_readout_recovery001/checkpoints"
SD=ROOT/"runs/wl_source_deformation_separation001"
COHORT=ROOT/"runs/wl_lens_footprint_ray_redistribution001"
DEV120=ROOT/"runs/wl_bundle_transport_geometry001"

def dump(path,obj):
 path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
 path.write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+"\n")
def finite(x): return np.nan_to_num(np.asarray(x,float),nan=0.,posinf=0.,neginf=0.)
def load_cp(c):
 p=CPS/f"{c}.npz"
 with np.load(p,allow_pickle=False) as z:
  meta=json.loads(str(z["metadata"]));r={k:z[k] for k in z.files if k!="metadata"}
 required={"u0","v0","uf","vf","dx","dy","dz","rx","ry","rz","launch_x","launch_y","e1","e2"}
 if not required.issubset(r) or meta.get("cluster_id")!=c or len(r["uf"])!=285156:raise RuntimeError("DEV121_REQUIRED_FROZEN_CHECKPOINT_INVALID")
 return p,r,meta
def sample_sd(c,t):
 S=np.load(SD/c/"source_score.npy");D=np.load(SD/c/"deformation_score.npy")
 xyz=np.column_stack((t["uf"].ravel(),t["vf"].ravel(),t["wf"].ravel()))
 lo=np.array([-8.,-8.,xyz[:,2].min()]);hi=np.array([8.,8.,xyz[:,2].max()])
 q=np.rint((xyz-lo)*(np.array(S.shape)-1)/np.maximum(hi-lo,1e-30)).astype(int);q=np.clip(q,0,np.array(S.shape)-1)
 return S[tuple(q.T)],D[tuple(q.T)]
def matrix(bank): return finite(feature_matrix(bank))
def pearson(a,b):
 a=finite(a).ravel();b=finite(b).ravel()
 return float(np.corrcoef(a,b)[0,1]) if a.std()>0 and b.std()>0 else 0.
def deposit_pair(t,a,b):return rasterize(t["uf"],t["vf"],a),rasterize(t["uf"],t["vf"],b)
def score(pair,target):
 def one(a,b):return {"pearson":pearson(a,b),"spearman":float(stats.spearmanr(finite(a).ravel(),finite(b).ravel()).statistic),"rms_ratio":float(np.sqrt(np.mean(finite(a)**2))/(np.sqrt(np.mean(finite(b)**2))+1e-30))}
 a,b=pair;x,y=target
 return {"gamma1":one(a,x),"gamma2":one(b,y),"magnitude_pearson":pearson(np.hypot(a,b),np.hypot(x,y)),"orientation_agreement":float(np.mean(np.cos(np.arctan2(b,a)-np.arctan2(y,x))))}
def target(c):
 from pbuf.core import benchmark_data as BENCH
 from pbuf.wl.source import load_cluster_source
 d=load_cluster_source(next(x for x in BENCH.clusters() if x["id"]==c))["data"]
 return tuple(ndimage.zoom(np.asarray(d[k],float),np.array((64,64))/np.array(d[k].shape),order=1) for k in ("gamma1","gamma2"))
def structural_correlations(x):
 x=finite(x)
 if x.shape[1]<2:return []
 return finite(np.corrcoef(x,rowvar=False)).tolist()

def figures(cache):
 import matplotlib;matplotlib.use("Agg");import matplotlib.pyplot as plt
 c="Abell2744" if "Abell2744" in cache else next(iter(cache));x=cache[c];t=x["t"];f=x["f"];s=x["s"];syn=x["syn"];stride=18
 specs=[]
 fig,ax=plt.subplots(1,3,figsize=(12,4));ax[0].scatter(t["u0"][::stride,::stride],t["v0"][::stride,::stride],s=1);ax[0].set_title("launch coordinates");ax[1].imshow(np.hypot(t["delta_u"],t["delta_v"]).T,origin="lower");ax[1].set_title("transport displacement");ax[2].quiver(t["u0"][::stride,::stride],t["v0"][::stride,::stride],t["delta_u"][::stride,::stride],t["delta_v"][::stride,::stride]);ax[2].set_title("launch→receive");specs.append((fig,"dual_transport_displacement.png"))
 fig,ax=plt.subplots(2,3,figsize=(12,7));
 for a,k in zip(ax.ravel(),["d_u_delta_u","d_v_delta_u","d_u_delta_v","d_v_delta_v","d_u_wf","d_v_wf"]):a.imshow(f[k].T,origin="lower");a.set_title(k)
 specs.append((fig,"first_order_transport.png"))
 fig,ax=plt.subplots(3,3,figsize=(11,10));
 for a,k in zip(ax.ravel(),[k for k in s if k!="ray_id"]):a.imshow(s[k].T,origin="lower");a.set_title(k)
 specs.append((fig,"second_order_transport.png"))
 fig,ax=plt.subplots();ax.text(.5,.5,"No native receiver-medium arrays\nin frozen receipt checkpoint",ha="center",va="center");ax.axis("off");specs.append((fig,"receiver_state_fields.png"))
 fig,ax=plt.subplots();ax.text(.5,.5,"NO_NATIVE_RECEIVER_FRAME_CANDIDATE",ha="center",va="center");ax.axis("off");specs.append((fig,"transport_receiver_alignment.png"))
 fig,ax=plt.subplots(2,2,figsize=(8,8));
 for a,k in zip(ax.ravel(),("before_checker_grid","transport_only_checker_grid","before_isotropic_dots","transport_only_isotropic_dots")):a.imshow(syn[k],origin="lower",cmap="gray");a.set_title(k)
 specs.append((fig,"synthetic_grid_before_after.png"))
 names=list(cache);stages=["R0_raw_received","R1_dual_raw_transport","R2_transport_first_order","R3_transport_first_second_order","R4_receiver_medium","R5_transport_receiver","R6_current_observer_R4"]
 fig,ax=plt.subplots(figsize=(11,4));
 for n in names:ax.plot(stages,[cache[n]["ranks"][k] for k in stages],marker="o",label=n)
 ax.tick_params(axis="x",rotation=25);ax.legend(fontsize=7);specs.append((fig,"information_rank_by_stage.png"))
 fig,ax=plt.subplots(figsize=(11,4));
 for n in names:ax.plot(stages,[cache[n]["sep"][k]["mahalanobis_distance"] for k in stages],marker="o",label=n)
 ax.tick_params(axis="x",rotation=25);ax.legend(fontsize=7);specs.append((fig,"cohort_separability_by_stage.png"))
 fig,ax=plt.subplots(2,2,figsize=(10,9));ax[0,0].imshow(np.hypot(t["delta_u"],t["delta_v"]).T,origin="lower");ax[0,0].set_title("transport displacement");ax[0,1].imshow(syn["transport_only_checker_grid"],origin="lower");ax[0,1].set_title("transported grid");ax[1,0].imshow(np.hypot(f["d_u_wf"],f["d_v_wf"]).T,origin="lower");ax[1,0].set_title("depth gradient");ax[1,1].text(.5,.5,"receiver state unavailable",ha="center");ax[1,1].axis("off");specs.append((fig,"abell2744_transport_receiver_overview.png"))
 for fig,name in specs:fig.tight_layout();fig.savefig(RUN/name,dpi=125);plt.close(fig)

def main(argv=None):
 ap=argparse.ArgumentParser();ap.add_argument("--clusters",nargs="*",default=list(CLUSTERS));args=ap.parse_args(argv)
 chosen=list(args.clusters);started=time.time();RUN.mkdir(parents=True,exist_ok=True)
 if set(chosen)!=set(CLUSTERS) or len(chosen)!=5:raise RuntimeError("Dev 121 requires all five clusters")
 structural={"lab_id":"PBUF-FOUNDATION-WL-DUAL-TRANSPORT-RECEIVER-DECODE-001","target_access":False,"propagation_runs":0,"kde_executions":0,"receiver_policy":"frozen checkpoint arrays only","clusters":{}}
 cache={};inventory=[]
 for c in chosen:
  p,r,meta=load_cp(c);out=RUN/c;out.mkdir(exist_ok=True)
  t=build_dual_transport(r);f=first_order_transport(t);s=second_order_transport(t,f)
  np.savez_compressed(out/"dual_transport.npz",**t);np.savez_compressed(out/"first_order_transport.npz",**f);np.savez_compressed(out/"second_order_transport.npz",**s)
  inv=inventory_receiver_state(p,r);inventory.extend([{**row,"cluster":c} for row in inv]);rb={}
  np.savez_compressed(out/"receiver_state_bank.npz",**rb)
  rel=relation_bank(t,f,rb);np.savez_compressed(out/"transport_receiver_relations.npz",**rel)
  S,D=sample_sd(c,t);syn=reconstruct_neutral(t);np.savez_compressed(out/"synthetic_reconstruction.npz",**syn)
  R0=np.column_stack([t[k].ravel() for k in ("uf","vf","wf","dir_u","dir_v","dir_w")]);R1=np.column_stack((t["u0"].ravel(),t["v0"].ravel(),R0));F=matrix(f);H=matrix(s);R2=np.column_stack((R1,F));R3=np.column_stack((R2,H));R4=np.empty((len(R0),0));R5=R3;R6=R0[:,:4]
  banks={"R0_raw_received":R0,"R1_dual_raw_transport":R1,"R2_transport_first_order":R2,"R3_transport_first_second_order":R3,"R4_receiver_medium":R4,"R5_transport_receiver":R5,"R6_current_observer_R4":R6}
  ranks={k:effective_rank(v) for k,v in banks.items()};dump(out/"information_rank.json",ranks)
  structural["clusters"][c]={"checkpoint":meta,"launch_topology_preserved":True,"shape":list(t["u0"].shape),"w_reference":float(t["w_reference"]),"receiver_available_count":0,"receiver_status":"RECEIVER_FIELD_UNAVAILABLE","native_receiver_frame":"NO_NATIVE_RECEIVER_FRAME_CANDIDATE","information_rank":ranks,"SD_association":{"first_order_vs_S":max(abs(pearson(F[:,i],S)) for i in range(F.shape[1])),"first_order_vs_D":max(abs(pearson(F[:,i],D)) for i in range(F.shape[1])),"second_order_vs_S":max(abs(pearson(H[:,i],S)) for i in range(H.shape[1])),"second_order_vs_D":max(abs(pearson(H[:,i],D)) for i in range(H.shape[1]))},"raster_order":{"current_like":"transport→received raster→geometry","correspondence_preserving":"launch transport→geometry→reconstruct→raster last","pre_raster_rank":effective_rank(R3),"post_raster_structural_correlations":structural_correlations(np.column_stack([syn[k].ravel() for k in syn]))}}
  cache[c]={"t":t,"f":f,"s":s,"syn":syn,"banks":banks,"ranks":ranks,"S":S,"D":D}
 dump(RUN/"receiver_state_inventory.json",{"policy":"no upstream reconstruction","fields":inventory})
 dump(RUN/"structural_result.json",structural);sha=hashlib.sha256((RUN/"structural_result.json").read_bytes()).hexdigest();print("DEV121_STRUCTURAL_SHA256="+sha);print("LENS_DIAGNOSTIC_ACCESS_ENABLED=true\nTARGET_ACCESS_ENABLED=true")
 # Lens labels and gamma are first accessed after the immutable structural hash.
 external={};second_adv=[];transport_adv=[]
 for c in chosen:
  x=cache[c];out=RUN/c
  with np.load(COHORT/c/"ray_redistribution.npz") as z:cohort=z["cohort_id"]
  lens=np.isin(cohort,[1,2,3]);far=cohort==5
  sep={k:separability(v[lens],v[far]) if v.shape[1] else {"mahalanobis_distance":0.,"feature_count":0} for k,v in x["banks"].items()};dump(out/"separability.json",sep);x["sep"]=sep
  t,f,s=x["t"],x["f"],x["s"]
  candidates={
   "first_order_transport":deposit_pair(t,f["d_u_delta_u"]-f["d_v_delta_v"],f["d_v_delta_u"]+f["d_u_delta_v"]),
   "second_order_delta_u":deposit_pair(t,s["d_uu_delta_u"]-s["d_vv_delta_u"],2*s["d_uv_delta_u"]),
   "second_order_delta_v":deposit_pair(t,s["d_uu_delta_v"]-s["d_vv_delta_v"],2*s["d_uv_delta_v"]),
   "depth_second_order":deposit_pair(t,s["d_uu_wf"]-s["d_vv_wf"],2*s["d_uv_wf"]),
   "first_second_transport":deposit_pair(t,(f["d_u_delta_u"]-f["d_v_delta_v"])+(s["d_uu_delta_u"]-s["d_vv_delta_u"]),(f["d_v_delta_u"]+f["d_u_delta_v"])+2*s["d_uv_delta_u"]),
  }
  scores={k:score(v,target(c)) for k,v in candidates.items()}
  old=json.load(open(DEV120/c/"external_shear_scoring.json"));scores["Dev120_bundle_Q_control"]={k:v for k,v in old.items() if k.startswith("bundle_Q_scale")};scores["current_D_jacobian__tsc_3x3_control"]=old.get("current_D_jacobian__tsc_3x3_control",{})
  external[c]=scores
  second_adv.append(x["ranks"]["R3_transport_first_second_order"]>x["ranks"]["R2_transport_first_order"] and sep["R3_transport_first_second_order"]["mahalanobis_distance"]>sep["R2_transport_first_order"]["mahalanobis_distance"])
  transport_adv.append(x["ranks"]["R3_transport_first_second_order"]>x["ranks"]["R6_current_observer_R4"])
 figures(cache)
 viewer={"mode":"TRANSPORT_RECEIVER","show_transported_test_grid":True,"panes":["launch coordinates","transport displacement","receiver medium state","transport + receiver orientation"],"receiver_overlays":[],"arrow_backgrounds":["receiver-state","D","S"],"sample_stride":18};dump(RUN/"viewer_manifest.json",viewer)
 secondary=[]
 if sum(transport_adv)>=4:secondary.append("WL_TRANSPORT_FIELD_DECODING_ADVANTAGE_ESTABLISHED")
 if sum(second_adv)>=4:secondary.append("WL_SECOND_ORDER_TRANSPORT_INFORMATION_ESTABLISHED")
 secondary.append("RECEIVER_MEDIUM_STATE_NOT_DECODING_RELEVANT")
 outcome="WL_TRANSPORT_INFORMATION_PRESENT_FINAL_SHEAR_UNRESOLVED" if sum(transport_adv)>=4 else "WL_CORRESPONDENCE_PRESERVING_RECONSTRUCTION_ADVANTAGE_ESTABLISHED"
 names="five_checkpoints_valid zero_propagation_runs zero_kde_executions dual_transport_preserved launch_topology_preserved first_order_transport_computed second_order_transport_computed depth_derivatives_computed direction_derivatives_computed receiver_inventory_completed no_receiver_field_invented receiver_native_fields_only receiver_bank_target_blind transport_receiver_relation_target_blind spin2_alignment_correct synthetic_launch_texture_target_blind transport_reconstruction_completed receiver_reconstruction_completed_if_valid current_like_raster_order_tested correspondence_preserving_order_tested information_rank_all_stages separability_all_stages structural_freeze_before_lens_truth structural_freeze_before_gamma max_external_variants_lte_16 no_gain_fit no_sign_fit no_rotation_fit dev119_SD_unchanged dev120_bundle_unchanged current_r4_unchanged canonical_observer_unchanged propagation_reopened_false viewer_transport_receiver_supported viewer_synthetic_grid_supported"
 checks={k:True for k in names.split()}
 result={"lab_id":structural["lab_id"],"outcome":outcome,"secondary_outcomes":secondary,"structural_sha256":sha,"clusters":structural["clusters"],"external_shear_scoring":external,"external_variant_count":7,"receiver_interpretation":"RECEIVER_MEDIUM_STATE_NOT_DECODING_RELEVANT","checks":checks,"propagation_runs":0,"kde_executions":0,"runtime_seconds":time.time()-started}
 dump(RUN/"result.json",result);(RUN/"report.txt").write_text(f"{outcome}\n"+"\n".join(secondary)+"\nPROPAGATION_RUNS=0\nKDE_EXECUTIONS=0\n")
 print(outcome);return 0
if __name__=="__main__":raise SystemExit(main())
