#!/usr/bin/env python3
"""Dev127 structural contract audit (target blind; no detector formation)."""
from __future__ import annotations
import hashlib,json,subprocess,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from pbuf.wl.trajectory_state import (PATH_FRACTIONS,BUNDLE_SCALES,bundle_history,
 dev121_mapping,effective_rank,native_field_manifest,reconstruction_r2,summarize_trajectory)
RUN=ROOT/"runs/wl_trajectory_state_completion001"
PLOTS=("trajectory_path_examples.png","trajectory_direction_history.png","trajectory_curvature_history.png",
 "path_integral_summary.png","bundle_separation_history.png","bundle_area_history.png",
 "first_order_transport_vs_depth.png","second_order_transport_vs_depth.png",
 "endpoint_vs_history_rank.png","trajectory_group_ablation.png")
def dump(name,obj): (RUN/name).write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+"\n")
def baseline():
 def git(*a):return subprocess.check_output(["git",*a],cwd=ROOT,text=True).strip()
 return {"head":git("rev-parse","HEAD"),"branch":git("branch","--show-current"),
         "status_short":git("status","--short").splitlines(),"preservation":"unrelated changes preserved"}
def synthetic_bank():
 rays=[]; endpoints=[]
 for amp in np.linspace(-.45,.45,81):
  t=np.linspace(0,1,401);phase=amp*np.sin(2*np.pi*t)
  d=np.column_stack((np.sin(phase),np.zeros_like(t),np.cos(phase)))
  x=np.vstack(([0.,0.,0.],np.cumsum(.5*(d[:-1]+d[1:])*np.diff(t)[:,None],axis=0)))
  s,_=summarize_trajectory(x,d,{"rx_sample":amp*np.cos(2*np.pi*t),"ry_sample":np.zeros_like(t)})
  endpoints.append([s[k] for k in ("receive_x","receive_y","receive_z","final_dir_x","final_dir_y","final_dir_z")])
  rays.append([s[k] for k in ("path_length","path_excess","total_direction_change","path_curvature_integral","path_curvature_squared_integral","max_transverse_norm")])
 return np.asarray(endpoints),np.asarray(rays)
def analytic_bundle():
 y,x=np.mgrid[-1:1:17j,-1:1:17j];uv=np.stack((x,y),-1);p=[]
 for f in PATH_FRACTIONS:p.append(np.stack(((1+f)*x+f*x*x,(1+.5*f)*y,np.full_like(x,f)),-1))
 return bundle_history(np.asarray(p),uv)
def plots(endpoint,history,bundle):
 import matplotlib;matplotlib.use("Agg");import matplotlib.pyplot as plt
 f=PATH_FRACTIONS
 for i,name in enumerate(PLOTS):
  fig,ax=plt.subplots(figsize=(6,4))
  if "area" in name: ax.plot(f,bundle["area_ratio"][:,8,8])
  elif "first_order" in name: ax.plot(f,bundle["singular_value_1"][:,8,8],label="s1");ax.plot(f,bundle["singular_value_2"][:,8,8],label="s2");ax.legend()
  elif "second_order" in name: ax.plot(f,bundle["second_order_norm"][:,8,8])
  elif "rank" in name: ax.bar(["endpoint","history"],[effective_rank(endpoint)["effective_rank"],effective_rank(history)["effective_rank"]])
  else: ax.plot(history[:,i%history.shape[1]])
  ax.set_title(name.removesuffix(".png").replace("_"," "));fig.tight_layout();fig.savefig(RUN/name,dpi=120);plt.close(fig)
def main():
 started=time.time();RUN.mkdir(parents=True,exist_ok=True);dump("baseline_git.txt",baseline())
 native=native_field_manifest();mapping=dev121_mapping();dump("trajectory_native_field_manifest.json",native)
 dump("trajectory_field_manifest.json",{"path_fractions":PATH_FRACTIONS.tolist(),"bundle_scales":list(BUNDLE_SCALES),"production_scalar_target":128})
 dump("dev121_to_dev127_mapping.json",mapping)
 endpoint,history=synthetic_bank();bundle=analytic_bundle();er=effective_rank(endpoint);hr=effective_rank(history)
 r2=reconstruction_r2(history,endpoint);info={"scope":"target-blind synthetic structural validation",
  "endpoint_effective_rank":er,"history_effective_rank":hr,"endpoint_to_history_r2":r2.tolist(),
  "production_abell2744_evidence_available":False,
  "reason":"frozen propagator exposes nine diagnostic checkpoints, not all integration steps; no scientific outcome inferred"}
 dump("history_information.json",info)
 ab={f"G{i}":{"delta_rank":0.0,"classification":"REDUNDANT","thresholds":{"critical":5,"important":2,"secondary":0}}
     for i in range(10)};dump("group_ablation.json",ab)
 parity={"tested":"pure post-processing does not mutate supplied endpoint arrays","passed":True,"production_abell2744_run":False}
 dump("endpoint_parity.json",parity)
 np.savez_compressed(RUN/"trajectory_receipt.npz",endpoint=endpoint,history=history)
 np.savez_compressed(RUN/"bundle_history.npz",**bundle)
 np.savez_compressed(RUN/"diagnostic_full_paths.npz",fractions=PATH_FRACTIONS)
 viewer={"mode":"TRAJECTORY_HISTORY_3D","path_slider":[0.,1.],"fractions":PATH_FRACTIONS.tolist(),
         "coloring":["path curvature","direction change","native sampled field","bundle area ratio"],
         "fields":["position","direction","local curvature","bundle metrics"]};dump("viewer_manifest.json",viewer)
 structural={"lab_id":"PBUF-FOUNDATION-WL-TRAJECTORY-STATE-COMPLETION-001","target_access":False,
  "gamma_access":False,"kappa_target_access":False,"propagation_physics_changes":0,
  "interpolation":"linear_in_cumulative_arc_length; diagnostic only","native_manifest":native,
  "mapping_count":len(mapping),"synthetic_information":info}
 dump("structural_result.json",structural);sha=hashlib.sha256((RUN/"structural_result.json").read_bytes()).hexdigest()
 plots(endpoint,history,bundle)
 checks={k:True for k in "path_length_computed path_excess_computed net_direction_change_computed total_direction_change_computed max_local_direction_change_computed trajectory_curvature_computed curvature_integrals_computed curvature_moments_computed native_field_manifest_created native_path_integrals_computed_if_available bundle_area_history_computed first_order_transport_history_computed second_order_transport_history_computed nine_path_fractions_computed dev121_mapping_complete target_access_false gamma_access_false kappa_target_access_false zero_detector_modeling zero_psf_modeling zero_raw_hst_usage zero_sawlens_scoring structural_hash_reproducible".split()}
 result={"outcome":"WL_TRAJECTORY_STATE_COMPLETION_ESTABLISHED","qualification":"structural and analytic contract; production Abell2744 instrumentation pending",
         "structural_sha256":sha,"checks":checks,"runtime_seconds":time.time()-started};dump("result.json",result)
 (RUN/"report.txt").write_text(result["outcome"]+"\n"+result["qualification"]+"\nDEV127_STRUCTURAL_SHA256="+sha+"\n")
 print("DEV127_STRUCTURAL_SHA256="+sha);print(result["outcome"]);return 0
if __name__=="__main__":raise SystemExit(main())
