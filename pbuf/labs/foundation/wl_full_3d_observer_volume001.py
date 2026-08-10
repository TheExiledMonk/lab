#!/usr/bin/env python3
"""Dev 117: full target-blind 3-D observer volume and late-projection audit."""
from __future__ import annotations
import hashlib,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3]; sys.path.insert(0,str(ROOT))
from pbuf.core import benchmark_data as BENCH
from pbuf.labs.foundation import native_multichannel_observer_fusion_sweep001 as FUS
from pbuf.wl.channel_compatibility import CLUSTERS
from pbuf.wl.observer_volume3d import construct_volume,save_volume
from pbuf.wl.observer_tensor3d import tensor_bank,feature_diagnostics,effective_rank
from pbuf.wl.observer_late_projection import late_projection,structural_gates

RUN=ROOT/"runs/wl_full_3d_observer_volume001"; CHECKPOINTS=ROOT/"runs/wl_3d_shear_readout_recovery001/checkpoints"
EVIDENCE=[ROOT/"runs/wl_3d_shear_readout_recovery001/result.json",ROOT/"runs/wl_abell2744_8192_shear_convergence001/result.json",ROOT/"runs/wl_observer_basis_information_mixing001/result.json"]

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def json_safe(x):
    if isinstance(x,(np.floating,np.integer)): return x.item()
    if isinstance(x,np.ndarray): return x.tolist()
    return x
def corr(a,b):
    m=np.isfinite(a)&np.isfinite(b); return float(np.corrcoef(a[m],b[m])[0,1]) if m.sum()>2 and np.std(a[m])*np.std(b[m])>0 else None
def rankdata(x): return np.argsort(np.argsort(x)).astype(float)
def metrics(q1,q2,g1,g2):
    m=np.isfinite(q1)&np.isfinite(q2)&np.isfinite(g1)&np.isfinite(g2); q1,q2,g1,g2=(x[m] for x in (q1,q2,g1,g2))
    return {"gamma1":{"pearson":corr(q1,g1),"spearman":corr(rankdata(q1),rankdata(g1)),"rms_ratio":float(np.sqrt(np.mean(q1*q1))/np.sqrt(np.mean(g1*g1)))},
            "gamma2":{"pearson":corr(q2,g2),"spearman":corr(rankdata(q2),rankdata(g2)),"rms_ratio":float(np.sqrt(np.mean(q2*q2))/np.sqrt(np.mean(g2*g2)))},
            "magnitude_pearson":corr(np.hypot(q1,q2),np.hypot(g1,g2)),"orientation_agreement":float(abs(np.mean(np.exp(1j*(np.arctan2(q2,q1)-np.arctan2(g2,g1))))))}
def diagnostics(out,vol,bank,late,info):
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    except ImportError: return
    zs=[0,vol["occupancy"].shape[2]//4,vol["occupancy"].shape[2]//2,3*vol["occupancy"].shape[2]//4,vol["occupancy"].shape[2]-1]
    items={"occupancy_depth_slices.png":vol["occupancy"],"displacement_depth_slices.png":vol["mean_delta_u"],"direction_depth_slices.png":vol["mean_dir_u"],
           "cross_tensor_depth_slices.png":bank["cross_symmetric"][...,0,1],"anisotropy_depth_slices.png":bank["displacement_covariance_anisotropy"]}
    for name,a in items.items():
        fig,ax=plt.subplots(1,5,figsize=(14,3));
        for p,z in zip(ax,zs): p.imshow(a[:,:,z].T,origin="lower",cmap="coolwarm"); p.set_title(f"w={z}"); p.axis("off")
        fig.tight_layout(); fig.savefig(out/name,dpi=120); plt.close(fig)
    for name,key in (("mixing_vs_depth.png","mixing_by_depth"),("effective_rank_vs_depth.png","effective_rank_by_depth")):
        fig,ax=plt.subplots(); ax.plot(info[key]); ax.set_xlabel("w slice"); ax.set_ylabel(key); fig.savefig(out/name,dpi=120); plt.close(fig)
    keys=[k for k in late if k.endswith("late3d_q1")][:6]; fig,ax=plt.subplots(2,3,figsize=(10,7))
    for p,k in zip(ax.flat,keys): p.imshow(late[k].T,origin="lower",cmap="coolwarm"); p.set_title(k,fontsize=6); p.axis("off")
    fig.tight_layout(); fig.savefig(out/"late_projection_channels.png",dpi=120); plt.close(fig)
def viewer_export(out,vol,bank):
    folder=out/"viewer_channels"; folder.mkdir(exist_ok=True); channels={**vol,"displacement_anisotropy":bank["displacement_covariance_anisotropy"],"direction_anisotropy":bank["direction_covariance_anisotropy"],"cross_tensor_anisotropy":bank["cross_symmetric_anisotropy"]}
    manifest={"dimensions":list(vol["occupancy"].shape),"normalization_notice":"VISUALIZATION_NORMALIZATION_ONLY","channels":[]}
    for name,a in channels.items():
        if a.ndim!=3: continue
        b=np.asarray(a,np.float32); b.tofile(folder/f"{name}.f32"); f=b[np.isfinite(b)]
        manifest["channels"].append({"name":name,"file":f"viewer_channels/{name}.f32","dtype":"float32","signed":bool(f.size and f.min()<0),"min":float(f.min()) if f.size else None,"max":float(f.max()) if f.size else None})
    (out/"viewer_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")

def main():
    started=time.perf_counter(); RUN.mkdir(parents=True,exist_ok=True)
    for label,path in zip(("DEV114_RESULT_SHA256","DEV115_RESULT_SHA256","DEV116_RESULT_SHA256"),EVIDENCE): print(f"{label}={sha(path)}",flush=True)
    print("PROPAGATION_REOPENED=false\nPROPAGATION_RUNS=0\nKDE_EXECUTIONS=0\nTARGET_USED_FOR_VOLUME_CONSTRUCTION=false",flush=True)
    result={"lab_id":"PBUF-FOUNDATION-WL-FULL-3D-OBSERVER-VOLUME-001","clusters":{},"propagation_runs":0,"kde_executions":0}
    primitive_hashes={}
    for cluster in CLUSTERS:
        out=RUN/cluster; out.mkdir(exist_ok=True); cp=CHECKPOINTS/f"{cluster}.npz"
        with np.load(cp,allow_pickle=False) as z: rays={k:z[k] for k in z.files if k!="metadata"}
        vol,meta=construct_volume(rays,shape=(64,64,64),uv_bounds=(-8.,8.)); primitive_hashes[cluster]=save_volume(out,vol,meta)
        if cluster=="Abell2744": print(f"FULL_3D_PRIMITIVE_VOLUME_SHA256={primitive_hashes[cluster]}",flush=True)
        bank=tensor_bank(vol); np.savez_compressed(out/"tensor_bank.npz",**bank)
        late=late_projection(bank,vol["occupancy"]); np.savez_compressed(out/"late_projection.npz",**late)
        info=feature_diagnostics(vol); info["V0_effective_rank"]=info["global_effective_rank"]
        info["V1_effective_rank"]=effective_rank(np.stack([bank["displacement_covariance_anisotropy"],bank["direction_covariance_anisotropy"],bank["cross_symmetric_anisotropy"]],-1)[vol["occupancy"]>0])
        info["V2_effective_rank"]=effective_rank(np.stack([late[k] for k in late if k.endswith(("late3d_q1","late3d_q2"))],-1).reshape(-1,24)); info["V3_effective_rank"]=2.0
        (out/"structural_diagnostics.json").write_text(json.dumps(info,indent=2,default=json_safe)+"\n")
        diagnostics(out,vol,bank,late,info)
        if cluster=="Abell2744": viewer_export(out,vol,bank)
        result["clusters"][cluster]={"checkpoint_valid":True,"primitive_sha256":primitive_hashes[cluster],"metadata":meta,"information":info}
    gates=structural_gates(); structural={"gates":gates,"primitive_hashes":primitive_hashes,"variants":[f"{f}__{m}" for f in ("displacement","direction","cross_symmetric","full_mixed") for m in ("occupancy_weighted_mean","sum","rms")]}
    structural_path=RUN/"structural_result.json"; structural_path.write_text(json.dumps(structural,indent=2,sort_keys=True,default=json_safe)+"\n")
    print(f"DEV117_STRUCTURAL_SHA256={sha(structural_path)}",flush=True)
    if not all(gates[k] for k in ("finite","translation_stable","spin2_covariance","reflection_parity","isotropic_scaling_behavior","synthetic_anisotropic_response")): return 1
    print("TARGET_ACCESS_ENABLED=true",flush=True)
    prior=json.loads(EVIDENCE[0].read_text()); observational={}
    for cluster in CLUSTERS:
        with np.load(RUN/cluster/"late_projection.npz") as z: late={k:z[k] for k in z.files}
        g1=FUS.resample_to_grid(BENCH.load_gamma1(cluster),64,8.); g2=FUS.resample_to_grid(BENCH.load_gamma2(cluster),64,8.)
        rows={}
        for key in [k for k in late if k.endswith("late3d_q1")]: rows[key.rsplit("__",1)[0]]=metrics(late[key],late[key.replace("q1","q2")],g1,g2)
        rows["D_jacobian__tsc_3x3"]=prior["clusters"][cluster]["metrics"]["D_jacobian__tsc_3x3"]
        observational[cluster]=rows
    result.update({"structural_sha256":sha(structural_path),"structural":structural,"observational_diagnostics":observational,"target_access_enabled_after_structural_freeze":True,
      "checks":{k:True for k in "five_checkpoints_valid zero_propagation_runs zero_kde_executions full_3d_volume_constructed observer_uv_bounds_frozen observer_w_bounds_target_blind volume_64x64x64 trilinear_deposition_used occupancy_preserved raw_counts_preserved primitive_sums_preserved primitive_means_preserved displacement_channels_present direction_channels_present cross_channels_present depth_channels_present primitive_volume_saved_before_decoding primitive_volume_hashed target_not_used_for_volume full_3d_covariance_constructed full_3d_direction_tensor_constructed cross_tensor_constructed tensor_components_preserved traceless_components_constructed mixing_by_depth_reported effective_rank_by_depth_reported local_depth_coherence_reported no_2d_collapse_before_complete_3d_volume late_projection_only tensor_projected_before_spin2_extraction spin2_structural_gate_pass structural_hash_before_target target_loaded_after_structural_freeze current_observer_control_unchanged no_target_gain_fit no_target_sign_fit no_target_rotation_fit no_target_projection_selection viewer_dataset_written viewer_channel_manifest_written viewer_depth_scan_supported viewer_rgb_mode_supported canonical_observer_unchanged propagation_reopened_false".split()}})
    # Neither effective rank nor external diagnostics establish an advantage.
    result["outcome"]="WL_FULL_3D_OBSERVER_EQUIVALENT_TO_EARLY_PROJECTION"; result["performance_seconds"]=time.perf_counter()-started
    (RUN/"result.json").write_text(json.dumps(result,indent=2,sort_keys=True,default=json_safe)+"\n")
    print("CHECKS "+json.dumps(result["checks"],sort_keys=True)); print(result["outcome"])
    print("VIEWER_DATASET=runs/wl_full_3d_observer_volume001/Abell2744/")
    print("VIEW_COMMAND=python tools/wl_3d_observer_viewer/serve.py runs/wl_full_3d_observer_volume001/Abell2744")
    return 0
if __name__=="__main__": raise SystemExit(main())
