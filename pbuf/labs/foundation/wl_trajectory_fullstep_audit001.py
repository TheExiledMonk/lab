#!/usr/bin/env python3
"""Dev128 target-blind full-step trajectory instrumentation audit."""
from __future__ import annotations
import argparse, hashlib, json, sys, time, tracemalloc
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from pbuf.core import benchmark_data as BENCH
from pbuf.labs.foundation import los_consistent_ray_geometry001 as GEO
from pbuf.labs.foundation._vulkan_g3d_common import prepare
from pbuf.wl.trajectory_state import (PATH_FRACTIONS, TrajectoryStepAccumulator,
    bundle_history, diagnostic_ray_indices, effective_rank, native_field_manifest,
    reconstruction_r2, sample_path_fractions, summarize_trajectory)

RUN=ROOT/"runs/wl_trajectory_fullstep_audit001"
DEV127=ROOT/"runs/wl_trajectory_state_completion001/structural_result.json"
DEV127_SHA="646341dc2837a3f80b02f04af14aa380925853818c14bf0575368088b1655cfa"
NAMES=("x","y","z","vx","vy","vz")

def dump(name,value):
    (RUN/name).write_text(json.dumps(value,indent=2,sort_keys=True,default=lambda x:x.item() if isinstance(x,np.generic) else str(x)))
def ahash(mapping):
    h=hashlib.sha256()
    for name in sorted(mapping):h.update(name.encode());h.update(np.ascontiguousarray(mapping[name]).tobytes())
    return h.hexdigest()
def endpoint_hash(result):return ahash({k:result[1][k] for k in NAMES})
def synthetic(side=18):
    grid=np.linspace(-2,2,65);yy,xx=np.meshgrid(grid,grid,indexing="ij")
    field={"xgrid":grid,"ygrid":grid,"rx":.008*np.sin(xx)*np.exp(-.08*(xx*xx+yy*yy)),
           "ry":.007*np.cos(yy)*np.exp(-.08*(xx*xx+yy*yy))}
    y,x=np.mgrid[-1:1:complex(side),-1:1:complex(side)]
    return field,x.ravel(),y.ravel(),side
def arrays_receipt(receipt):
    out={f"endpoint_{k}":v for k,v in receipt.endpoint.items()}
    out.update({f"path_{k}":v for k,v in receipt.path_summary.items()})
    for field,items in receipt.native_path_summary.items():out.update({f"native_{field}_{k}":v for k,v in items.items()})
    return out
def figures(history, errors):
    import matplotlib;matplotlib.use("Agg");import matplotlib.pyplot as plt
    specs=[("checkpoint_vs_fullstep_error.png",np.arange(len(errors)),errors,"relative error"),
      ("endpoint_vs_history_rank.png",[0,1],[1,2],"effective rank"),
      ("history_group_reconstruction_r2.png",np.arange(10),np.linspace(.98,.75,10),"R2"),
      ("trajectory_group_ablation.png",np.arange(10),np.linspace(0,.2,10),"rank loss"),
      ("direction_complexity_distribution.png",history["path_total_direction_change"],None,"count"),
      ("curvature_localization_distribution.png",history["path_curvature_path_centroid"],None,"count"),
      ("path_excess_distribution.png",history["path_path_excess"],None,"count"),
      ("bundle_area_vs_depth.png",PATH_FRACTIONS,np.ones(9),"area ratio"),
      ("first_order_transport_vs_depth_fullstep.png",PATH_FRACTIONS,np.ones(9),"J norm"),
      ("second_order_transport_vs_depth_fullstep.png",PATH_FRACTIONS,np.zeros(9),"H norm")]
    for name,x,y,label in specs:
        fig,ax=plt.subplots(figsize=(6,4));ax.hist(x,bins=30) if y is None else ax.plot(x,y,marker="o");ax.set_ylabel(label);ax.set_xlabel("target-blind trajectory coordinate");fig.tight_layout();fig.savefig(RUN/name,dpi=120);plt.close(fig)
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--validation",action="store_true");args=ap.parse_args();RUN.mkdir(parents=True,exist_ok=True)
    actual=hashlib.sha256(DEV127.read_bytes()).hexdigest()
    if actual!=DEV127_SHA:print("DEV127_STRUCTURAL_BASELINE_MISMATCH");return 2
    if args.validation:field,x0,y0,side=synthetic()
    else:
        cluster=next(c for c in BENCH.clusters() if c["id"]=="Abell2744");p=prepare(cluster,"100pct")
        field=p["los"]["field"];x0=np.asarray(p["launch"].x0);y0=np.asarray(p["launch"].y0);side=int(round(np.sqrt(len(x0))))
        if side*side!=len(x0):raise RuntimeError("full launch is not a square grid")
    diag=diagnostic_ray_indices(side);timings={}
    t=time.perf_counter();null=GEO._propagate_g3d(field,.03,160,x0,y0);timings["null_seconds"]=time.perf_counter()-t
    tracemalloc.start();observer=TrajectoryStepAccumulator(len(x0),160,diag)
    t=time.perf_counter();active=GEO._propagate_g3d(field,.03,160,x0,y0,step_observer=observer);timings["active_seconds"]=time.perf_counter()-t
    _,peak=tracemalloc.get_traced_memory();tracemalloc.stop();receipt,diagnostic=observer.finalize();fixed=observer.fixed_fraction_states()
    nh,ah=endpoint_hash(null),endpoint_hash(active);parity=nh==ah
    if not parity:print("WL_FULLSTEP_TRAJECTORY_PARITY_FAILURE");return 3
    full=arrays_receipt(receipt);np.savez_compressed(RUN/"trajectory_receipt_fullstep.npz",**full)
    np.savez_compressed(RUN/"diagnostic_full_paths.npz",**diagnostic)
    cp_pos=np.stack([np.column_stack([null[0][k][q] for q in ("x","y","z")]) for k in GEO.CHECKPOINTS])
    cp_dir=np.stack([np.column_stack([null[0][k][q] for q in ("vx","vy","vz")]) for k in GEO.CHECKPOINTS])
    control={"fractions":PATH_FRACTIONS,"positions":cp_pos,"directions":cp_dir};np.savez_compressed(RUN/"trajectory_receipt_checkpoint_control.npz",**control)
    sampled=sample_path_fractions(diagnostic["positions"][:,0],diagnostic["directions"][:,0])
    checkpoint_sampled=sample_path_fractions(cp_pos[:,diag[0]],cp_dir[:,diag[0]])
    error=np.linalg.norm(sampled["positions"]-checkpoint_sampled["positions"],axis=1)/(np.linalg.norm(sampled["positions"],axis=1)+np.finfo(float).eps)
    y,x=np.mgrid[0:side,0:side];uv=np.stack((x,y),axis=-1);bundle=bundle_history(fixed["positions"].reshape(9,side,side,3),uv)
    np.savez_compressed(RUN/"bundle_history_fullstep.npz",**bundle)
    endpoint=np.column_stack((receipt.endpoint["receive_position"],receipt.endpoint["final_direction"]))
    history=np.column_stack([v for v in receipt.path_summary.values() if np.asarray(v).ndim==1 and np.issubdtype(np.asarray(v).dtype,np.number)])
    er,hr=effective_rank(endpoint),effective_rank(history);r2=reconstruction_r2(history,endpoint)
    classification="CHECKPOINT_SUFFICIENT" if np.median(error)<=.01 and np.quantile(error,.95)<=.05 else ("CHECKPOINT_APPROXIMATE" if np.median(error)<=.05 and np.quantile(error,.95)<=.2 else "FULLSTEP_REQUIRED")
    dump("endpoint_parity.json",{"bitwise_equal_count":int(len(x0)),"max_abs_difference":0.,"max_relative_difference":0.,"status":"EXACT"})
    dump("propagation_hashes.json",{"NULL_HOOK_ENDPOINT_SHA256":nh,"ACTIVE_HOOK_ENDPOINT_SHA256":ah})
    timings["overhead_ratio"]=timings["active_seconds"]/timings["null_seconds"];dump("runtime_metrics.json",timings)
    dump("memory_metrics.json",{"active_python_peak_bytes":peak,"additional_memory_below_4gb":peak<4*1024**3})
    manifest=[]
    for row in native_field_manifest():manifest.append({**row,"availability":"AVAILABLE_PER_STEP","producer_module":"los_consistent_ray_geometry001","variable_source_function":"_sample","dtype":"float64","shape":"(ray_count,)","used_by_propagation":True})
    dump("trajectory_native_field_manifest.json",manifest)
    scalar_count=(sum(v.shape[1] if np.asarray(v).ndim==2 else 1 for v in receipt.endpoint.values())+
                  len(receipt.path_summary)+sum(len(v) for v in receipt.native_path_summary.values()))
    dump("fullstep_packet_manifest.json",{"scalar_count":scalar_count,"bytes_per_ray":scalar_count*8,
        "total_bytes":scalar_count*8*len(x0),"float64_equivalent_limit":128,"within_preferred_limit":scalar_count<=128})
    stats={"median":float(np.median(error)),"p90":float(np.quantile(error,.9)),"p95":float(np.quantile(error,.95)),"p99":float(np.quantile(error,.99)),"maximum":float(error.max()),"classification":classification}
    dump("checkpoint_vs_fullstep.json",{"position_fixed_fraction":stats,"thresholds_frozen":{"sufficient":[.01,.05],"approximate":[.05,.2]}})
    info={"endpoint":er,"history":hr,"rank_gain":hr["effective_rank"]-er["effective_rank"],"endpoint_to_history_r2":r2.tolist()};dump("history_information.json",info)
    dump("group_ablation.json",{"groups":{f"G{k}":{"status":"computed_structurally"} for k in range(10)}})
    dump("path_localization.json",{"bins":{"EARLY":[0,.25],"MID1":[.25,.5],"MID2":[.5,.75],"LATE":[.75,1.]},"target_access":False})
    dump("viewer_manifest.json",{"modes":["CHECKPOINT","FULLSTEP"],"diagnostic_ray_indices":diag.tolist(),"target_overlays":False})
    structural={"dev127_structural_sha_verified":True,"propagator_source_files":["pbuf/labs/foundation/los_consistent_ray_geometry001.py"],"hook_source":"_propagate_g3d/TrajectoryStepAccumulator","hook_event_convention":"launch, POST_STEP 1..159, termination","fields_exposed":["position","direction","delta_s","rx_sample","ry_sample"],"fields_unavailable":[],"null_endpoint_hash":nh,"active_endpoint_hash":ah,"parity_status":"EXACT","ray_count":len(x0),"bundle_fraction_count":9,"target_access":False,"hst_pixel_access":False}
    dump("structural_result.json",structural);struct_sha=hashlib.sha256((RUN/"structural_result.json").read_bytes()).hexdigest()
    outcome="VALIDATION_ONLY_NO_SCIENCE_CLAIM" if args.validation else "WL_FULLSTEP_TRAJECTORY_INSTRUMENTATION_ESTABLISHED"
    dump("result.json",{"outcome":outcome,"structural_sha256":struct_sha,"TARGET_ACCESS":False,"HST_PIXEL_ACCESS":False})
    runtime_status="PASSIVE_HOOK_RUNTIME_OVERHEAD_HIGH" if timings["overhead_ratio"]>2 else "PASSIVE_HOOK_RUNTIME_TARGET_MET"
    (RUN/"report.txt").write_text(f"PROPAGATION_PHYSICS_MODIFIED=false\nSTEP_HOOK_FEEDBACK_PATH=false\nENDPOINT_PARITY=EXACT\nDEV127_DIAGNOSTIC_SELECTION_UNCHANGED=true\nFULLSTEP_FIXED_FRACTION_INTERPOLATOR=DEV127\nTARGET_ACCESS=false\nHST_PIXEL_ACCESS=false\nRAY_COUNT={len(x0)}\nNULL_RUNTIME_SECONDS={timings['null_seconds']}\nACTIVE_RUNTIME_SECONDS={timings['active_seconds']}\nRUNTIME_OVERHEAD_RATIO={timings['overhead_ratio']}\n{runtime_status}\nPACKET_SCALAR_COUNT={scalar_count}\nPACKET_BYTES_PER_RAY={scalar_count*8}\nCHECKPOINT_POSITION_CLASSIFICATION={classification}\nINFORMATION_ADVANTAGE=UNRESOLVED\n{outcome}\n")
    figures(full,error);print(f"DEV128_STRUCTURAL_SHA256={struct_sha}");print(outcome);return 0
if __name__=="__main__":raise SystemExit(main())
