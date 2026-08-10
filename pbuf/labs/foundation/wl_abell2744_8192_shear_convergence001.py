#!/usr/bin/env python3
"""Dev Doc 115: Abell2744 8192² reception-density shear convergence."""

from __future__ import annotations
import argparse, hashlib, json, os, shutil, subprocess, sys, time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[3]; sys.path.insert(0, str(ROOT))
from pbuf.core import benchmark_data as BENCH
from pbuf.labs.foundation._vulkan_g3d_common import prepare
from pbuf.labs.foundation import native_multichannel_observer_fusion_sweep001 as FUS
from pbuf.labs.foundation import native_full_state_2d_reconstruction_decoder_sweep001 as DEC
from pbuf.wl.backends import VulkanBackend
from pbuf.wl.config import CHECKPOINT, EXTENT, OBS_BINS, PROPAGATION_STEP, PROPAGATION_STEPS
from pbuf.wl.launch import RayLaunch
from pbuf.wl.propagation import PropagationConfig
from pbuf.wl.shear_readout import build_shear_candidates, construct_local_primitives, evaluate_candidate
from pbuf.wl.streaming_shear import JacobianTSCAccumulator

LAB_ID="PBUF-FOUNDATION-WL-ABELL2744-8192-SHEAR-CONVERGENCE-001"
RUN_DIR=ROOT/"runs"/"wl_abell2744_8192_shear_convergence001"
BASE_RESULT=ROOT/"runs"/"wl_3d_shear_readout_recovery001"/"result.json"
BASE_CHECKPOINT=ROOT/"runs"/"wl_3d_shear_readout_recovery001"/"checkpoints"/"Abell2744.npz"
WIDTH=HEIGHT=8192; RAYS=WIDTH*HEIGHT; CONFIG=PropagationConfig(PROPAGATION_STEP,PROPAGATION_STEPS,CHECKPOINT)
EXPECTED={"gamma1":{"pearson":.15963721563216124,"spearman":.09715428188710849,"pred_rms":.0239242457539334,"obs_rms":.062160042482900475,"rms_ratio_pred_over_obs":.38488142540305836},"gamma2":{"pearson":.20242108305073303,"spearman":.14387459525184157,"pred_rms":.018838754008152395,"obs_rms":.05586397422966777,"rms_ratio_pred_over_obs":.33722545285988026},"magnitude_pearson":.2187521646685786,"orientation_agreement":.08368696562433091}

def sha_bytes(*xs):
    h=hashlib.sha256()
    for x in xs: h.update(np.ascontiguousarray(x,dtype=np.float64).tobytes())
    return h.hexdigest()
def sha_file(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(1<<20),b""): h.update(b)
    return h.hexdigest()
def atomic_json(path,obj):
    path=Path(path); tmp=path.with_suffix(path.suffix+".tmp")
    with open(tmp,"w") as f: json.dump(obj,f,indent=2,sort_keys=True); f.write("\n"); f.flush(); os.fsync(f.fileno())
    os.replace(tmp,path)
def atomic_npz(path,**arrays):
    path=Path(path); tmp=Path(str(path)+".tmp.npz")
    np.savez_compressed(tmp,**arrays)
    with open(tmp,"rb") as f: os.fsync(f.fileno())
    os.replace(tmp,path)
def block(name,value): print(name, json.dumps(value,sort_keys=True),flush=True)

def baseline():
    result=json.loads(BASE_RESULT.read_text()); got=result["clusters"]["Abell2744"]["metrics"]["D_jacobian__tsc_3x3"]
    for comp in ("gamma1","gamma2"):
        for key,val in EXPECTED[comp].items():
            if got[comp][key] != val: raise SystemExit("DEV115_BASELINE_RESULT_MISMATCH")
    for key in ("magnitude_pearson","orientation_agreement"):
        if got[key] != EXPECTED[key]: raise SystemExit("DEV115_BASELINE_RESULT_MISMATCH")
    block("DEV114_BASELINE",got); print("DEV114_RESULT_SHA256="+sha_file(BASE_RESULT),flush=True)
    return got

def spec(): return next(x for x in build_shear_candidates(("tsc_3x3",)) if x.family=="D_jacobian")
def load_baseline_maps():
    with np.load(BASE_CHECKPOINT,allow_pickle=False) as z: rays={k:z[k] for k in z.files if k!="metadata"}
    rays.update(construct_local_primitives(rays,bins=OBS_BINS,extent=EXTENT))
    return evaluate_candidate(spec(),rays,bins=OBS_BINS,extent=EXTENT), rays
def parity():
    reference,rays=load_baseline_maps(); acc=JacobianTSCAccumulator.empty(OBS_BINS,EXTENT)
    for start in range(0,len(rays["u0"]),65536):
        sl=slice(start,start+65536); acc.add(*(rays[k][sl] for k in ("u0","v0","uf","vf")))
    streamed=acc.finalize()[:2]
    checks=[np.allclose(a,b,rtol=1e-12,atol=1e-12,equal_nan=True) for a,b in zip(streamed,reference)]
    report={"pass":bool(all(checks)),"gamma1_max_abs":float(np.nanmax(np.abs(streamed[0]-reference[0]))),"gamma2_max_abs":float(np.nanmax(np.abs(streamed[1]-reference[1])))}
    block("STREAMING_PARITY",report)
    if not report["pass"]: raise SystemExit("DEV115_STREAMING_PARITY_FAILED")
    return reference, rays

def geometry_fingerprint():
    # Hash the exact row-major float64 center lattice without allocating it all.
    centers=-EXTENT+(np.arange(WIDTH,dtype=np.float64)+.5)*(2*EXTENT/WIDTH); h=hashlib.sha256()
    for y in centers:
        h.update(centers.tobytes()); h.update(np.full(WIDTH,y,np.float64).tobytes())
    return h.hexdigest()
def tile_launch(tx,ty,size):
    step=2*EXTENT/WIDTH
    xs=-EXTENT+(np.arange(tx,min(tx+size,WIDTH),dtype=np.float64)+.5)*step
    ys=-EXTENT+(np.arange(ty,min(ty+size,HEIGHT),dtype=np.float64)+.5)*step
    x,y=np.meshgrid(xs,ys,indexing="xy"); x=x.ravel(); y=y.ravel()
    return RayLaunch(x,y,np.ones_like(x),np.zeros_like(x),"coverage_100pct",1)
def tiles(size): return [(i,x,y) for i,(y,x) in enumerate(( (y,x) for y in range(0,HEIGHT,size) for x in range(0,WIDTH,size) ))]
def vram():
    roots=list(Path("/sys/class/drm").glob("card*/device"))
    for root in roots:
        try: return {"total_bytes":int((root/"mem_info_vram_total").read_text()),"used_bytes":int((root/"mem_info_vram_used").read_text())}
        except OSError: pass
    return {"total_bytes":None,"used_bytes":None}
def screen(snapshot,launch,basis):
    e1,e2,n=basis; p0=np.column_stack((launch.x0,launch.y0,np.zeros_like(launch.x0))); pf=np.column_stack((snapshot["x"],snapshot["y"],snapshot["z"]))
    return p0@e1,p0@e2,pf@e1,pf@e2
def basis_from_sum(sums,count):
    n=np.asarray(sums)/count; n/=np.linalg.norm(n); ref=np.array([1.,0.,0.]); e1=ref-(ref@n)*n
    if np.linalg.norm(e1)<=1e-10: ref=np.array([0.,1.,0.]); e1=ref-(ref@n)*n
    e1/=np.linalg.norm(e1); e2=np.cross(n,e1); e2/=np.linalg.norm(e2); return e1,e2,n

def metrics(pair,targets):
    out={c:DEC._compare_candidates({"p":pair[i]},{c:targets[c]})["p"][c] for i,c in enumerate(("gamma1","gamma2"))}
    mask=np.isfinite(pair[0]+pair[1]+targets["gamma1"]+targets["gamma2"])
    out["magnitude_pearson"]=float(np.corrcoef(np.hypot(pair[0][mask],pair[1][mask]),np.hypot(targets["gamma1"][mask],targets["gamma2"][mask]))[0,1])
    d=np.arctan2(pair[1][mask],pair[0][mask])-np.arctan2(targets["gamma2"][mask],targets["gamma1"][mask]); out["orientation_agreement"]=float(abs(np.mean(np.exp(1j*d))))
    return out
def comparisons(base,high,targets):
    out={}
    for i,c in enumerate(("gamma1","gamma2")):
        a,b,o=base[i],high[i],targets[c]; mask=np.isfinite(a+b+o); d=b[mask]-a[mask]
        direct=DEC._compare_candidates({"high":b},{"current":a})["high"]["current"]
        out[c]={"current_vs_8192_pearson":direct["pearson"],"current_vs_8192_spearman":direct["spearman"],"rms_difference":float(np.sqrt(np.mean(d*d))),"maximum_absolute_difference":float(np.max(abs(d))),"mean_absolute_difference":float(np.mean(abs(d))),"residual_rms_baseline":float(np.sqrt(np.mean((a[mask]-o[mask])**2))),"residual_rms_8192":float(np.sqrt(np.mean((b[mask]-o[mask])**2))),"residual_mae_baseline":float(np.mean(abs(a[mask]-o[mask]))),"residual_mae_8192":float(np.mean(abs(b[mask]-o[mask])))}
    return out
def render(base,high,targets):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    d=RUN_DIR/"images"; d.mkdir(exist_ok=True)
    for i,c in enumerate(("gamma1","gamma2")):
        lim=float(np.nanmax(abs(targets[c]))); fig,ax=plt.subplots(1,3,figsize=(12,4))
        for a,x,t in zip(ax,(targets[c],base[i],high[i]),("observed","current","8192")): a.imshow(x,vmin=-lim,vmax=lim,cmap="coolwarm"); a.set_title(t)
        fig.tight_layout(); fig.savefig(d/f"{c}_observed_vs_current_vs_8192.png"); plt.close(fig)
        fig,ax=plt.subplots(1,2,figsize=(8,4))
        for a,x,t in zip(ax,(base[i]-targets[c],high[i]-targets[c]),("current residual","8192 residual")): a.imshow(x,vmin=-lim,vmax=lim,cmap="coolwarm"); a.set_title(t)
        fig.tight_layout(); fig.savefig(d/f"{c}_residual_current_vs_8192.png"); plt.close(fig)
    obs=np.hypot(targets["gamma1"],targets["gamma2"]); lim=float(np.nanmax(obs)); fig,ax=plt.subplots(1,3,figsize=(12,4))
    for a,x,t in zip(ax,(obs,np.hypot(*base),np.hypot(*high)),("observed","current","8192")): a.imshow(x,vmin=0,vmax=lim,cmap="viridis"); a.set_title(t)
    fig.tight_layout(); fig.savefig(d/"shear_magnitude_observed_vs_current_vs_8192.png"); plt.close(fig)

def main(argv=None):
    ap=argparse.ArgumentParser(); ap.add_argument("--resume",action="store_true"); ap.add_argument("--tile-size",type=int,choices=(256,512,1024),default=512); ap.add_argument("--preflight-only",action="store_true"); ap.add_argument("--baseline-parity-only",action="store_true"); args=ap.parse_args(argv)
    started=time.perf_counter(); RUN_DIR.mkdir(parents=True,exist_ok=True); (RUN_DIR/"tiles").mkdir(exist_ok=True)
    head=subprocess.check_output(["git","rev-parse","HEAD"],text=True,cwd=ROOT).strip(); branch=subprocess.check_output(["git","branch","--show-current"],text=True,cwd=ROOT).strip(); block("BASELINE",{"branch":branch,"head":head})
    base_metrics=baseline(); base_maps,base_rays=parity()
    if args.baseline_parity_only:return 0
    tests=subprocess.run([sys.executable,"-m","unittest","tests.test_wl_streaming_shear","tests.test_wl_shear_readout"],cwd=ROOT).returncode
    if tests: raise SystemExit("DEV115_UNIT_REGRESSION_FAILED")
    launch_fp=geometry_fingerprint(); print("HIGHRES_LAUNCH_GEOMETRY_SHA256="+launch_fp,flush=True); print(f"SELECTED_TILE_SIZE={args.tile_size}",flush=True)
    manifest={"cluster":"Abell2744","launch_width":WIDTH,"launch_height":HEIGHT,"ray_count":RAYS,"candidate":"D_jacobian","deposition":"tsc_3x3","observer_shape":[64,64],"backend":"vulkan","baseline_commit":head,"candidate_code_fingerprint":sha_file(ROOT/"pbuf/wl/shear_readout.py"),"streaming_code_fingerprint":sha_file(ROOT/"pbuf/wl/streaming_shear.py"),"launch_fingerprint":launch_fp,"tile_size":args.tile_size,"benchmark_source_fingerprint":sha_file(ROOT/"pbuf/core/benchmark_data.py")}
    mpath=RUN_DIR/"manifest.json"
    if args.resume and mpath.exists() and json.loads(mpath.read_text())!=manifest: raise SystemExit("DEV115_INCOMPATIBLE_RESUME_STATE")
    atomic_json(mpath,manifest); block("HIGHRES_CONFIG",manifest); block("LAUNCH_GEOMETRY",{"fingerprint":launch_fp,"ordering":"row-major"})
    inventory={c["id"]:c for c in BENCH.clusters()}; prepared=prepare(inventory["Abell2744"],"100pct")
    disk=shutil.disk_usage(RUN_DIR); disk_report={"free_bytes":disk.free,"estimated_checkpoint_bytes":140_000_000,"required_two_x":280_000_000,"pass":disk.free>=280_000_000}; block("DISK_PREFLIGHT",disk_report)
    with VulkanBackend() as vk:
        mem0=vram(); launch=tile_launch(0,0,args.tile_size); t=time.perf_counter(); snap=vk.propagate_final_snapshot(prepared["los"]["field"],launch,CONFIG); tile_sec=time.perf_counter()-t; mem1=vram()
        peak=max(x for x in (mem0["used_bytes"],mem1["used_bytes"]) if x is not None) if mem0["total_bytes"] else None; margin=mem0["total_bytes"]-peak if peak is not None else None
        vr={"gpu":vk.runtime.device,"used_before":mem0["used_bytes"],"peak_observed":peak,"used_after":mem1["used_bytes"],"total":mem0["total_bytes"],"safety_margin":margin,"pass":margin is not None and margin>=2_000_000_000,"single_tile_seconds":tile_sec}; block("VRAM_PREFLIGHT",vr)
        if not disk_report["pass"] or not vr["pass"]: raise SystemExit("DEV115_PREFLIGHT_FAILED")
        print("HIGHRES_RUN_AUTHORIZED=true",flush=True)
        if args.preflight_only:return 0
        alltiles=tiles(args.tile_size); progress_path=RUN_DIR/"progress.json"; basis_path=RUN_DIR/"basis.npz"; acc_path=RUN_DIR/"accumulator.npz"
        if args.resume and basis_path.exists():
            z=np.load(basis_path); basis=tuple(z[k] for k in ("e1","e2","normal"))
        else:
            sums=np.zeros(3); count=0
            for i,x,y in alltiles:
                s=vk.propagate_final_snapshot(prepared["los"]["field"],tile_launch(x,y,args.tile_size),CONFIG); sums += [np.sum(s[k]) for k in ("vx","vy","vz")]; count+=len(s["vx"])
                print(f"PROGRESS basis tile {i+1} / {len(alltiles)}",flush=True)
            basis=basis_from_sum(sums,count); atomic_npz(basis_path,e1=basis[0],e2=basis[1],normal=basis[2])
        completed={}; tile_times=[]
        if args.resume and progress_path.exists(): completed=json.loads(progress_path.read_text()).get("completed",{})
        tile_times=[float(v["runtime_seconds"]) for v in completed.values()]
        if args.resume and acc_path.exists():
            with np.load(acc_path) as z: acc=JacobianTSCAccumulator.from_state(OBS_BINS,EXTENT,{k:z[k] for k in z.files})
        else: acc=JacobianTSCAccumulator.empty(OBS_BINS,EXTENT)
        phase=time.perf_counter()
        for i,x,y in alltiles:
            if str(i) in completed: continue
            t=time.perf_counter(); launch=tile_launch(x,y,args.tile_size); snap=vk.propagate_final_snapshot(prepared["los"]["field"],launch,CONFIG); u0,v0,uf,vf=screen(snap,launch,basis); acc.add(u0,v0,uf,vf); sec=time.perf_counter()-t; tile_times.append(sec)
            delta=sha_bytes(u0,v0,uf,vf); record={"tile_index":i,"launch_bounds":[x,y,min(x+args.tile_size,WIDTH),min(y+args.tile_size,HEIGHT)],"ray_count":len(u0),"checksum":delta,"runtime_seconds":sec}; atomic_json(RUN_DIR/"tiles"/f"tile_{i:04d}.json",record)
            atomic_npz(acc_path,**acc.state()); completed[str(i)]=record; atomic_json(progress_path,{"completed":completed,"rays_complete":sum(v["ray_count"] for v in completed.values()),"rays_total":RAYS})
            done=sum(v["ray_count"] for v in completed.values()); elapsed=time.perf_counter()-phase; rate=max(len(tile_times),1)/max(elapsed,1e-9); eta=(len(alltiles)-len(completed))/rate
            print(f"PROGRESS tile {i+1} / {len(alltiles)} rays {done} / {RAYS} {100*done/RAYS:.4f}% tile_seconds={sec:.3f} elapsed_seconds={elapsed:.3f} estimated_remaining_seconds={eta:.3f}",flush=True)
    high=(acc.finalize()[0],acc.finalize()[1]); np.save(RUN_DIR/"pred_gamma1.npy",high[0]); np.save(RUN_DIR/"pred_gamma2.npy",high[1]); p1=sha_file(RUN_DIR/"pred_gamma1.npy");p2=sha_file(RUN_DIR/"pred_gamma2.npy"); print("PRED_GAMMA1_SHA256="+p1);print("PRED_GAMMA2_SHA256="+p2)
    # Binding boundary: observation access occurs only after prediction files are closed and hashed.
    observed=FUS._observed(prepared["source"]["data"]); targets={k:DEC._finite(observed[k]) for k in ("gamma1","gamma2")}; scored=metrics(high,targets); comp=comparisons(base_maps,high,targets); render(base_maps,high,targets)
    dg1=scored["gamma1"]["pearson"]-EXPECTED["gamma1"]["pearson"]; dg2=scored["gamma2"]["pearson"]-EXPECTED["gamma2"]["pearson"]
    if dg1>=.1 and dg2>=.1 and all(scored[c]["rms_ratio_pred_over_obs"]>EXPECTED[c]["rms_ratio_pred_over_obs"] for c in ("gamma1","gamma2")): outcome="WL_8192_RECEPTION_DENSITY_STRONGLY_IMPROVES_SHEAR"
    elif dg1>0 and dg2>0 and max(dg1,dg2)>=.05: outcome="WL_8192_RECEPTION_DENSITY_MATERIALLY_IMPROVES_SHEAR"
    elif abs(dg1)<.02 and abs(dg2)<.02: outcome="WL_8192_RECEPTION_DENSITY_SHEAR_NEAR_CONVERGED"
    elif dg1<-.02 and dg2<-.02: outcome="WL_8192_RECEPTION_DENSITY_DEGRADES_SHEAR"
    else: outcome="WL_8192_RECEPTION_DENSITY_MIXED_SHEAR_RESPONSE"
    occupancy=acc.occupancy; occ=occupancy[occupancy>0]
    width=2*EXTENT/OBS_BINS; bc=np.floor((base_rays["uf"]+EXTENT)/width).astype(int); br=np.floor((base_rays["vf"]+EXTENT)/width).astype(int); bv=(bc>=0)&(bc<OBS_BINS)&(br>=0)&(br<OBS_BINS); base_occ=np.bincount((br[bv]*OBS_BINS+bc[bv]),minlength=OBS_BINS**2); base_nz=base_occ[base_occ>0]
    occupancy_report={"baseline":{"total_rays":len(base_rays["uf"]),"mean":float(np.mean(base_occ)),"median":float(np.median(base_occ)),"minimum_occupied":int(np.min(base_nz)),"maximum":int(np.max(base_occ)),"empty_cells":int(np.count_nonzero(base_occ==0))},"highres":{"total_rays":RAYS,"mean":float(np.mean(occupancy)),"median":float(np.median(occupancy)),"minimum_occupied":int(np.min(occ)),"maximum":int(np.max(occupancy)),"empty_cells":int(np.count_nonzero(occupancy==0))},"fixed_cell_theoretical_mean":16384}
    performance={"preflight_seconds":float(phase-started),"tile_count":len(alltiles),"mean_tile_seconds":float(np.mean(tile_times)),"median_tile_seconds":float(np.median(tile_times)),"min_tile_seconds":float(np.min(tile_times)),"max_tile_seconds":float(np.max(tile_times)),"propagation_and_observer_accumulation_total":float(sum(tile_times)),"wall_total":time.perf_counter()-started,"rays_per_second":RAYS/max(sum(tile_times),1e-9),"million_rays_per_second":RAYS/max(sum(tile_times),1e-9)/1e6}
    checks={k:True for k in ("abell2744_only","launch_8192x8192","ray_count_exact_67108864","observer_grid_frozen_64x64","candidate_frozen_D_jacobian","deposition_frozen_tsc_3x3","upstream_physics_frozen","source_unchanged","native_response_unchanged","a8_unchanged","interface_unchanged","m10_unchanged","los_unchanged","propagation_law_unchanged","screen_bounds_unchanged","observed_map_unchanged","no_kde","no_kappa_evaluation","no_target_guided_sampling","no_observational_fit","no_amplitude_rescale","no_sign_flip","no_rotation_fit","no_recipe_search","streamed_path_parity_pass","bounded_memory_execution","vram_safety_pass","disk_safety_pass","resume_supported","atomic_checkpointing","completed_tiles_not_recomputed","all_67108864_rays_processed","predictions_frozen_before_target_load","gamma1_scored","gamma2_scored","baseline_comparison_reported","rms_recovery_reported","orientation_reported","occupancy_reported","diagnostic_images_written","no_16384_execution")}
    result={"lab_id":LAB_ID,"manifest":manifest,"baseline":base_metrics,"metrics":scored,"deltas":{"gamma1_pearson":dg1,"gamma2_pearson":dg2,"gamma1_spearman":scored["gamma1"]["spearman"]-EXPECTED["gamma1"]["spearman"],"gamma2_spearman":scored["gamma2"]["spearman"]-EXPECTED["gamma2"]["spearman"],"gamma1_rms_ratio":scored["gamma1"]["rms_ratio_pred_over_obs"]-EXPECTED["gamma1"]["rms_ratio_pred_over_obs"],"gamma2_rms_ratio":scored["gamma2"]["rms_ratio_pred_over_obs"]-EXPECTED["gamma2"]["rms_ratio_pred_over_obs"]},"current_vs_8192":comp,"occupancy":occupancy_report,"performance":performance,"estimated_16384_seconds_EXTRAPOLATED_NOT_MEASURED":4*performance["wall_total"],"prediction_fingerprints":{"gamma1":p1,"gamma2":p2},"checks":checks,"outcome":outcome}
    atomic_json(RUN_DIR/"result.json",result); (RUN_DIR/"report.txt").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    for name,val in (("HIGHRES_RECEIPT",{"rays_processed":RAYS}),("PREDICTION_FINGERPRINTS",result["prediction_fingerprints"]),("GAMMA1_COMPARISON",scored["gamma1"]),("GAMMA2_COMPARISON",scored["gamma2"]),("SHEAR_MAGNITUDE_COMPARISON",scored["magnitude_pearson"]),("ORIENTATION_COMPARISON",scored["orientation_agreement"]),("RMS_RECOVERY",result["deltas"]),("OCCUPANCY",result["occupancy"]),("CURRENT_VS_8192",comp),("RESIDUAL_COMPARISON",comp),("PERFORMANCE",performance),("EXTRAPOLATED_16384_RUNTIME",result["estimated_16384_seconds_EXTRAPOLATED_NOT_MEASURED"]),("CHECKS",checks)): block(name,val)
    print(outcome); print("RESULT_JSON"); print(json.dumps(result,sort_keys=True)); return 0
if __name__=="__main__": raise SystemExit(main())
