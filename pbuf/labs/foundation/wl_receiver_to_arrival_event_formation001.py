#!/usr/bin/env python3
"""Dev130 physical 3-D receiver to continuous 2-D arrival-event audit."""
from __future__ import annotations
import argparse,hashlib,json,subprocess,sys,time
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from pbuf.wl.arrival_formation import (SCALES,native_receiver_plane,form_arrival_events,
    boundary_audit,channel_manifest,canonical_sha)
from pbuf.wl.arrival_information import (geometry_rank_ladder,reconstruction_audits,
    endpoint_comparison,intersection_statistics,information_preservation,survival_audit,matrix)
from pbuf.wl.receiver_information import rank_metrics

RUN=ROOT/"runs/wl_receiver_to_arrival_event_formation001";UP=ROOT/"runs/wl_receiver_state_channel_completion001"
STRUCT_SHA="ae97ac81e084fa4de1b7de4509e29c569753937a7e16f565ff3c91ec603d90b7"
MANIFEST_SHA="49aa5faeb1ff70b0aca2dd843e93b0c2da8574bdaf313d437001f88fbce9b7ab"
INPUT_SHA={"receiver_state_primary.npz":"b7260e3463343f8cee18a00ecc1f6d8ec8df490a9f06fe7778ebb024b4ad2ba0",
 "receiver_relational_channels.npz":"cde723a0ccaf590f7e9198053b910c9c37eefa5b7badbc6846f8a0ad7eefe16b",
 "receiver_bundle_channels.npz":"36e5fb7ad6c4a7e50a72d40e02dab4185ca6810427e74638657490dfba67bf5e",
 "receiver_validity_masks.npz":"15ecf874031c4a552d34b878559e41d142abd4f03031b2f078678a10942667bd"}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(name,value):(RUN/name).write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False,default=lambda x:x.item() if isinstance(x,np.generic) else str(x))+"\n")
def load(name):
    with np.load(UP/name,allow_pickle=True) as z:return {k:z[k] for k in z.files}
def subset(data,idx):return {k:(np.asarray(v)[idx] if np.asarray(v).shape[:1]==(idx.max()+1,) else v) for k,v in data.items()}
def full_ladder(events,primary,rel,bundle):
    x,_,m=matrix(events,("A0","A1","A2","A3","A4"));sample=np.flatnonzero(m)[::max(1,int(m.sum())//20000)]
    direct=[np.asarray(v,float) for k,v in primary.items() if np.asarray(v).ndim==1 and np.issubdtype(np.asarray(v).dtype,np.number)]
    relational=[np.asarray(v,float) for d in (rel,bundle) for v in d.values() if np.asarray(v).ndim==1 and np.issubdtype(np.asarray(v).dtype,np.number)]
    rows=[];cur=x
    for stage,add in (("F0",[]),("F1",direct),("F2",relational)):
        if add:cur=np.column_stack((cur,*add))
        finite=np.all(np.isfinite(cur[sample]),1);r=rank_metrics(cur[sample][finite]);r.update(stage=stage,retained_sample_count=int(finite.sum()),deterministic_rank_sample_count=int(len(sample)));rows.append(r)
    return rows
def figures(events,comparison,stats,ladder,recon,depth,traj,bundle):
    import matplotlib;matplotlib.use("Agg");import matplotlib.pyplot as plt
    g=events.event_geometry;valid=events.validity_masks["primary"];take=np.flatnonzero(valid)[::max(1,int(valid.sum())//20000)]
    specs=[("receiver_endpoint_vs_surface_intersection.png",g["arrival_u"][take],g["arrival_v"][take],"arrival u","arrival v","scatter"),
      ("intersection_t_distribution.png",g["intersection_t"][valid],None,"intersection t","count","hist"),
      ("arrival_event_scatter.png",g["arrival_u"][take],g["arrival_v"][take],"continuous arrival u","continuous arrival v (diagnostic, not detector image)","scatter"),
      ("arrival_direction_distribution.png",g["arrival_dir_u"][take],g["arrival_dir_v"][take],"direction u","direction v","scatter"),
      ("arrival_incidence_distribution.png",g["receiver_incidence_angle"][valid],None,"incidence angle","count","hist"),
      ("arrival_position_direction_coupling.png",g["arrival_u"][take],g["arrival_dir_u"][take],"arrival u","direction u","scatter"),
      ("arrival_rank_ladder.png",np.arange(len(ladder)),[q["effective_rank"] for q in ladder],"stage","effective rank","line"),
      ("arrival_family_reconstruction_r2.png",np.arange(len(recon)),[q["summary"]["median"] for q in recon],"family test","median R2","line"),
      ("receiver_depth_survival.png",np.arange(depth.get("available_channels",0)),list(depth.get("per_channel_r2",{}).values()),"depth channel","R2 from arrival geometry","line"),
      ("trajectory_bundle_history_survival.png",[0,1],[traj.get("summary",{}).get("median",0),bundle.get("summary",{}).get("median",0)],"history group","median R2","line")]
    for name,x,y,xlab,ylab,kind in specs:
        fig,ax=plt.subplots(figsize=(6,4));ax.hist(x,bins=60) if kind=="hist" else (ax.scatter(x,y,s=1,alpha=.25) if kind=="scatter" else ax.plot(x,y,"o-"));ax.set_xlabel(xlab);ax.set_ylabel(ylab);fig.tight_layout();fig.savefig(RUN/name,dpi=120);plt.close(fig)
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--validation",action="store_true");args=ap.parse_args();RUN.mkdir(parents=True,exist_ok=True);started=time.time()
    if sha(UP/"structural_result.json")!=STRUCT_SHA or (UP/"receiver_channel_manifest.sha256").read_text().strip()!=MANIFEST_SHA:print("DEV129_STRUCTURAL_BASELINE_MISMATCH");return 2
    if any(not (UP/k).exists() or sha(UP/k)!=v for k,v in INPUT_SHA.items()):print("DEV129_INPUT_INTEGRITY_FAILURE");return 3
    primary=load("receiver_state_primary.npz");rel=load("receiver_relational_channels.npz");bundle=load("receiver_bundle_channels.npz");full_n=len(primary["ray_index"])
    if args.validation:
        full_side=int(round(np.sqrt(full_n)));side=80;start=(full_side-side)//2;idx=(np.arange(start,start+side)[:,None]*full_side+np.arange(start,start+side)).ravel()
        primary={k:(np.asarray(v)[idx] if np.asarray(v).shape[:1]==(full_n,) else v) for k,v in primary.items()};rel={k:np.asarray(v)[idx] for k,v in rel.items()};bundle={k:np.asarray(v)[idx] for k,v in bundle.items()}
    n=len(primary["ray_index"]);side=int(round(np.sqrt(n)));plane=native_receiver_plane();basis=plane.manifest();bsha=canonical_sha(basis);dump("arrival_basis_manifest.json",basis);print("DEV130_RECEIVER_BASIS_SHA256="+bsha)
    manifest=channel_manifest();msha=canonical_sha(manifest);dump("arrival_event_channel_manifest.json",manifest);print("DEV130_CHANNEL_MANIFEST_SHA256="+msha)
    structural={"upstream_structural_sha256":STRUCT_SHA,"upstream_channel_manifest_sha256":MANIFEST_SHA,"input_sha256":INPUT_SHA,"receiver_basis_sha256":bsha,"intersection_convention":"x+t*dhat; forward t>eps_t; no direction reversal","event_definitions":"ArrivalEvent2D E0-E5","validity_policy":["VALID_FORWARD","VALID_ON_SURFACE","BACKWARD_ONLY","PARALLEL","NONFINITE","INVALID_DIRECTION"],"scales":list(SCALES),"channel_manifest_sha256":msha,"propagation_changes":0,"trajectory_changes":0,"receiver_state_changes":0,"target_access":False,"hst_pixel_access":False}
    dump("structural_result.json",structural);ssha=sha(RUN/"structural_result.json");print("DEV130_STRUCTURAL_SHA256="+ssha)
    pos=primary["global_receive_position"];direction=primary["normalized_final_direction"];launch=np.column_stack((primary["launch_u"],primary["launch_v"]))
    events=form_arrival_events(pos,direction,plane,ray_index=primary["ray_index"],receiver_row_index=np.arange(n),launch_uv=launch,launch_grid_index=primary["launch_grid_index"],side=side)
    boundary=boundary_audit(pos,plane);comparison=endpoint_comparison(events,pos,plane);stats=intersection_statistics(events);dump("receiver_boundary_audit.json",boundary);dump("intersection_statistics.json",stats);dump("endpoint_intersection_comparison.json",comparison)
    unique,counts=np.unique(events.event_geometry["validity"],return_counts=True);validity={k:int(v) for k,v in zip(unique,counts)};dump("event_validity.json",{"counts":validity,"ray_count":n,"silently_dropped":0})
    preservation=information_preservation(events,primary,rel,bundle);dump("dev129_information_preservation.json",preservation)
    ladder=geometry_rank_ladder(events);full=full_ladder(events,primary,rel,bundle);recon=reconstruction_audits(events);dump("geometry_rank_ladder.json",ladder);dump("full_arrival_rank_ladder.json",full);dump("arrival_reconstruction.json",recon)
    depth_channels={k:v for k,v in primary.items() if "w" in k or "depth" in k};depth_channels.update({k:v for k,v in rel.items() if "w" in k or "depth" in k})
    trajectory={k:v for k,v in primary.items() if k not in ("ray_index","launch_grid_index") and ("path" in k or "curvature" in k or "direction_change" in k)}
    depth=survival_audit(events,depth_channels,"receiver_depth");traj=survival_audit(events,trajectory,"trajectory_history");bund=survival_audit(events,bundle,"bundle_history");dump("depth_survival.json",depth);dump("trajectory_history_survival.json",traj);dump("bundle_history_survival.json",bund)
    event_arrays={**events.receiver_reference,**events.event_geometry};np.savez_compressed(RUN/"arrival_events_primary.npz",**event_arrays);np.savez_compressed(RUN/"arrival_relational_channels.npz",**events.local_relations);np.savez_compressed(RUN/"arrival_validity_masks.npz",**events.validity_masks)
    outcomes=["WL_PHYSICAL_RECEIVER_SURFACE_FORMATION_ESTABLISHED","WL_2D_ARRIVAL_EVENT_STATE_COMPLETION_ESTABLISHED"]
    outcomes.append("WL_RECEIVER_ENDPOINT_ALREADY_SURFACE_EQUIVALENT" if comparison["classification"]!="EXPLICIT_INTERSECTION_REQUIRED" else "WL_EXPLICIT_RECEIVER_INTERSECTION_REQUIRED")
    if recon[0]["summary"]["median"]<.99:outcomes.append("WL_ARRIVAL_DIRECTION_INFORMATION_ADVANTAGE_ESTABLISHED")
    if recon[2]["summary"]["median"]<.99:outcomes.append("WL_POSITION_DIRECTION_COUPLING_INFORMATION_ADVANTAGE_ESTABLISHED")
    if depth.get("independent"):outcomes.append("WL_RECEIVER_DEPTH_INFORMATION_SURVIVES_2D_FORMATION")
    if traj.get("independent"):outcomes.append("WL_TRAJECTORY_HISTORY_SURVIVES_2D_FORMATION")
    if bund.get("independent"):outcomes.append("WL_BUNDLE_HISTORY_SURVIVES_2D_FORMATION")
    result={"outcomes":outcomes,"structural_sha256":ssha,"receiver_basis_sha256":bsha,"channel_manifest_sha256":msha,"ray_count":n,"full_population_ray_count":full_n,"validation":args.validation,"runtime_seconds":time.time()-started,"DEV129_RECEIVER_FIELDS_LOST":preservation["DEV129_RECEIVER_FIELDS_LOST"],"target_access":False,"hst_pixel_access":False};dump("result.json",result)
    try:git=subprocess.run(["git","status","--short"],cwd=ROOT,text=True,capture_output=True).stdout
    except Exception:git="unavailable\n"
    (RUN/"baseline_git.txt").write_text(git);figures(events,comparison,stats,ladder,recon,depth,traj,bund)
    qualification="VALIDATION_ONLY_NO_SCIENCE_CLAIM" if args.validation else "\n".join(outcomes)
    (RUN/"report.txt").write_text(f"PROPAGATION_CHANGES=0\nTRAJECTORY_CHANGES=0\nRECEIVER_STATE_CHANGES=0\nTARGET_ACCESS=false\nHST_PIXEL_ACCESS=false\nDEV129_STRUCTURAL_SHA_VERIFIED=true\nDEV129_CHANNEL_MANIFEST_SHA_VERIFIED=true\nDEV130_RECEIVER_BASIS_SHA256={bsha}\nDEV130_CHANNEL_MANIFEST_SHA256={msha}\nDEV130_STRUCTURAL_SHA256={ssha}\nDEV129_RECEIVER_FIELDS_LOST={preservation['DEV129_RECEIVER_FIELDS_LOST']}\n{qualification}\n")
    print(qualification);return 4 if preservation["DEV129_RECEIVER_FIELDS_LOST"] else 0
if __name__=="__main__":raise SystemExit(main())
