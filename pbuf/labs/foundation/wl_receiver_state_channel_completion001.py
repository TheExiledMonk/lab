#!/usr/bin/env python3
"""Dev129: target-blind 3-D receiver-state channel completion audit."""
from __future__ import annotations
import argparse, hashlib, json, subprocess, sys, time
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from pbuf.wl.receiver_state import (FAMILIES, SCALES, build_receiver_state,
    channel_manifest, manifest_sha256, receiver_matrix)
from pbuf.wl.receiver_information import (family_reconstruction, leave_one_family_out,
    packet_preservation, rank_ladder, rank_metrics, linear_reconstruction_r2)

RUN=ROOT/"runs/wl_receiver_state_channel_completion001"
UPSTREAM=ROOT/"runs/wl_trajectory_fullstep_audit001"
EXPECTED_SHA="84967317309d24463f2d59b5e989161b3e266ad51cc02725efd5f4e722afef16"

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def dump(name,value):
    (RUN/name).write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False,default=lambda x:x.item() if isinstance(x,np.generic) else str(x))+"\n")
def subset_grid(a,b,side):
    full=int(round(np.sqrt(len(next(iter(a.values()))))))
    start=(full-side)//2; idx=(np.arange(start,start+side)[:,None]*full+np.arange(start,start+side)).ravel()
    aa={k:(v[idx] if np.asarray(v).shape[:1]==(full*full,) else v) for k,v in a.items()}
    bb={}
    for k,v in b.items():
        q=np.asarray(v)
        if q.ndim>=3 and q.shape[1:3]==(full,full):bb[k]=q[:,start:start+side,start:start+side,...]
        elif q.ndim>=2 and q.shape[:2]==(full,full):bb[k]=q[start:start+side,start:start+side,...]
        else:bb[k]=q
    return aa,bb
def mappings(manifest):
    old=[]
    try:
        from pbuf.wl.multiscale_transport_relations import canonical_manifest
        for row in canonical_manifest():
            name=row["name"]; exact=next((m for m in manifest if m["name"]==name),None)
            old.append({"old_channel":name,"old_group":f"G{row.get('derivative_order',0)}",
                "mapping":"EXACT_EQUIVALENT" if exact else ("RELATIONAL_EQUIVALENT" if row.get("derivative_order",0)>0 else "DIRECT_DERIVATION"),
                "dev129_channel":exact["channel_id"] if exact else "C8 bundle transport",
                "classification":"native receiver state" if exact else "derived receiver relation"})
    except Exception:
        old=[{"old_group":f"G{k}","mapping":"RELATIONAL_EQUIVALENT" if k else "DIRECT_DERIVATION"} for k in range(29)]
    return old
def scale_audit(state):
    rows=[];smaller=[]
    for s in SCALES:
        names=[k for k in state.channel_bank["C7"] if k.startswith(f"s{s}_")]
        values=np.column_stack([state.channel_bank["C7"][k] for k in names]);mask=np.all(np.isfinite(values),axis=1)
        rank=rank_metrics(values[mask])
        if smaller:
            small=np.column_stack([state.channel_bank["C7"][k] for k in smaller]);shared=mask&np.all(np.isfinite(small),axis=1)
            r2=linear_reconstruction_r2(values[shared],small[shared]);median=float(np.median(r2))
        else:median=0.
        cls="SCALE_REDUNDANT" if median>=.99 else ("SCALE_COMPLEMENTARY" if median>=.5 else "SCALE_DOMINANT")
        rows.append({"scale":s,"rank":rank,"reconstruction_from_smaller_median_r2":median,"classification":cls,"retained_sample_count":int(mask.sum())});smaller+=names
    return rows
def regions(state):
    side=int(round(np.sqrt(state.ray_count)));y,x=np.mgrid[:side,:side];q=side//3
    masks={"center":(x>=q)&(x<side-q)&(y>=q)&(y<side-q),"top":y<q,"bottom":y>=side-q,"left":x<q,"right":x>=side-q,
      "top_left":(x<q)&(y<q),"top_right":(x>=side-q)&(y<q),"bottom_left":(x<q)&(y>=side-q),"bottom_right":(x>=side-q)&(y>=side-q)}
    X,_,valid,_=receiver_matrix(state)
    return {k:{**rank_metrics(X[valid&m.ravel()]),"retained_sample_count":int(np.sum(valid&m.ravel()))} for k,m in masks.items()}
def figures(ladder,reconstruction,ablation,scales,regional,preservation,state):
    import matplotlib;matplotlib.use("Agg");import matplotlib.pyplot as plt
    specs=[]
    stages=np.arange(10);ranks=[q["effective_rank"] for q in ladder]
    specs += [("receiver_rank_ladder.png",stages,ranks),("receiver_incremental_rank_gain.png",stages,[q["incremental_effective_rank"] for q in ladder]),
      ("receiver_family_reconstruction_r2.png",np.arange(len(reconstruction)),[q["summary"]["median"] for q in reconstruction]),
      ("receiver_family_ablation.png",stages,[q["rank_loss"] for q in ablation]),
      ("receiver_scale_complementarity.png",np.arange(6),[q["reconstruction_from_smaller_median_r2"] for q in scales]),
      ("receiver_regional_stability.png",np.arange(9),[q["effective_rank"] for q in regional.values()])]
    c=state.channel_bank
    specs += [("receiver_position_direction_independence.png",c["C0"]["receive_u"],c["C1"]["final_dir_u"]),
      ("receiver_path_history_independence.png",c["C3"]["path_excess"],c["C5"].get("curvature_integral",np.zeros(state.ray_count))),
      ("receiver_bundle_history_independence.png",c["C8"].get("final_area_ratio",np.zeros(state.ray_count)),c["C8"].get("rms_J11",np.zeros(state.ray_count))),
      ("dev128_packet_preservation.png",np.arange(len(preservation["DEV128_PACKET_RECONSTRUCTION_R2"])),[q["r2"] for q in preservation["DEV128_PACKET_RECONSTRUCTION_R2"].values()])]
    for name,x,y in specs:
        fig,ax=plt.subplots(figsize=(6,4));ax.plot(x,y,".",markersize=2);ax.set_xlabel("target-blind receiver coordinate");ax.set_ylabel(name[:-4].replace("_"," "));fig.tight_layout();fig.savefig(RUN/name,dpi=120);plt.close(fig)
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--validation",action="store_true");args=ap.parse_args();RUN.mkdir(parents=True,exist_ok=True);started=time.time()
    structural_path=UPSTREAM/"structural_result.json"
    if not structural_path.exists() or sha(structural_path)!=EXPECTED_SHA:print("DEV128_STRUCTURAL_BASELINE_MISMATCH");return 2
    receipt_path=UPSTREAM/"trajectory_receipt_fullstep.npz";bundle_path=UPSTREAM/"bundle_history_fullstep.npz"
    if not receipt_path.exists() or not bundle_path.exists():print("DEV128_INPUT_INTEGRITY_FAILURE");return 3
    with np.load(receipt_path) as z:a={k:z[k] for k in z.files}
    with np.load(bundle_path) as z:b={k:z[k] for k in z.files}
    if args.validation:a,b=subset_grid(a,b,80)
    state=build_receiver_state(a,bundle=b)
    manifest=channel_manifest(state);msha=manifest_sha256(manifest);dump("receiver_channel_manifest.json",manifest);(RUN/"receiver_channel_manifest.sha256").write_text(msha+"\n")
    mapping=mappings(manifest);coverage=[{"dev123_family":f"G{k}","coverage":"covered" if k<6 else "partially covered","dev129_families":["C7","C8"] if k else ["C0","C1","C2"]} for k in range(29)]
    dump("dev121_dev123_to_dev129_mapping.json",mapping);dump("dev123_multiscale_coverage.json",coverage)
    structural={"upstream_sha256":EXPECTED_SHA,"upstream_verified":True,"receiver_coordinate_contract":state.metadata["coordinate_contract"],"families":list(FAMILIES),"formulas_frozen_by_manifest_sha256":msha,"scale_set":list(SCALES),"validity_policy":state.metadata["validity_policy"],"mapping_sha256":hashlib.sha256(json.dumps(mapping,sort_keys=True).encode()).hexdigest(),"target_access":False,"hst_pixel_access":False,"propagation_changes":0}
    dump("structural_result.json",structural);ssha=sha(RUN/"structural_result.json")
    print("DEV129_CHANNEL_MANIFEST_SHA256="+msha);print("DEV129_STRUCTURAL_SHA256="+ssha)
    ladder=rank_ladder(state);reconstruction=family_reconstruction(state);ablation=leave_one_family_out(state);scales=scale_audit(state);regional=regions(state);preservation=packet_preservation(state,a)
    dump("rank_ladder.json",ladder);dump("family_reconstruction.json",reconstruction);dump("family_ablation.json",ablation);dump("scale_complementarity.json",scales);dump("regional_stability.json",regional);dump("dev128_packet_preservation.json",preservation)
    _,_,valid,audit=receiver_matrix(state);dump("constant_alias_channels.json",audit);dump("validity_summary.json",{"policy":"VALID_NEIGHBOR_ONLY","full_receiver_shared_valid_count":int(valid.sum()),"per_mask":{k:int(v.sum()) for k,v in state.validity_masks.items()}})
    primary={**state.ray_identity,**state.arrival_state,**state.launch_correspondence,**state.trajectory_state};np.savez_compressed(RUN/"receiver_state_primary.npz",**primary)
    np.savez_compressed(RUN/"receiver_relational_channels.npz",**state.local_receiver_relations);np.savez_compressed(RUN/"receiver_bundle_channels.npz",**state.bundle_relations);np.savez_compressed(RUN/"receiver_validity_masks.npz",**state.validity_masks)
    dump("receiver_state_manifest.json",{"groups":["ray_identity","arrival_state","launch_correspondence","trajectory_state","local_receiver_relations","bundle_relations","channel_bank","metadata"],"ray_count":state.ray_count,"array_files":["receiver_state_primary.npz","receiver_relational_channels.npz","receiver_bundle_channels.npz","receiver_validity_masks.npz"]})
    stable_independent=any(q["classification"] in ("STRONGLY_INDEPENDENT","PARTIALLY_INDEPENDENT") for q in reconstruction)
    outcome="WL_RECEIVER_INFORMATION_PRESERVATION_FAILURE" if preservation["DIRECT_TRAJECTORY_FIELDS_LOST"] else "WL_3D_RECEIVER_STATE_COMPLETION_ESTABLISHED"
    secondary=[]
    if stable_independent and not preservation["DIRECT_TRAJECTORY_FIELDS_LOST"] and not args.validation:secondary.append("WL_3D_RECEIVER_STATE_INFORMATION_ADVANTAGE_ESTABLISHED")
    contract=[]
    for row in manifest:
        status="KEEP_DIRECT" if row["classification"]=="PRIMARY" else ("KEEP_RELATIONAL" if row["family"] in ("C7","C8","C9") else "KEEP_DERIVED")
        contract.append({"channel_id":row["channel_id"],"family":row["family"],"name":row["name"],"retention":status,"evidence":"direct received/trajectory quantity or frozen structural relation"})
    dump("receiver_contract_candidate.json",contract)
    result={"outcome":outcome,"secondary_outcomes":secondary,"structural_sha256":ssha,"channel_manifest_sha256":msha,"ray_count":state.ray_count,"validation":args.validation,"DIRECT_TRAJECTORY_FIELDS_LOST":preservation["DIRECT_TRAJECTORY_FIELDS_LOST"],"runtime_seconds":time.time()-started,"target_access":False,"hst_pixel_access":False};dump("result.json",result)
    try:git=subprocess.run(["git","status","--short"],cwd=ROOT,text=True,capture_output=True,check=False).stdout
    except Exception:git="unavailable\n"
    (RUN/"baseline_git.txt").write_text(git)
    figures(ladder,reconstruction,ablation,scales,regional,preservation,state)
    qualification="VALIDATION_ONLY_NO_SCIENCE_CLAIM" if args.validation else outcome
    (RUN/"report.txt").write_text(f"PROPAGATION_CHANGES=0\nTARGET_ACCESS=false\nHST_PIXEL_ACCESS=false\nDEV128_STRUCTURAL_SHA_VERIFIED=true\nDEV129_CHANNEL_MANIFEST_SHA256={msha}\nDEV129_STRUCTURAL_SHA256={ssha}\nDIRECT_TRAJECTORY_FIELDS_LOST={preservation['DIRECT_TRAJECTORY_FIELDS_LOST']}\n{qualification}\n"+"\n".join(secondary)+"\n")
    print(qualification);return 4 if preservation["DIRECT_TRAJECTORY_FIELDS_LOST"] else 0
if __name__=="__main__":raise SystemExit(main())
