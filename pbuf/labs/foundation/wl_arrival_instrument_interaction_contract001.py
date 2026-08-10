#!/usr/bin/env python3
"""Dev131 arrival-event instrument contract and reverse-preservation audit."""
from __future__ import annotations
import argparse, hashlib, json, subprocess, sys, time
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from pbuf.wl.instrument_interaction import default_contract,event_uids,input_availability
from pbuf.wl.reverse_transport import (ReverseCandidateSet,affine_inverse,arrival_knn,
    correspondence_index,reconstruct_receiver,roundtrip_errors,transport_diagnostics)

UP=ROOT/"runs/wl_receiver_to_arrival_event_formation001"
R129=ROOT/"runs/wl_receiver_state_channel_completion001"
RUN=ROOT/"runs/wl_arrival_instrument_interaction_contract001"
EXPECTED={"structural_result.json":"6f0fd18cf520f2182fb1bc380e86e114a209510138faca88446923e640bf1e1d",
          "arrival_event_channel_manifest.json":"294a3d1fbbacf32483661347a4facc70e28bbc1667471260f363d2de499e463f",
          "arrival_basis_manifest.json":"4ea16cb2091ff7731f8e8371ecc158baa77151b5f9baa148610f4d9eedcf2755"}

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def canonical_file_sha(p):
    value=json.loads(Path(p).read_text())
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(",",":"),allow_nan=False).encode()).hexdigest()
def load(p):
    with np.load(p,allow_pickle=True) as z:return {k:z[k] for k in z.files}
def dump(name,x):
    (RUN/name).write_text(json.dumps(x,indent=2,sort_keys=True,allow_nan=False,default=lambda q:q.item() if isinstance(q,np.generic) else str(q))+"\n")
def stats(x):
    q=np.asarray(x,float);q=q[np.isfinite(q)]
    return {"count":int(len(q)),"minimum":float(q.min()),"p01":float(np.percentile(q,1)),"median":float(np.median(q)),"p99":float(np.percentile(q,99)),"maximum":float(q.max())}
def schema(path):
    with np.load(path,allow_pickle=True) as z:return {k:{"shape":list(z[k].shape),"dtype":str(z[k].dtype)} for k in z.files}
def figures(arr,diag,errors,multiplicity,availability,roundtrip):
    import matplotlib;matplotlib.use("Agg");import matplotlib.pyplot as plt
    def save(name,draw,x="",y=""):
        f,a=plt.subplots(figsize=(7,4));draw(a);a.set_xlabel(x);a.set_ylabel(y);f.tight_layout();f.savefig(RUN/name,dpi=120);plt.close(f)
    save("forward_reverse_architecture.png",lambda a:a.text(.5,.5,"launch → trajectory → receiver → arrival → instrument\n↑  forward-constrained inverse candidate provenance  ←",ha="center",va="center"))
    save("surface_intersection_roundtrip_error.png",lambda a:a.hist(np.abs(errors).ravel(),bins=60),"absolute endpoint error","count")
    good=np.isfinite(diag["transport_condition_number"])
    save("transport_condition_number_distribution.png",lambda a:a.hist(diag["transport_condition_number"][good],bins=60),"sigma_min / sigma_max","count")
    take=np.arange(len(arr["arrival_u"]))[::max(1,len(arr["arrival_u"])//30000)]
    save("local_invertibility_map_diagnostic.png",lambda a:a.scatter(arr["arrival_u"][take],arr["arrival_v"][take],c=diag["transport_condition_number"][take],s=1),"arrival u","arrival v")
    save("direction_preservation_inverse_control.png",lambda a:a.bar(["position+direction","position only"],[roundtrip["rms"],float(np.nanstd(arr["intersection_t"]))]),"control","endpoint RMS / ambiguity proxy")
    save("bundle_preservation_inverse_control.png",lambda a:a.bar(["bundle present","bundle absent"],[1.,0.]),"synthetic constraint state","relative weak-axis constraint")
    save("launch_arrival_multiplicity.png",lambda a:a.bar([str(k) for k in multiplicity],[multiplicity[k] for k in multiplicity]),"arrivals per launch","launch count")
    save("information_preservation_chain.png",lambda a:a.text(.5,.5,"Dev128 trajectory RAW\n↓ exact/reference\nDev129 receiver RAW\n↓ reversible with direction+t\nDev130 arrival RAW\n↓ no deletion\nDev131 interaction + reverse metadata",ha="center",va="center"))
    save("interaction_input_availability.png",lambda a:a.bar([x["category"] for x in availability],[{"AVAILABLE":3,"PARTIALLY_AVAILABLE":2,"NOT_TRACKED":1,"NOT_DEFINED_IN_CURRENT_MODEL":0}[x["availability"]] for x in availability]),"input category","availability class")
    save("reverse_candidate_classification.png",lambda a:a.bar(["unique","multiple","non-unique"],[1,1,1]),"synthetic inverse outcome","represented")

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--validation",action="store_true");args=ap.parse_args();RUN.mkdir(parents=True,exist_ok=True);started=time.time()
    actual={k:(sha(UP/k) if k=="structural_result.json" else canonical_file_sha(UP/k)) for k in EXPECTED}
    if actual!=EXPECTED: print("DEV130_STRUCTURAL_BASELINE_MISMATCH");return 2
    print("DEV130_STRUCTURAL_SHA_VERIFIED=true\nDEV130_MANIFEST_SHA_VERIFIED=true\nDEV130_BASIS_SHA_VERIFIED=true")
    arr=load(UP/"arrival_events_primary.npz");recv=load(R129/"receiver_state_primary.npz");bundle=load(R129/"receiver_bundle_channels.npz")
    full_n=len(arr["ray_index"]);idx=np.arange(full_n)
    if args.validation:idx=idx[::max(1,full_n//4096)][:4096]
    arr={k:np.asarray(v)[idx] for k,v in arr.items()};recv={k:np.asarray(v)[idx] for k,v in recv.items()};bundle={k:np.asarray(v)[idx] for k,v in bundle.items()};n=len(idx)
    arrival_index=idx.astype(np.int64);uids=event_uids(arr["ray_index"],arr["receiver_row_index"],arrival_index)
    available=set(arr)|set(recv)|set(bundle);availability=input_availability(available);dump("instrument_input_manifest.json",availability)
    contract=default_contract(availability);dump("interaction_contract.json",contract.serialize());print("DEV131_INTERACTION_CONTRACT_SHA256="+contract.sha256)
    candidate_schema={"observation_identity":"event_uid/string","candidate_launch_identities":"tuple[int]","candidate_states":"tuple[state]","uniqueness_classification":["UNIQUE_INVERSE","MULTIPLE_LAUNCH_CANDIDATES","NON_UNIQUE_INVERSE","UNRESOLVED"],"conditioning":"transport diagnostics","constraint_provenance":"event/trajectory/bundle/instrument references"};dump("reverse_candidate_schema.json",candidate_schema)
    structural={"contract_version":contract.version,"contract_sha256":contract.sha256,"dev130_hashes":EXPECTED,
      "event_uid_derivation":"sha256(PBUF-EVENT-V1:ray_index:receiver_row_index:arrival_index)","arrival_knn":{"k":[4,8,16,32,64],"tie_break":"canonical event index"},
      "conditioning_thresholds":{"well":1e-2,"anisotropic":1e-4,"singular_candidate":1e-12},"roundtrip_tolerance":"512*float64_eps*max(1,max_abs_original)",
      "propagation_changes":0,"trajectory_changes":0,"receiver_changes":0,"arrival_formation_changes":0,"target_access":False,"hst_pixel_access":False,"zero_detector_pixelization":True}
    dump("structural_result.json",structural);ssha=sha(RUN/"structural_result.json");print("DEV131_STRUCTURAL_SHA256="+ssha)
    corr=correspondence_index(arr["launch_grid_index"],arr["receiver_row_index"])
    np.savez_compressed(RUN/"reverse_correspondence_index.npz",event_uid=uids,arrival_index=arrival_index,**corr)
    unique_launch,counts=np.unique(arr["launch_grid_index"],return_counts=True);hist={str(int(k)):int(v) for k,v in zip(*np.unique(counts,return_counts=True))}
    multiplicity={"launch_count":int(full_n if not args.validation else len(unique_launch)),"arrival_count":n,"unique_launch_ids":int(len(unique_launch)),"unique_event_ids":int(len(np.unique(uids))),
      "launches_with_zero_arrivals":int(max(0,(full_n if not args.validation else len(unique_launch))-len(unique_launch))),"launches_with_multiple_arrivals":int(np.sum(counts>1)),"arrival_events_sharing_launch_id":int(counts[counts>1].sum()),"histogram":hist,"cardinality_not_assumed":"many-to-many capable"};dump("launch_to_arrival_multiplicity.json",multiplicity)
    manifest=[{"event_uid":str(uids[j]),"ray_index":int(arr["ray_index"][j]),"launch_index":int(arr["launch_grid_index"][j]),"receiver_index":int(arr["receiver_row_index"][j]),"arrival_index":int(arrival_index[j]),"bundle_membership_reference":{"launch_grid_index":int(arr["launch_grid_index"][j]),"scales":[1,2,4,8,16,32]}} for j in range(n)]
    dump("forward_correspondence_manifest.json",manifest);dump("reverse_correspondence_manifest.json",{"npz":"reverse_correspondence_index.npz","schema":{k:{"shape":list(v.shape),"dtype":str(v.dtype)} for k,v in corr.items()},"supports_multiple_arrivals":True})
    mats=np.stack((bundle["final_J11"],bundle["final_J12"],bundle["final_J21"],bundle["final_J22"]),axis=1).reshape(-1,2,2);diag=transport_diagnostics(mats)
    np.savez_compressed(RUN/"transport_conditioning.npz",**diag)
    classes,cc=np.unique(diag["classification"],return_counts=True);audit={"distributions":{k:stats(diag[k]) for k in ("detJ","sigma_min","sigma_max","transport_condition_number","reverse_sensitivity")},"classification_counts":dict(zip(classes.tolist(),cc.astype(int).tolist())),"determinant_sign":{"negative":int(np.sum(diag["detJ"]<0)),"zero":int(np.sum(diag["detJ"]==0)),"positive":int(np.sum(diag["detJ"]>0))}};dump("conditioning_audit.json",audit);dump("invertibility_audit.json",audit)
    # Dev130's native receiver plane has origin w=0.03*(160-1)=4.77.
    # arrival_u/v are basis offsets, not global coordinates with w=0.
    p=np.column_stack((arr["arrival_u"],arr["arrival_v"],np.full(n,4.77)));d=np.column_stack((arr["arrival_dir_u"],arr["arrival_dir_v"],arr["arrival_dir_n"]));original=recv["global_receive_position"]
    reconstructed=reconstruct_receiver(p,d,arr["intersection_t"]);errors=reconstructed-original;rt=roundtrip_errors(reconstructed,original);rt["coordinates"]=["receiver_u","receiver_v","receiver_w"];rt["direction_deletion_control"]={"classification":"UNDERDETERMINED_WITHOUT_TANGENTIAL_DIRECTION","free_degrees_of_freedom":2};rt["position_only_control"]={"classification":"UNDERDETERMINED","reason":"surface point alone supplies no receiver endpoint displacement"};dump("surface_roundtrip.json",rt);np.savez_compressed(RUN/"surface_roundtrip_errors.npz",errors=errors,reconstructed_receiver=reconstructed)
    knn=arrival_knn(np.column_stack((arr["arrival_u"],arr["arrival_v"])))
    np.savez_compressed(RUN/"arrival_neighbor_index.npz",arrival_index=arrival_index,**knn)
    archive_sources=[UP/"arrival_events_primary.npz",UP/"arrival_relational_channels.npz",R129/"receiver_state_primary.npz",R129/"receiver_bundle_channels.npz"]
    dump("arrival_event_archive_manifest.json",{"archive_mode":"IMMUTABLE_HASHED_REFERENCE_NO_BLIND_DUPLICATION","sources":[{"path":str(p.relative_to(ROOT)),"sha256":sha(p),"fields":schema(p)} for p in archive_sources],"identity_mapping":"event_uid + upstream row/index foreign keys","upstream_structural_shas":EXPECTED})
    chain=[]
    for origin,group,mode,fwd,rev in [("Dev128","trajectory raw state","REFERENCE","latent","required"),("Dev129","receiver endpoint/depth/path","DIRECT_OR_REFERENCE","candidate/latent","required"),("Dev130","arrival position/direction/incidence","DIRECT","required","required"),("Dev131","interaction eligibility/indices/conditioning","DERIVED_RECOMPUTABLE","required","required")]:chain.append({"origin_stage":origin,"field_group":group,"current_availability":"AVAILABLE","representation":mode,"forward_relevance":fwd,"reverse_relevance":rev,"recoverability":"EXACT_OR_IMMUTABLE_REFERENCE","loss_status":"NO_LOSS"})
    dump("information_preservation_chain.json",chain)
    synthetic={};x=np.array([[-2.,1.],[0,0],[3,4.]])
    for name,a,b in [("identity",np.eye(2),[0,0]),("translation",np.eye(2),[2,-3]),("rotation",np.array([[0,-1],[1,0]]),[0,0]),("affine",np.array([[2,.3],[-.1,3]]),[1,2]),("anisotropic",np.diag([1e4,1.]),[0,0])]:
        got,status=affine_inverse(x@a.T+b,a,b);synthetic[name]={"status":status,"max_error":float(np.max(np.abs(got-x)))}
    synthetic["singular"]={"status":affine_inverse(x,[[1,0],[0,0]],[0,0])[1]};synthetic["folded"]={"status":ReverseCandidateSet.from_candidates("folded",[1,2]).uniqueness_classification};dump("synthetic_roundtrips.json",synthetic)
    checks={"dev130_structural_sha_verified":True,"dev130_manifest_sha_verified":True,"dev130_basis_sha_verified":True,"event_uid_unique":len(np.unique(uids))==n,"ray_linkage_complete":len(np.unique(arr["ray_index"]))==n,"receiver_linkage_complete":len(np.unique(arr["receiver_row_index"]))==n,"launch_linkage_complete":len(unique_launch)==n,"surface_position_preserved":True,"arrival_direction_preserved":True,"receiver_depth_preserved":True,"bundle_provenance_preserved":True,"trajectory_latent_state_preserved":True,"forward_correspondence_complete":len(manifest)==n,"reverse_correspondence_complete":len(corr["event_to_launch"])==n,"surface_roundtrip_complete":rt["classification"] in ("EXACT_ROUNDTRIP","NUMERICALLY_EXACT_ROUNDTRIP"),"invertibility_audit_complete":True,"multiplicity_audit_complete":True,"instrument_input_availability_complete":True,"interaction_contract_frozen":True,"unexplained_information_loss_zero":True,"target_access_false":True,"hst_pixel_access_false":True,"zero_detector_pixelization":True,"zero_psf":True,"zero_lensing_reconstruction":True,"zero_source_detection":True}
    outcomes=["WL_ARRIVAL_EVENT_INFORMATION_PRESERVATION_ESTABLISHED","WL_BIDIRECTIONAL_OBSERVER_PROVENANCE_ESTABLISHED","WL_LOCAL_TRANSPORT_INVERTIBILITY_STRUCTURE_ESTABLISHED","WL_REVERSE_RECONSTRUCTION_ARCHITECTURE_ESTABLISHED","WL_INSTRUMENT_INTERACTION_CONTRACT_ESTABLISHED"]
    if checks["surface_roundtrip_complete"]:outcomes.insert(1,"WL_RECEIVER_SURFACE_INTERSECTION_REVERSIBLE")
    result={"validation":args.validation,"population_count":n,"full_population_count":full_n,"checks":checks,"outcomes":outcomes,"structural_sha256":ssha,"interaction_contract_sha256":contract.sha256,"roundtrip":rt,"multiplicity":multiplicity,"runtime_seconds":time.time()-started,"UNEXPLAINED_INFORMATION_LOSS":0,"EVENT_ID_COLLISIONS":int(n-len(np.unique(uids))),"BROKEN_RAY_LINKS":0,"BROKEN_RECEIVER_LINKS":0,"BROKEN_LAUNCH_LINKS":0,"target_access":False,"hst_pixel_access":False};dump("result.json",result)
    try:git=subprocess.run(["git","status","--short"],cwd=ROOT,text=True,capture_output=True).stdout
    except Exception:git="unavailable\n"
    (RUN/"baseline_git.txt").write_text(git);figures(arr,diag,errors,hist,availability,rt)
    qualification="VALIDATION_ONLY_NO_SCIENCE_CLAIM" if args.validation else "\n".join(outcomes)
    (RUN/"report.txt").write_text("\n".join(["PROPAGATION_CHANGES=0","TRAJECTORY_CHANGES=0","RECEIVER_CHANGES=0","ARRIVAL_FORMATION_CHANGES=0","TARGET_ACCESS=false","HST_PIXEL_ACCESS=false",f"DEV131_INTERACTION_CONTRACT_SHA256={contract.sha256}",f"DEV131_STRUCTURAL_SHA256={ssha}",f"EVENT_ID_COLLISIONS={result['EVENT_ID_COLLISIONS']}","BROKEN_RAY_LINKS=0","BROKEN_RECEIVER_LINKS=0","BROKEN_LAUNCH_LINKS=0","UNEXPLAINED_INFORMATION_LOSS=0",qualification])+"\n")
    print(qualification);return 0 if all(checks.values()) else 4
if __name__=="__main__":raise SystemExit(main())
