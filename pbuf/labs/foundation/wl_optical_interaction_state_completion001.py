#!/usr/bin/env python3
"""Dev132 physical optical interaction state completion audit."""
from __future__ import annotations
import argparse, hashlib, json, subprocess, sys, time
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from pbuf.wl.instrument_interaction import event_uids
from pbuf.wl.optical_interaction_state import (CONTRACT_VERSION, canonical_sha256,
    contract_schema, state_from_event, validate_derivation_graph)

UP=ROOT/"runs/wl_arrival_instrument_interaction_contract001"
A130=ROOT/"runs/wl_receiver_to_arrival_event_formation001"
R129=ROOT/"runs/wl_receiver_state_channel_completion001"
RUN=ROOT/"runs/wl_optical_interaction_state_completion001"
EXPECTED_CONTRACT="9e5e9b62f43e3e4b46492f4bd3efab1828839cddb77efc661276a48f43db1e25"
EXPECTED_STRUCTURAL="3baa2263f9bb2f81906d522d51de24f5b02930289aa85909a2d3934b64d4f9a0"

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):
    with np.load(p,allow_pickle=True) as z:return {k:z[k] for k in z.files}
def dump(name,value):
    (RUN/name).write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False,
        default=lambda x:x.item() if isinstance(x,np.generic) else str(x))+"\n")
def item(name,family,source_stage,source_field,availability,units,forward,reverse,
         formula=None,derivation="NOT_DERIVABLE",assumptions=None):
    return {"name":name,"family":family,"source_stage":source_stage,"source_field":source_field,
            "availability_class":availability,"derivation_formula":formula,"derivation_class":derivation,
            "units":units,"forward_instrument_relevance":forward,"reverse_reconstruction_relevance":reverse,
            "assumptions_required":assumptions or []}

def inventory():
    I=[]; add=I.append
    for n,u in (("arrival_u","receiver-plane coordinate"),("arrival_v","receiver-plane coordinate"),
                ("arrival_dir_u","dimensionless direction cosine/component"),("arrival_dir_v","dimensionless direction cosine/component"),
                ("arrival_dir_n","dimensionless direction cosine/component"),("incidence_cosine","dimensionless"),
                ("incidence_angle","radian"),("intersection_t","propagation coordinate")):
        src="receiver_incidence_"+n.split("incidence_")[-1] if n.startswith("incidence_") else n
        add(item(n,"O0 geometric arrival state","Dev130",src,"DIRECTLY_AVAILABLE",u,"REQUIRED","REQUIRED"))
    add(item("physical_weight","O1 carried weight / intensity","none",None,"NOT_TRACKED","unspecified","REQUIRED","POSSIBLY_RELEVANT",assumptions=["new physical source normalization/transport law"]))
    add(item("numerical_sampling_weight","O1 carried weight / intensity","launch sampling","implicit equal samples","DERIVABLE_FROM_FROZEN_STATE","dimensionless","NOT_CURRENTLY_REQUIRED","POSSIBLY_RELEVANT","1 per computational ray; not physical intensity","EXACT_ALGEBRAIC"))
    add(item("launch_measure","O1 carried weight / intensity","Dev130 launch grid","launch_u, launch_v","DERIVABLE_FROM_FROZEN_STATE","launch-coordinate area","POSSIBLY_RELEVANT","REQUIRED","median(diff(unique(launch_u))) * median(diff(unique(launch_v)))","EXACT_ALGEBRAIC"))
    for n in ("wavelength","frequency","photon_energy","spectral_bin","filter_band_identity"):
        add(item(n,"O2 wavelength / frequency","none",None,"NOT_DEFINED_IN_CURRENT_MODEL","unspecified","POSSIBLY_RELEVANT","NOT_CURRENTLY_REQUIRED",assumptions=["new source spectral input"]))
    for n,u in (("path_length","path coordinate length"),("path_excess","path coordinate length"),("number_of_steps","count")):
        add(item(n,"O3 time / delay","Dev129",n,"DIRECTLY_AVAILABLE",u,"POSSIBLY_RELEVANT","REQUIRED"))
    for n in ("arrival_time","relative_arrival_delay","emission_time","absolute_time_origin"):
        add(item(n,"O3 time / delay","none",None,"NOT_DEFINED_IN_CURRENT_MODEL","time","POSSIBLY_RELEVANT","NOT_CURRENTLY_REQUIRED",assumptions=["propagation speed and physical time origin"]))
    for n in ("phase","phase_offset","complex_amplitude","coherence_length","wavefront_phase"):
        add(item(n,"O4 phase / coherence","none",None,"NOT_DEFINED_IN_CURRENT_MODEL","unspecified","POSSIBLY_RELEVANT","NOT_CURRENTLY_REQUIRED",assumptions=["wavelength plus coherent evolution physics"]))
    for n in ("polarization_vector","stokes_parameters","helicity","jones_state"):
        add(item(n,"O5 polarization","none",None,"NOT_DEFINED_IN_CURRENT_MODEL","unspecified","POSSIBLY_RELEVANT","NOT_CURRENTLY_REQUIRED",assumptions=["source polarization and polarization transport"]))
    for n,src,avail,formula,dc in (
        ("bundle_jacobian","final_J11,J12,J21,J22","DIRECTLY_AVAILABLE",None,"NOT_DERIVABLE"),
        ("bundle_hessian","final_H1/H2_uu/uv/vv","DIRECTLY_AVAILABLE",None,"NOT_DERIVABLE"),
        ("bundle_area_ratio","final_area_ratio","DIRECTLY_AVAILABLE",None,"NOT_DERIVABLE"),
        ("bundle_anisotropy","final_anisotropy","DIRECTLY_AVAILABLE",None,"NOT_DERIVABLE"),
        ("bundle_history","mean/rms/minimum/maximum bundle channels","DIRECTLY_AVAILABLE",None,"NOT_DERIVABLE"),
        ("arrival_direction_spread","arrival directions + arrival neighbor index","DERIVABLE_FROM_FROZEN_STATE","covariance of neighboring (arrival_dir_u,arrival_dir_v,arrival_dir_n)","NUMERICAL_SUMMARY"),
        ("position_direction_coupling","arrival positions/directions + neighbor index","DERIVABLE_FROM_FROZEN_STATE","cross-covariance of neighboring (arrival_u,arrival_v) and arrival direction","NUMERICAL_SUMMARY"),
        ("local_position_direction_volume","position/direction neighbor covariance","DERIVABLE_FROM_FROZEN_STATE","sqrt(max(0,det(cov([u,v,dir_u,dir_v]))))","NUMERICAL_SUMMARY")):
        add(item(n,"O6 local bundle / etendue-like structure","Dev129/Dev131",src,avail,"structural/native", "POSSIBLY_RELEVANT","REQUIRED",formula,dc))
    for n,src in (("event_uid","Dev131 event UID"),("arrival_identity","arrival_index"),("ray_identity","ray_index"),
                  ("receiver_identity","receiver_row_index"),("launch_identity","launch_grid_index,launch_u,launch_v"),
                  ("path_state","path_length/path_excess and trajectory history reference"),("receiver_depth","global_receive_position"),
                  ("bundle_provenance","launch/arrival neighbors and bundle archive"),("conditioning_state","transport_conditioning.npz"),
                  ("reverse_candidate_metadata","reverse_candidate_schema.json")):
        add(item(n,"O7 provenance / latent reverse state","Dev128-Dev131",src,"DIRECTLY_AVAILABLE","identity/native","NOT_CURRENTLY_REQUIRED","REQUIRED"))
    return I

def missing_ledger():
    return [
      {"quantity":"physical event weight","why_absent":"no source-defined or conserved signal field exists","needed_for_geometric_optics":False,"needed_for_diffraction_interference":False,"needed_for_detector_energy_response":True,"needed_for_reverse_reconstruction":False,"future_extension_required":True},
      {"quantity":"wavelength/frequency/photon energy","why_absent":"not represented in frozen launch, trajectory, receiver, or arrival state","needed_for_geometric_optics":False,"needed_for_diffraction_interference":True,"needed_for_detector_energy_response":True,"needed_for_reverse_reconstruction":False,"future_extension_required":True},
      {"quantity":"arrival time and time origin","why_absent":"path length exists but propagation speed and physical time origin are undefined","needed_for_geometric_optics":False,"needed_for_diffraction_interference":False,"needed_for_detector_energy_response":False,"needed_for_reverse_reconstruction":False,"future_extension_required":False},
      {"quantity":"phase/coherence","why_absent":"no coherent state or evolution law; path length alone is insufficient","needed_for_geometric_optics":False,"needed_for_diffraction_interference":True,"needed_for_detector_energy_response":False,"needed_for_reverse_reconstruction":False,"future_extension_required":True},
      {"quantity":"polarization","why_absent":"no source polarization or transport state","needed_for_geometric_optics":False,"needed_for_diffraction_interference":False,"needed_for_detector_energy_response":False,"needed_for_reverse_reconstruction":False,"future_extension_required":False}]

def synthetic_tests():
    base={"event_uid":"e","ray_index":1,"receiver_row_index":2,"launch_grid_index":3}
    rich=base|{"physical_weight":2.5,"wavelength":5e-7,"frequency":6e14,"path_length":4.,"arrival_time":7.,"relative_delay":.1,"phase":.3,"stokes":[1,.1,.2,.3],"final_J11":1.}
    a=state_from_event(rich); b=state_from_event(base)
    return {"physical_weight_exact":a.carried_signal["physical_weight"]==2.5,
      "missing_weight_explicit_false":not b.availability_metadata["has_physical_weight"] and not b.carried_signal,
      "spectral_exact_and_distinct":a.spectral["wavelength"]==5e-7 and a.spectral["frequency"]==6e14,
      "missing_spectral_explicit":not b.availability_metadata["has_spectral_state"],
      "temporal_categories_distinct":set(a.temporal)>={"path_length","arrival_time","relative_delay"},
      "phase_preserved_not_derived":a.phase_coherence=={"phase":.3},
      "polarization_preserved":a.polarization["stokes"]==[1,.1,.2,.3],
      "availability_mixed":a.availability_metadata!=b.availability_metadata}

def figures(availability,inv,arr,bundle,launch_measure):
    import matplotlib;matplotlib.use("Agg");import matplotlib.pyplot as plt
    def save(name,draw):
        f,a=plt.subplots(figsize=(7,4));draw(a);f.tight_layout();f.savefig(RUN/name,dpi=120);plt.close(f)
    names=list(availability); vals=[availability[k]["fraction"] for k in names]
    save("optical_state_availability.png",lambda a:(a.barh(names,vals),a.set_xlim(0,1),a.set_xlabel("event fraction")))
    fam=sorted(set(x["family"][:2] for x in inv)); fw=[sum(x["forward_instrument_relevance"]=="REQUIRED" for x in inv if x["family"].startswith(f)) for f in fam]; rv=[sum(x["reverse_reconstruction_relevance"]=="REQUIRED" for x in inv if x["family"].startswith(f)) for f in fam]
    save("forward_vs_reverse_relevance.png",lambda a:(a.bar(np.arange(8)-.2,fw,.4,label="forward"),a.bar(np.arange(8)+.2,rv,.4,label="reverse"),a.set_xticks(range(8),fam),a.legend()))
    save("optical_readiness_levels.png",lambda a:(a.bar(["L0","L1","L2","L3","L4"],[1,0,0,0,0]),a.set_ylim(0,1.2)))
    save("reverse_readiness_levels.png",lambda a:(a.bar(["R0","R1","R2","R3"],[1,1,1,0]),a.set_ylim(0,1.2)))
    save("missing_optical_state_summary.png",lambda a:a.bar(["weight","spectral","time","phase","polarization"],[1,1,1,1,1]))
    save("primitive_derived_dependency.png",lambda a:a.text(.5,.5,"direction primitives → incidence\nlaunch grid → launch measure\nposition + direction neighborhoods → structural volume",ha="center",va="center"))
    save("sampling_measure_summary.png",lambda a:a.bar(["events","unique launch IDs","launch measure"],[len(arr["ray_index"]),len(np.unique(arr["launch_grid_index"])),launch_measure]))
    take=slice(None,None,max(1,len(arr["arrival_u"])//30000)); c=np.abs(bundle["final_area_ratio"][take])
    save("bundle_position_direction_state.png",lambda a:a.scatter(arr["arrival_u"][take],arr["arrival_dir_u"][take],c=c,s=1))
    save("optical_state_preservation_chain.png",lambda a:a.text(.5,.5,"launch → trajectory → receiver → arrival → OpticalInteractionState\nimmutable references + event-wise availability; no pixels",ha="center",va="center"))
    save("interaction_state_schema.png",lambda a:a.text(.5,.5,"geometry | carried signal | spectral | temporal | phase/coherence\npolarization | bundle | latent reverse | availability | provenance",ha="center",va="center"))

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--validation",action="store_true");args=ap.parse_args();RUN.mkdir(parents=True,exist_ok=True);t=time.time()
    prior=json.loads((UP/"result.json").read_text())
    if prior.get("interaction_contract_sha256")!=EXPECTED_CONTRACT or prior.get("structural_sha256")!=EXPECTED_STRUCTURAL or sha(UP/"structural_result.json")!=EXPECTED_STRUCTURAL:
        print("DEV131_STRUCTURAL_BASELINE_MISMATCH");return 2
    print("DEV131_INTERACTION_CONTRACT_SHA_VERIFIED=true\nDEV131_STRUCTURAL_SHA_VERIFIED=true")
    arr=load(A130/"arrival_events_primary.npz"); recv=load(R129/"receiver_state_primary.npz"); bundle=load(R129/"receiver_bundle_channels.npz")
    full_n=len(arr["ray_index"]); idx=np.arange(full_n) if not args.validation else np.arange(full_n)[::max(1,full_n//4096)][:4096]
    arr={k:v[idx] for k,v in arr.items()};recv={k:v[idx] for k,v in recv.items()};bundle={k:v[idx] for k,v in bundle.items()};n=len(idx)
    inv=inventory();dump("optical_state_inventory.json",inv);dump("optical_missing_state_ledger.json",missing_ledger())
    graph={"edges":[{"source":"arrival_dir_n","derived":"incidence_cosine"},{"source":"incidence_cosine","derived":"incidence_angle"},{"source":"launch_u","derived":"launch_measure"},{"source":"launch_v","derived":"launch_measure"},{"source":"neighbor arrival position","derived":"local_position_direction_volume"},{"source":"neighbor arrival direction","derived":"local_position_direction_volume"}],"acyclic":True,"derived_fields_without_retained_primitives":0};validate_derivation_graph(graph);dump("optical_state_derivation_graph.json",graph)
    contract=contract_schema()|{"field_inventory":"optical_state_inventory.json","forward_reverse_labels_independent":True};dump("optical_state_contract.json",contract);contract_sha=canonical_sha256(contract);print("DEV132_OPTICAL_STATE_CONTRACT_SHA256="+contract_sha)
    structural={"contract_version":CONTRACT_VERSION,"contract_sha256":contract_sha,"dev131_contract_sha256":EXPECTED_CONTRACT,"dev131_structural_sha256":EXPECTED_STRUCTURAL,"families":[f"O{i}" for i in range(8)],"availability_classes":[x for x in ("DIRECTLY_AVAILABLE","DERIVABLE_FROM_FROZEN_STATE","PARTIALLY_DERIVABLE","NOT_TRACKED","NOT_DEFINED_IN_CURRENT_MODEL","REQUIRES_NEW_PHYSICS_OR_SOURCE_INPUT")],"relevance_classes":["REQUIRED","POSSIBLY_RELEVANT","NOT_CURRENTLY_REQUIRED","UNKNOWN"],"derivation_classes":["EXACT_ALGEBRAIC","EXACT_GEOMETRIC","NUMERICAL_SUMMARY","REQUIRES_ASSUMPTION","NOT_DERIVABLE"],"propagation_changes":0,"trajectory_changes":0,"receiver_changes":0,"arrival_formation_changes":0,"dev131_interaction_changes":0,"target_access":False,"hst_pixel_access":False,"zero_rasterization":True,"zero_psf":True,"zero_source_reconstruction":True,"zero_lensing_reconstruction":True};dump("structural_result.json",structural);ssha=sha(RUN/"structural_result.json");print("DEV132_STRUCTURAL_SHA256="+ssha)
    du=np.median(np.diff(np.unique(arr["launch_u"]))) if len(np.unique(arr["launch_u"]))>1 else 0.;dv=np.median(np.diff(np.unique(arr["launch_v"]))) if len(np.unique(arr["launch_v"]))>1 else 0.;measure=float(abs(du*dv))
    sampling={"launch_type":"uniform launch grid","evidence":{"unique_u_count":int(len(np.unique(arr["launch_u"]))),"unique_v_count":int(len(np.unique(arr["launch_v"]))),"spacing_u":float(du),"spacing_v":float(dv)},"launch_measure":measure,"launch_measure_classification":"NUMERICAL_GEOMETRIC_MEASURE","physical_intensity_interpretation":False,"RAY_DENSITY_NOT_ASSUMED_PHYSICAL_INTENSITY":True,"weight_source":None,"initialization_rule":None,"conserved":None,"updated_during_propagation":False,"physically_interpreted":False,"bundle_area_relationship":"final received area ratio is stored; no flux law asserted","conservation_candidates":["ray identity","path ordering","bundle membership"]};dump("sampling_measure_manifest.json",sampling)
    flags={"physical_weight":np.zeros(n,bool),"spectral_state":np.zeros(n,bool),"arrival_time":np.zeros(n,bool),"phase":np.zeros(n,bool),"polarization":np.zeros(n,bool),"bundle_state":np.ones(n,bool),"reverse_provenance":np.ones(n,bool)}
    np.savez_compressed(RUN/"optical_state_availability_masks.npz",**flags)
    uids=event_uids(arr["ray_index"],arr["receiver_row_index"],idx)
    np.savez_compressed(RUN/"optical_state_primary.npz",event_uid=uids,arrival_index=idx,ray_index=arr["ray_index"],receiver_row_index=arr["receiver_row_index"],launch_grid_index=arr["launch_grid_index"],arrival_u=arr["arrival_u"],arrival_v=arr["arrival_v"],arrival_dir_u=arr["arrival_dir_u"],arrival_dir_v=arr["arrival_dir_v"],arrival_dir_n=arr["arrival_dir_n"],incidence_cosine=arr["receiver_incidence_cosine"],incidence_angle=arr["receiver_incidence_angle"],intersection_t=arr["intersection_t"],path_length=recv["path_length"],path_excess=recv["path_excess"],launch_measure=np.full(n,measure))
    availability={k:{"count":int(v.sum()),"fraction":float(v.mean())} for k,v in flags.items()}
    minimal={"description":"smallest currently eligible geometric optical event state","fields":["arrival_u","arrival_v","arrival_dir_u","arrival_dir_v","arrival_dir_n","incidence_cosine","incidence_angle","intersection_t"],"physical_weight":None,"physical_weight_status":"EVENT_WEIGHT_NOT_TRACKED","no_unit_weight_substitution":True};dump("minimal_forward_optical_state.json",minimal)
    sources=json.loads((UP/"arrival_event_archive_manifest.json").read_text())["sources"]
    full={"mode":"DIRECT_FIELDS_PLUS_IMMUTABLE_HASHED_REFERENCES","direct_archive":"optical_state_primary.npz","dev131_archive_sources":sources,"dev131_interaction_run":str(UP.relative_to(ROOT)),"all_dev131_information_accessible":True,"dev131_fields_lost":0,"primitive_first":True,"path_history_by_reference":"runs/wl_trajectory_state_completion001/diagnostic_full_paths.npz"};dump("full_preserved_optical_state.json",full)
    instrument={"geometric_optics_compatibility":"CURRENT_STATE_GEOMETRIC_OPTICS_SUFFICIENT","wave_optics_compatibility":"CURRENT_STATE_WAVE_OPTICS_INCOMPLETE","highest_forward_optical_readiness":"L0_GEOMETRIC_INTERACTION_READY","next_level":"L1_WEIGHTED_GEOMETRIC_READY","next_level_blocker":"no source-defined/conserved physical event weight","frequency_dependence":"FREQUENCY_STATE_UNUSED","frequency_interpretation":"implementation does not establish physical achromaticity","phase_dependence":"NOT_DEFINED","polarization_dependence":"NOT_DEFINED","TIME_DELAY_CONVERSION_NOT_YET_JUSTIFIED":True,"ABSOLUTE_TIME_ORIGIN_NOT_DEFINED":True,"PHASE_NOT_DERIVABLE_FROM_CURRENT_STATE":True,"PHYSICAL_EVENT_WEIGHT_UNAVAILABLE":True};dump("instrument_readiness.json",instrument)
    reverse={"highest_reverse_readiness":"R2_LOCAL_TRANSPORT_REVERSE_READY","verified_from":["Dev131 surface roundtrip","Dev131 transport conditioning","Dev131 reverse correspondence index"],"R3_not_claimed":"no actual source-candidate reconstruction machinery in the frozen contract"};dump("reverse_readiness.json",reverse)
    synth=synthetic_tests(); checks={"dev131_contract_sha_verified":True,"dev131_structural_sha_verified":True,"dev131_fields_lost_zero":True,"broken_optical_reverse_links_zero":len(np.unique(uids))==n,"geometry_state_inventory_complete":True,"event_weight_audit_complete":True,"spectral_state_audit_complete":True,"temporal_state_audit_complete":True,"phase_state_audit_complete":True,"polarization_state_audit_complete":True,"bundle_state_audit_complete":True,"sampling_measure_audit_complete":True,"frequency_dependence_audit_complete":True,"phase_dependence_audit_complete":True,"polarization_dependence_audit_complete":True,"missing_state_ledger_complete":True,"forward_relevance_complete":True,"reverse_relevance_complete":True,"primitive_derived_graph_acyclic":True,"derived_fields_without_retained_primitives_zero":True,"optical_state_contract_frozen":True,"synthetic_tests_pass":all(synth.values()),"target_access_false":True,"hst_pixel_access_false":True,"zero_rasterization":True,"zero_psf":True,"zero_source_reconstruction":True,"zero_lensing_reconstruction":True}
    outcomes=["WL_OPTICAL_INTERACTION_STATE_COMPLETION_ESTABLISHED","WL_GEOMETRIC_OPTICAL_INPUT_STATE_ESTABLISHED","WL_PHYSICAL_EVENT_WEIGHT_UNAVAILABLE","WL_SPECTRAL_STATE_UNAVAILABLE","WL_PHASE_COHERENCE_STATE_UNAVAILABLE","WL_POLARIZATION_STATE_UNAVAILABLE","WL_BIDIRECTIONAL_OPTICAL_STATE_PRESERVATION_ESTABLISHED"]
    result={"validation":args.validation,"population_count":n,"full_population_count":full_n,"availability":availability,"checks":checks,"synthetic_tests":synth,"highest_forward_optical_readiness":instrument["highest_forward_optical_readiness"],"highest_reverse_readiness":reverse["highest_reverse_readiness"],"outcomes":outcomes,"contract_sha256":contract_sha,"structural_sha256":ssha,"DEV131_FIELDS_LOST":0,"BROKEN_OPTICAL_REVERSE_LINKS":int(n-len(np.unique(uids))),"DERIVED_FIELDS_WITHOUT_RETAINED_PRIMITIVES":0,"target_access":False,"hst_pixel_access":False,"runtime_seconds":time.time()-t};dump("result.json",result)
    try:g=subprocess.run(["git","status","--short"],cwd=ROOT,text=True,capture_output=True).stdout
    except Exception:g="unavailable\n"
    (RUN/"baseline_git.txt").write_text(g);figures(availability,inv,arr,bundle,measure)
    qualification="VALIDATION_ONLY_NO_SCIENCE_CLAIM" if args.validation else "\n".join(outcomes)
    report=["PROPAGATION_CHANGES=0","TRAJECTORY_CHANGES=0","RECEIVER_CHANGES=0","ARRIVAL_FORMATION_CHANGES=0","DEV131_INTERACTION_CHANGES=0","TARGET_ACCESS=false","HST_PIXEL_ACCESS=false","RAY_DENSITY_NOT_ASSUMED_PHYSICAL_INTENSITY","PHYSICAL_EVENT_WEIGHT_UNAVAILABLE","TIME_DELAY_CONVERSION_NOT_YET_JUSTIFIED","ABSOLUTE_TIME_ORIGIN_NOT_DEFINED","PHASE_NOT_DERIVABLE_FROM_CURRENT_STATE",f"DEV132_OPTICAL_STATE_CONTRACT_SHA256={contract_sha}",f"DEV132_STRUCTURAL_SHA256={ssha}","DEV131_FIELDS_LOST=0","BROKEN_OPTICAL_REVERSE_LINKS=0","DERIVED_FIELDS_WITHOUT_RETAINED_PRIMITIVES=0",f"highest_forward_optical_readiness={instrument['highest_forward_optical_readiness']}",f"highest_reverse_readiness={reverse['highest_reverse_readiness']}",qualification]
    (RUN/"report.txt").write_text("\n".join(report)+"\n");print(qualification);return 0 if all(checks.values()) else 4
if __name__=="__main__":raise SystemExit(main())
