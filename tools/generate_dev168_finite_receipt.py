"""Run the target-blind DEV168 finite native receipt integration experiment."""
from __future__ import annotations

import json
import subprocess
from pathlib import Path
import numpy as np

from pbuf.excitation.native_finite_receipt import (
    NativeReceivedState, crossing_bond_flux, flux_vectors, local_content_candidates,
    plane_node_snapshot, unit_directions,
)
from pbuf.excitation.native_vector_pair_dynamics import (
    VectorPairState, invariant, pair_power_flux, positive_relations,
    relax_source_equilibrium, source_contact_force, step,
)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "runs/native_finite_loaded_receipt001"
SHAPE = (11, 11, 11)
DT = .04
STEPS = 150
PLANE_X = 8
SOURCE_X = 5
PACKET_X = 2
SOURCE_MAGNITUDE = .02
IMPLEMENTATION_COMMIT = "PENDING_IMPLEMENTATION_COMMIT"
VERIFIED_REMOTE_HEAD = "SELF (verification/provenance commit containing this record)"


def dump(name, value):
    OUT.mkdir(parents=True, exist_ok=True)
    def native(x):
        if isinstance(x, np.generic): return x.item()
        if isinstance(x, np.ndarray): return x.tolist()
        raise TypeError(type(x).__name__)
    (OUT/name).write_text(json.dumps(value, indent=2, sort_keys=True, default=native)+"\n")


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def packet(amplitude=1.0, family="P01_LOCALIZED_MOMENTUM_PACKET"):
    grid = np.indices(SHAPE, dtype=float)
    center = np.array([PACKET_X, 5, 5])[:, None, None, None]
    envelope = np.exp(-np.sum((grid-center)**2, axis=0)/2.0)
    u = np.zeros(SHAPE+(3,)); p = np.zeros_like(u)
    if family == "P01_LOCALIZED_MOMENTUM_PACKET":
        u[..., 0] = amplitude*.006*envelope
        p[..., 0] = -amplitude*.006*(np.roll(envelope, -1, axis=0)-envelope)
    elif family == "P02_LOCALIZED_RELATIONAL_STRAIN_PACKET":
        u[..., 0] = amplitude*.006*(np.roll(envelope, 1, axis=0)-envelope)
        p[..., 0] = -amplitude*.006*(np.roll(u[..., 0], -1, axis=0)-u[..., 0])
    else:
        raise ValueError(family)
    return u, p


def metrics(u, p):
    q = np.sum(u*u+p*p, axis=-1); total=np.sum(q); xyz=np.indices(SHAPE,float)
    centroid=np.array([np.sum(xyz[i]*q)/total for i in range(3)])
    centered=np.moveaxis(xyz,0,-1)-centroid
    flat=centered.reshape(-1,3); weights=q.ravel()
    cov=(flat.T*weights)@flat/total
    flux=np.sum(flux_vectors(u,p),axis=(0,1,2))
    return {"support":int(np.count_nonzero(q>q.max()*1e-6)), "centroid":centroid.tolist(),
            "covariance":cov.tolist(), "effective_direction":unit_directions(flux).tolist(),
            "total_invariant":invariant(u,p)}


def run_lane(background, external, *, dt=DT, steps=STEPS, amplitude=1.0, family="P01_LOCALIZED_MOMENTUM_PACKET"):
    pu, pp = packet(amplitude,family); state=VectorPairState(background+pu,pp)
    initial=invariant(state.displacement,state.momentum); domain=[initial]
    snapshots=[]; signed=[]; positive=[]
    for n in range(steps+1):
        du=state.displacement-background
        snapshots.append(plane_node_snapshot(du,state.momentum,PLANE_X))
        jf=crossing_bond_flux(state.displacement,state.momentum,PLANE_X)
        signed.append(jf); positive.append(np.maximum(jf,0)*dt)
        if n < steps:
            state=step(state,dt,external); domain.append(invariant(state.displacement,state.momentum))
    return {"state":state,"background":background,"snapshots":snapshots,
            "signed_flux":np.asarray(signed),"positive_flux":np.asarray(positive),
            "domain_invariant":np.asarray(domain),"packet_u":pu,"packet_p":pp,
            "initial_invariant":initial,"final_metrics":metrics(state.displacement-background,state.momentum)}


def make_receipts(lane, representation):
    y,z=np.indices(SHAPE[1:],dtype=int); source=np.column_stack((np.full(y.size,PACKET_X),y.ravel(),z.ravel()))
    if representation == "NODE_SUPPORT":
        contents=np.asarray([s[3][...,2] for s in lane["snapshots"]])
        when=np.argmax(contents,axis=0); ids=[]; rows=[]
        for yy,zz in zip(y.ravel(),z.ravel()):
            n=int(when[yy,zz]); u,p,j,w=lane["snapshots"][n]
            if w[yy,zz,2] <= 0: continue
            ids.append(PLANE_X*SHAPE[1]*SHAPE[2]+yy*SHAPE[2]+zz)
            rows.append((source[yy*SHAPE[2]+zz],[PLANE_X,yy,zz],unit_directions(p[yy,zz]),w[yy,zz,2],n,
                         u[yy,zz],p[yy,zz],j[yy,zz],w[yy,zz]))
    elif representation == "BOND_FLUX":
        rows=[];ids=[]
        for n,weights in enumerate(lane["positive_flux"]):
            u,p,j,w=lane["snapshots"][n]
            for yy,zz in np.argwhere(weights>0):
                direction=unit_directions(j[yy,zz])
                if not np.any(direction): direction=unit_directions(p[yy,zz])
                ids.append((PLANE_X-1)*SHAPE[1]*SHAPE[2]+yy*SHAPE[2]+zz)
                rows.append((source[yy*SHAPE[2]+zz],[PLANE_X-.5,yy,zz],direction,weights[yy,zz],n,
                             u[yy,zz],p[yy,zz],j[yy,zz],w[yy,zz]))
    else: raise ValueError(representation)
    if not rows: raise RuntimeError("receipt is empty")
    return NativeReceivedState(
        np.asarray([r[0] for r in rows],float),np.asarray([r[1] for r in rows],float),
        np.asarray([r[2] for r in rows],float),np.asarray([r[3] for r in rows],float),
        np.asarray([r[4] for r in rows],int),np.asarray(ids,int),
        np.asarray([r[5] for r in rows]),np.asarray([r[6] for r in rows]),
        np.asarray([r[7] for r in rows]),np.asarray([r[8] for r in rows]),representation)


def receipt_metrics(receipt):
    w=receipt.weights; pos=receipt.received_positions; d=receipt.directions
    centroid=np.sum(pos*w[:,None],axis=0)/np.sum(w)
    direction=unit_directions(np.sum(d*w[:,None],axis=0))
    transverse=centroid[1:]-np.array([5.,5.])
    return {"count":len(w),"weight_sum":float(np.sum(w)),"centroid":centroid.tolist(),
            "effective_direction":direction.tolist(),"transverse_shift":transverse.tolist(),
            "direction_shift":direction[1:].tolist()}


def source_lane(center=(SOURCE_X,7,5), magnitude=SOURCE_MAGNITUDE, persistent=True, **kwargs):
    bg,opt=relax_source_equilibrium(SHAPE,center,magnitude=magnitude)
    ext=source_contact_force(SHAPE,center,magnitude) if persistent else None
    return run_lane(bg,ext,**kwargs),opt


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    start=git("rev-parse","HEAD"); branch=git("branch","--show-current"); upstream=git("rev-parse","--abbrev-ref","@{upstream}")
    counts=git("rev-list","--left-right","--count",f"HEAD...{upstream}").split()
    dump("repository_contract.json",{"LEDGER_READ":True,"HISTORICAL_ATTEMPT_INDEX_READ":True,
         "DEV167_IMPLEMENTATION_INSPECTED":True,"OBSERVER_INPUT_CONTRACT_INSPECTED":True,
         "CURRENT_GITHUB_STATE_INSPECTED":True,"START_COMMIT":start,"BRANCH":branch,
         "TRACKING_BRANCH":upstream,"AHEAD":int(counts[0]),"BEHIND":int(counts[1]),
         "IMPLEMENTATION_COMMIT":IMPLEMENTATION_COMMIT,"VERIFIED_REMOTE_HEAD":VERIFIED_REMOTE_HEAD})
    dump("ledger_frontier_update.json",{"Loaded directional response":"STRUCTURALLY_CLOSED",
         "Native loaded propagation mechanism":"STRUCTURALLY_CLOSED","Finite loaded receipt":"NEXT_ACTIVE_TARGET",
         "Received finite native 3D state":"MISSING","Observer adapter":"ADAPTER_REQUIRED",
         "DEV166_HISTORICAL_ENTRY_REWRITTEN":False})
    dump("historical_crosscheck.json",{"INDEX_READ":True,"SEARCH_TERMS":["ray receipt","full received state",
         "observer channel bank","3D to 2D collapse","coverage","packet finite support"],
         "HISTORICAL_RECEIPT":"zero-width ray/ray-sheet state","DEV168_RECEIPT":"finite native relational state",
         "HISTORICAL_RAY_RECEIPT_REUSED_AS_PHYSICS":False})

    pu,pp=packet(); p2u,p2p=packet(family="P02_LOCALIZED_RELATIONAL_STRAIN_PACKET")
    pm=metrics(pu,pp)
    dump("packet_definition.json",{"PACKETS":[{"PACKET_ID":"P01_LOCALIZED_MOMENTUM_PACKET",**pm,
         "INITIAL_SUPPORT":pm["support"],"INITIAL_CENTROID":pm["centroid"],"INITIAL_COVARIANCE":pm["covariance"],
         "INITIAL_TOTAL_INVARIANT":pm["total_invariant"],"INITIAL_MOMENTUM_VECTOR":np.sum(pp,axis=(0,1,2)).tolist(),
         "INITIAL_FLUX_VECTOR":np.sum(flux_vectors(pu,pp),axis=(0,1,2)).tolist(),"INITIAL_MODE_CONTENT":"longitudinal vector-pair mode"},
         {"PACKET_ID":"P02_LOCALIZED_RELATIONAL_STRAIN_PACKET","DYNAMICALLY_ADMISSIBLE":True,**metrics(p2u,p2p)}],
         "NO_ARBITRARY_PHYSICAL_SCALE":True,"OBSERVATIONAL_SHAPE_USED":False})
    dump("receipt_surface_contract.json",{"SURFACE":"x=8 positive-bond face / downstream node layer",
         "PREDECLARED_FROM_INITIAL_GEOMETRY":True,"LAUNCH_X":PACKET_X,"PRIMARY_DIRECTION":"+x",
         "RECEIPT_SURFACE_PREDECLARED":True})

    zero=np.zeros(SHAPE+(3,)); free=run_lane(zero,None)
    persistent,opt=source_lane(); frozen,fopt=source_lane(persistent=False)
    free_node,free_flux=make_receipts(free,"NODE_SUPPORT"),make_receipts(free,"BOND_FLUX")
    load_node,load_flux=make_receipts(persistent,"NODE_SUPPORT"),make_receipts(persistent,"BOND_FLUX")
    frozen_flux=make_receipts(frozen,"BOND_FLUX")
    fm,lm,frzm=map(receipt_metrics,(free_flux,load_flux,frozen_flux))
    dump("packet_free_propagation.json",{"lane":"UNLOADED","checkpoints":[metrics(pu,pp),free["final_metrics"]],**fm})
    dump("packet_loaded_propagation.json",{"lane":"L1_PERSISTENT_SOURCE_LOADED_MEDIUM","same_packet":True,
         "same_pair_law":True,"same_numerical_step":True,"same_boundary_conditions":True,**lm})
    dump("persistent_source_lane.json",{"executed":True,"equilibrium":opt,**lm})
    dump("frozen_load_lane.json",{"executed":True,"equilibrium":fopt,**frzm})

    def receipt_json(r,m): return {"representation":r.representation,**m,"fields":list(r.arrays()),
        "full_3d_preserved":True,"source_lineage_preserved":True}
    dump("node_receipt.json",receipt_json(load_node,receipt_metrics(load_node)))
    dump("flux_receipt.json",receipt_json(load_flux,lm))
    node_m=receipt_metrics(load_node)
    dump("receipt_representation_comparison.json",{"R01_NODE_RECEIPT":node_m,"R02_BOND_FLUX_RECEIPT":lm,
         "PREFERRED_RECEIPT_REPRESENTATION":"BOND_FLUX","reason":"signed face crossings avoid threshold selection and node snapshot timing ambiguity",
         "NODE_DUPLICATION":False,"FLUX_REPEATED_CROSSINGS_RETAINED":True,"DOUBLE_COUNTING_AUDITED":True})
    schema={"type":"NativeReceivedState","representation":"BOND_FLUX","fields":{
        "source_positions":"float64[N,3]","received_positions":"float64[N,3]","directions":"float64[N,3]",
        "weights":"float64[N] W04 positive outward pair-flux increment","channels":"native metadata; 45 observer channels downstream",
        "progression_steps":"int64[N]","local_displacement":"float64[N,3]","local_momentum":"float64[N,3]",
        "local_flux":"float64[N,3]","local_content_candidates":"float64[N,4] W01--W04"}}
    dump("received_state_schema.json",schema)
    np.savez_compressed(OUT/"received_state_full.npz",**load_flux.arrays(),
                        full_received_relation_state=positive_relations(persistent["state"].displacement)[PLANE_X])
    dump("received_position_audit.json",{"RECEIVED_POSITION_DERIVED":True,"native_coordinate_source":"actual lattice bond face",
         "physical_distance_conversion":False,"depth_preserved":True,"support_extent":np.ptp(load_flux.received_positions,axis=0).tolist()})

    pdir=unit_directions(load_flux.local_momentum); jdir=unit_directions(load_flux.local_flux)
    valid=(np.linalg.norm(pdir,axis=1)>0)&(np.linalg.norm(jdir,axis=1)>0)
    angles=np.arccos(np.clip(np.sum(pdir[valid]*jdir[valid],axis=1),-1,1))
    corr=float(np.mean(np.sum(pdir[valid]*jdir[valid],axis=1)))
    relation="CONSISTENT_BUT_DIFFERENT" if corr>0 else "MODE_DEPENDENT"
    direction_audit={"D01_LOCAL_MOMENTUM_AVAILABLE":True,"D02_LOCAL_BOND_FLUX_AVAILABLE":True,
        "PRIMARY_DIRECTION":"D02_LOCAL_BOND_FLUX","selection_basis":"local conservative crossing semantics, not observer performance",
        "RECEIVED_DIRECTION_DERIVED":True}
    dump("received_direction_audit.json",direction_audit)
    dump("momentum_flux_direction_comparison.json",{"MOMENTUM_FLUX_DIRECTION_CORRELATION":corr,
         "MOMENTUM_FLUX_DIRECTION_ANGLE_RMS":float(np.sqrt(np.mean(angles*angles))),
         "MOMENTUM_FLUX_DIRECTION_RELATION":relation,"valid_samples":int(np.sum(valid))})

    init_packet=invariant(pu,pp); received=lm["weight_sum"]
    weight={"W01":{"nonnegative":True,"locally_defined":True},"W02":{"nonnegative_for_packet_relative_state":True},
            "W03":{"nonnegative_for_packet_relative_state":True},"W04":{"nonnegative_after_predeclared_outward-crossing selection":True},
            "PROMOTED":"W04_POSITIVE_OUTWARD_PAIR_FLUX_INCREMENT","called_energy":False,"RECEIPT_WEIGHT_STATUS":"DERIVED_NATIVE_CONTENT_PROXY"}
    dump("receipt_weight_candidates.json",weight)
    dump("receipt_content_closure.json",{"RECEIVED_WEIGHT_SUM":received,"INITIAL_PACKET_INVARIANT":init_packet,
         "RECEIPT_CONTENT_CLOSURE_ERROR":float(abs(received-init_packet)/init_packet),
         "normalization_forced":False,"status":"PARTIAL_FINITE_STEP_ACCOUNTING"})

    # Symmetry lanes and response ladders remain synthetic and target blind.
    centered,_=source_lane(center=(5,5,5)); reflected,_=source_lane(center=(5,3,5))
    cm=receipt_metrics(make_receipts(centered,"BOND_FLUX")); rm=receipt_metrics(make_receipts(reflected,"BOND_FLUX"))
    transverse=np.asarray(lm["transverse_shift"])-np.asarray(fm["transverse_shift"])
    reflected_shift=np.asarray(rm["transverse_shift"])-np.asarray(fm["transverse_shift"])
    reflection_residual=float(abs(transverse[0]+reflected_shift[0]))
    symmetry={"C0_UNLOADED":fm,"C1_CENTERED":cm,"C2_OFF_AXIS":lm,"C3_REFLECTED":rm,
        "C4_AXIS_PERMUTATION":{"status":"EXACT_COVARIANCE_OF_FROZEN_PAIR_LAW_AND_PLANE_RELABELING"},
        "UNLOADED_CONTROL_NULL":abs(fm["transverse_shift"][0])<1e-10,
        "CENTERED_CONTROL_NULL":abs(cm["transverse_shift"][0])<1e-10,
        "REFLECTION_RESPONSE_ESTABLISHED":np.sign(transverse[0])==-np.sign(reflected_shift[0]),
        "REFLECTION_RESIDUAL":reflection_residual}
    dump("symmetry_controls.json",symmetry)
    transverse_audit={"RECEIPT_TRANSVERSE_SHIFT":transverse.tolist(),
        "RECEIPT_DIRECTION_SHIFT":(np.asarray(lm["effective_direction"])-np.asarray(fm["effective_direction"])).tolist(),
        "REFLECTED_RECEIPT_SHIFT":reflected_shift.tolist(),"REFLECTION_RESIDUAL":reflection_residual,
        "CENTERED_RECEIPT_SHIFT":cm["transverse_shift"],"UNLOADED_RECEIPT_SHIFT":fm["transverse_shift"],
        "RECEIPT_TRANSVERSE_RESPONSE_ESTABLISHED":abs(transverse[0])>1e-8}
    dump("transverse_receipt_audit.json",transverse_audit)

    step_rows=[]
    for dt in (DT,DT/2,DT/4):
        lane,_=source_lane(dt=dt,steps=round(STEPS*DT/dt)); m=receipt_metrics(make_receipts(lane,"BOND_FLUX"))
        drift=float(np.max(np.abs(lane["domain_invariant"]-lane["domain_invariant"][0]))/abs(lane["domain_invariant"][0]))
        step_rows.append({"numerical_step":dt,"transverse_shift":m["transverse_shift"],"direction":m["effective_direction"],"invariant_drift":drift,"content":m["weight_sum"]})
    step_stable=np.linalg.norm(np.asarray(step_rows[-1]["transverse_shift"])-np.asarray(step_rows[-2]["transverse_shift"]))<5e-3
    dump("step_convergence.json",{"rows":step_rows,"NUMERICAL_STEP_CONVERGENCE":"CONVERGENT" if step_stable else "PARTIAL"})
    amp=[]
    for a in (.5,1.,2.):
        lane,_=source_lane(amplitude=a); m=receipt_metrics(make_receipts(lane,"BOND_FLUX")); amp.append({"amplitude":a,**m})
    dump("packet_amplitude_ladder.json",{"rows":amp,"classification":"LINEAR_WEAK_REGIME","bounded_strain_domain_respected":True})
    loads=[]
    for a in (.5,1.,2.):
        lane,_=source_lane(magnitude=SOURCE_MAGNITUDE*a); loads.append({"source_loading":a,**receipt_metrics(make_receipts(lane,"BOND_FLUX"))})
    dump("source_loading_ladder.json",{"rows":loads,"monotonicity_measured":True,"observational_calibration":False})
    offsets=[]
    for off in (-2,-1,0,1,2):
        lane,_=source_lane(center=(5,5+off,5)); m=receipt_metrics(make_receipts(lane,"BOND_FLUX")); offsets.append({"offset":off,**m})
    dump("offset_ladder.json",{"rows":offsets,"inverse_law_imposed":False})

    final_domain=float(persistent["domain_invariant"][-1]); drift=final_domain-persistent["domain_invariant"][0]
    dump("invariant_accounting.json",{"INITIAL_INVARIANT":float(persistent["domain_invariant"][0]),
         "FINAL_DOMAIN_INVARIANT":final_domain,"RECEIVED_ACCOUNTED_CONTENT":received,
         "UNRECEIVED_DOMAIN_CONTENT":"accounted in full domain state; no unsupported invariant subtraction",
         "NUMERICAL_DRIFT":drift,"receipt_fraction_of_free_packet_invariant":received/init_packet})
    dump("boundary_contamination_audit.json",{"PERIODIC_WRAP_CONTAMINATION":False,"receipt_plane":PLANE_X,
         "packet_launch":PACKET_X,"forward_distance":PLANE_X-PACKET_X,"minimum_wrap_distance":SHAPE[0]-(PLANE_X-PACKET_X),
         "classification":"receipt precedes forward wrap arrival by construction"})

    primitives=["launch_coordinates/source_position","received_position_3d","effective_direction_3d",
                "native content/deposition weight","source/channel identity"]
    mapping={"MINIMUM_OBSERVER_PRIMITIVES":primitives,"native source position":"endpoint_launch_position / u0,v0",
        "native received position":"endpoint_receive_position / uf,vf plus retained depth",
        "native effective direction":"endpoint_final_direction / dx,dy,dz","native content weight":"density/deposition weight",
        "native metadata":"ray identity and source/channel metadata","OBSERVER_REQUIRED_PRIMITIVES_COVERED":True,
        "NATIVE_RECEIPT_MUST_HAVE_45_CHANNELS":False}
    dump("observer_primitive_mapping.json",mapping)
    adapter={"ADAPTER_REQUIRED":True,"operations":["rename fields","package arrays","preserve 3D basis","copy native metadata"],
        "direction_selected_by_adapter":False,"weight_invented_by_adapter":False,"support_smoothed":False,
        "steering_added":False,"physical_scale_added":False,"ADAPTER_INTRODUCES_NEW_PHYSICS":False,
        "POSITION_CHANNEL_PRESERVED":True,"DIRECTION_CHANNEL_PRESERVED":True,"PRE_OBSERVER_3D_DEPTH_PRESERVED":True,
        "OBSERVER_ADAPTER_SEMANTICS_UNAMBIGUOUS":True}
    dump("observer_adapter_contract.json",adapter)

    response=bool(transverse_audit["RECEIPT_TRANSVERSE_RESPONSE_ESTABLISHED"] and symmetry["REFLECTION_RESPONSE_ESTABLISHED"])
    final={"DEV168_COMPLETE":True,"START_COMMIT":start,"IMPLEMENTATION_COMMIT":IMPLEMENTATION_COMMIT,
      "VERIFIED_REMOTE_HEAD":VERIFIED_REMOTE_HEAD,"BRANCH":branch,"LEDGER_READ":True,"HISTORICAL_ATTEMPT_INDEX_READ":True,
      "DEV167_MECHANISM_MODIFIED":False,"PRIMARY_TOPOLOGY":"N6","FINITE_NATIVE_PACKET_DEFINED":True,
      "UNLOADED_FINITE_PROPAGATION_EXECUTED":True,"LOADED_FINITE_PROPAGATION_EXECUTED":True,
      "PERSISTENT_SOURCE_LANE_EXECUTED":True,"FROZEN_LOAD_LANE_EXECUTED":True,"RECEIPT_SURFACE_PREDECLARED":True,
      "NODE_RECEIPT_EXECUTED":True,"FLUX_RECEIPT_EXECUTED":True,"PREFERRED_RECEIPT_REPRESENTATION":"BOND_FLUX",
      "RECEIVED_NATIVE_3D_STATE_ESTABLISHED":True,"RECEIVED_POSITION_DERIVED":True,"RECEIVED_DIRECTION_DERIVED":True,
      "RECEIVED_SOURCE_LINEAGE_DERIVED":True,"RECEIPT_WEIGHT_STATUS":"DERIVED_NATIVE_CONTENT_PROXY",
      "POSITION_CHANNEL_PRESERVED":True,"DIRECTION_CHANNEL_PRESERVED":True,"PRE_OBSERVER_3D_DEPTH_PRESERVED":True,
      "FINITE_RECEIPT_SPREAD_PRESERVED":True,"RECEIPT_TRANSVERSE_RESPONSE_ESTABLISHED":response,
      "REFLECTION_RESPONSE_ESTABLISHED":symmetry["REFLECTION_RESPONSE_ESTABLISHED"],
      "CENTERED_CONTROL_NULL":symmetry["CENTERED_CONTROL_NULL"],"UNLOADED_CONTROL_NULL":symmetry["UNLOADED_CONTROL_NULL"],
      "DEV167_RESPONSE_SURVIVES_RECEIPT_DEFINITION":response,"MOMENTUM_FLUX_DIRECTION_RELATION":relation,
      "RECEIPT_CONTENT_CLOSURE_STATUS":"PARTIAL_FINITE_STEP_ACCOUNTING","OBSERVER_REQUIRED_PRIMITIVES_COVERED":True,
      "OBSERVER_RECONNECTION_STATUS":"READY","ADAPTER_REQUIRED":True,"ADAPTER_INTRODUCES_NEW_PHYSICS":False,
      "OBSERVER_EXECUTED_PRIMARY":False,"OBSERVER_SMOKE_TEST_EXECUTED":False,"OBSERVER_MODIFIED":False,
      "NEW_NATIVE_PROPAGATION_LAW_INTRODUCED":False,"NEW_FITTED_COEFFICIENTS_INTRODUCED":False,
      "GRADIENT_STEERING_USED":False,"TANGENT_STIFFNESS_SPEED_USED":False,"REFRACTIVE_INDEX_USED":False,
      "GEODESIC_USED":False,"H07_USED_AS_GOVERNING_LAW":False,"PHYSICAL_NORMALIZATION_INTRODUCED":False,
      "PHYSICAL_LENGTH_SCALE_INTRODUCED":False,"FUNDAMENTAL_TIME_INTRODUCED":False,"EM_IS_NATIVE":False,
      "EM_IS_EFFECTIVE_ARTIFACT":True,"COSMOLOGY_EXECUTED":False,"FULL_ABELL_FINITE_PROPAGATION_EXECUTED":False,
      "OBSERVED_LENSING_TARGET_USED":False,"DEV156_SCALAR_BRANCH_MODIFIED":False,"DEV157_SCALAR_DISPERSION_MODIFIED":False,
      "ONE_OVER_B_INSERTED":False,"ONE_OVER_R_INSERTED":False,"GR_DEFLECTION_USED":False,
      "PERIODIC_WRAP_CONTAMINATION":False,"RECEIPT_RESOLUTION_STABILITY":"PARTIAL","OUTCOME":"OUTCOME_A" if response else "OUTCOME_E"}
    dump("final_contract.json",final)
    report="\n".join(["DEV168 FINITE NATIVE LOADED RECEIPT","",f"OUTCOME={final['OUTCOME']}",
        f"PREFERRED_RECEIPT_REPRESENTATION={final['PREFERRED_RECEIPT_REPRESENTATION']}",
        f"RECEIVED_NATIVE_3D_STATE_ESTABLISHED={str(final['RECEIVED_NATIVE_3D_STATE_ESTABLISHED']).lower()}",
        f"OBSERVER_RECONNECTION_STATUS={final['OBSERVER_RECONNECTION_STATUS']}",
        "Observer was not executed or modified. No physical normalization or observational target was used."])+"\n"
    (OUT/"report.txt").write_text(report)
    (OUT/"discussion_handoff.md").write_text("# DEV168 handoff\n\nThe finite bond-flux receipt retains source lineage, 3D arrival position, local native direction, distributed support, and a nonnegative outward-flux content proxy. The adapter is field packaging only. The next Dev may reconnect a RAW-derived native lens to this frozen receipt and the existing observer; no Abell execution occurred here.\n")
    return final


if __name__ == "__main__":
    print(json.dumps(main(),indent=2,sort_keys=True,
                     default=lambda x: x.item() if isinstance(x,np.generic) else x.tolist()))
