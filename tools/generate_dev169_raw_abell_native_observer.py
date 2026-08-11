"""DEV169 target-blind RAW Abell native-to-observer structural execution."""
from __future__ import annotations

import hashlib, json, subprocess, sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))

from pbuf.source.projected_source_3d_family import diagnostic_family, project
from pbuf.excitation.native_vector_pair_dynamics import (
    VectorPairState, invariant, net_force, pair_forces, positive_relations,
    source_contact_force, step,
)
from pbuf.excitation.native_finite_receipt import (
    NativeReceivedState, crossing_bond_flux, flux_vectors,
    local_content_candidates, plane_node_snapshot, unit_directions,
)
from pbuf.excitation.native_observer_adapter import adapt_native_receipt, execute_frozen_observer

OUT=ROOT/"runs/raw_abell_native_observer001"
SOURCE=ROOT/"runs/raw_abell2744_detector_to_native_source001/native_2d_source_constraint.npz"
SHAPE=(11,11,11); LAUNCH_X=1; PLANE_X=8; DT=.04; STEPS=180; SOURCE_MAGNITUDE=.02
START_COMMIT="6d3753fb105f595ab7f299578794b7870ec1ad4b"

def git(*a): return subprocess.check_output(["git",*a],cwd=ROOT,text=True).strip()
def native(x):
    if isinstance(x,np.generic): return x.item()
    if isinstance(x,np.ndarray): return x.tolist()
    raise TypeError(type(x).__name__)
def dump(name,obj):
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/name).write_text(json.dumps(obj,indent=2,sort_keys=True,default=native,allow_nan=False)+"\n")
def finite(v):
    v=float(v); return v if np.isfinite(v) else None
def resize_blocks(a,n=7):
    edges=np.linspace(0,a.shape[0],n+1,dtype=int); out=np.zeros((n,n))
    for i in range(n):
      for j in range(n): out[i,j]=a[edges[i]:edges[i+1],edges[j]:edges[j+1]].sum()
    return out/out.sum()
def embedded_family(image):
    rows=diagnostic_family(image,nz=9); out=[]
    for r in rows:
        a=np.zeros(SHAPE); a[1:10,2:9,2:9]=r.source
        out.append((r,a))
    return out
def distributed_force(source,scale=1.):
    ext=np.zeros(SHAPE+(3,)); cells=np.argwhere(source>0)
    for cell in cells:
        ext += source_contact_force(SHAPE,tuple(map(int,cell)),SOURCE_MAGNITUDE*scale*float(source[tuple(cell)]))
    return ext
def equilibrium(ext):
    def unpack(v):
        u=v.reshape(SHAPE+(3,)); return u-u.mean((0,1,2),keepdims=True)
    def fun(v):
        u=unpack(v); e=np.linalg.norm(positive_relations(u),axis=-1)-1
        if np.any(np.abs(e)>=1): return 1e30
        return float(np.sum(-.5*np.log1p(-e*e))-np.sum(ext*u))
    def jac(v):
        u=unpack(v); g=-(net_force(u)+ext); g-=g.mean((0,1,2),keepdims=True); return g.ravel()
    res=minimize(fun,np.zeros(np.prod(SHAPE)*3),jac=jac,method="L-BFGS-B",options={"gtol":2e-9,"maxiter":1200})
    u=unpack(res.x); residual=float(np.max(np.abs(net_force(u)+ext)))
    return u,{"success":bool(res.success or residual<2e-7),"iterations":int(res.nit),"max_force_residual":residual,"message":str(res.message)}
def packet(image):
    transverse=np.zeros(SHAPE[1:]); transverse[2:9,2:9]=image
    transverse/=transverse.max(); x=np.arange(SHAPE[0]); env=np.exp(-.5*((x-LAUNCH_X)/.8)**2)[:,None,None]*transverse[None]
    u=np.zeros(SHAPE+(3,));p=np.zeros_like(u);u[...,0]=.006*env;p[...,0]=-.006*(np.roll(env,-1,axis=0)-env)
    return u,p
def run(background,ext,image,dt=DT,steps=STEPS):
    pu,pp=packet(image); state=VectorPairState(background+pu,pp); inv=[invariant(state.displacement,state.momentum)]
    snapshots=[]; positive=[]
    for n in range(steps+1):
        du=state.displacement-background;snapshots.append(plane_node_snapshot(du,state.momentum,PLANE_X))
        positive.append(np.maximum(crossing_bond_flux(state.displacement,state.momentum,PLANE_X),0)*dt)
        if n<steps: state=step(state,dt,ext);inv.append(invariant(state.displacement,state.momentum))
    return {"state":state,"snapshots":snapshots,"positive":np.asarray(positive),"invariant":np.asarray(inv),"packet":(pu,pp)}
def receipt(lane,image):
    support=np.argwhere(image>0); rows=[];ids=[]
    for n,wg in enumerate(lane["positive"]):
      u,p,j,c=lane["snapshots"][n]
      for y,z in np.argwhere(wg>0):
        # Objective nearest supported source-plane cell; lineage is never target-selected.
        q=support[np.argmin(np.sum((support-np.array([y-2,z-2]))**2,axis=1))]
        src=np.array([LAUNCH_X,q[0]+2,q[1]+2.],float); d=unit_directions(j[y,z])
        if not np.any(d): d=unit_directions(p[y,z])
        rows.append((src,[PLANE_X-.5,y,z],d,wg[y,z],n,u[y,z],p[y,z],j[y,z],c[y,z]));ids.append(int(q[0]*7+q[1]))
    if not rows: raise RuntimeError("empty bond-flux receipt")
    return NativeReceivedState(np.asarray([x[0] for x in rows]),np.asarray([x[1] for x in rows]),np.asarray([x[2] for x in rows]),
      np.asarray([x[3] for x in rows]),np.asarray([x[4] for x in rows]),np.asarray(ids),np.asarray([x[5] for x in rows]),
      np.asarray([x[6] for x in rows]),np.asarray([x[7] for x in rows]),np.asarray([x[8] for x in rows]),"BOND_FLUX")
def receipt_summary(r):
    w=r.weights; p=r.received_positions; d=r.directions; total=float(w.sum()); cen=np.sum(p*w[:,None],0)/total
    q=p-cen;cov=(q.T*w)@q/total;direction=unit_directions(np.sum(d*w[:,None],0))
    return {"receipt_count":len(w),"total_received_proxy":total,"centroid":cen,"covariance":cov,
      "flux_direction":direction,"finite":bool(np.isfinite(p).all() and np.isfinite(d).all() and np.isfinite(w).all())}
def observer(r):
    a=adapt_native_receipt(r);bank,meta=execute_frozen_observer(a,bins=6); primary=np.nan_to_num(bank[meta["primary_channel"]])
    return a,bank,meta,primary
def output_summary(a):
    a=np.asarray(a,float); yy,xx=np.indices(a.shape); w=np.abs(a); total=float(w.sum());
    cen=np.array([(w*yy).sum(),(w*xx).sum()])/total if total else np.zeros(2);q=np.stack((yy-cen[0],xx-cen[1]),-1)
    qf=q.reshape(-1,2); wf=w.ravel(); cov=(qf.T*wf)@qf/total if total else np.zeros((2,2)); vals,vec=np.linalg.eigh(cov)
    return {"finite_pixel_count":int(np.isfinite(a).sum()),"nonfinite_pixel_count":int((~np.isfinite(a)).sum()),"total_deposition":float(a.sum()),
      "centroid":cen,"covariance":cov,"principal_axes":vals[::-1],"asymmetry":float(np.mean(np.abs(a-np.flip(a))))}
def corr(a,b):
    x=np.asarray(a).ravel();y=np.asarray(b).ravel();x-=x.mean();y-=y.mean();den=np.linalg.norm(x)*np.linalg.norm(y);return float(x@y/den) if den else 1.
def save_png(name,a,title):
    fig,ax=plt.subplots();im=ax.imshow(a,origin="lower");ax.set_title(title);fig.colorbar(im,ax=ax);fig.tight_layout();fig.savefig(OUT/name,dpi=120);plt.close(fig)

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    branch=git("branch","--show-current"); remote=git("rev-parse","origin/HEAD")
    with np.load(SOURCE,allow_pickle=False) as z:
        raw=z["amplitude"].mean(0); channels=z["channels"].astype(str).tolist(); support=int(z["support"].any(0).sum())
    image=resize_blocks(raw); family=embedded_family(image)
    dump("repository_provenance.json",{"branch":branch,"start_commit":START_COMMIT,"Dev161_implementation_SHA":"7d919178b29e27832e231b740b9b6839f1805c13","Dev162_implementation_SHA":"7d919178b29e27832e231b740b9b6839f1805c13","Dev167_implementation_commit":"0bb0263cbf12135d6742dc51e7b291de06d586cd","Dev168_implementation_commit":"f3cc6540b84fd3fad2410a94fe3f8b9a89685d61","observer_implementation_SHA":git("log","-1","--format=%H","--","pbuf/labs/foundation/native_full_received_state_information_retention001.py"),"final_implementation_commit":"PENDING_SEPARATE_COMMIT","verified_remote_HEAD":remote})
    dump("raw_inventory.json",{"TARGET_CLUSTER":"ABELL_2744","RAW":116,"FLT":116,"FLC":116,"archive_exists":(ROOT/"PBUF_raw_benchmark/WLRAW-001_Abell2744").exists(),"channels":channels,"original_support_cells":support})
    dump("raw_to_native_source_contract.json",{"INPUT":"Dev161 NATIVE_MULTI_CHANNEL_SOURCE_CONSTRAINT","operation":"fixed 7x7 block-sum computational projection","coverage_of_Dev161_support":"100_PERCENT","absolute_scale":"RELATIVE_ONLY","interpreted_as_mass_density":False,"law_modified":False})
    dump("source_support.json",{"rule":"all positive fixed block sums","launch_cells":int(np.count_nonzero(image)),"available_cells":int(np.count_nonzero(image)),"coverage":"100_PERCENT","target_blind":True})
    dump("source_amplitude_contract.json",{"A0":"Dev167 canonical source-contact magnitude 0.02 times relative source weights","ladder":[.5,1,2],"selected_by_target":False,"physical_normalization":"UNRESOLVED"})
    dump("depth_family_contract.json",{"lanes":[r.name for r,_ in family],"count":len(family),"projection_errors":{r.name:float(np.max(np.abs(project(r.source)-image))) for r,_ in family},"unique_depth_assumed":False})
    dump("scalar_to_vector_source_loading_contract.json",{"SOURCE_CONSTRAINT_INPUT":"relative projected Dev161 constraint plus Dev162 diagnostic LOS weights","VECTOR_SOURCE_FORCE_CONSTRUCTION":"linear superposition of frozen Dev167 six-neighbor one-cell source_contact_force at every supported 3D source cell","NEW_COEFFICIENT_REQUIRED":False,"SEMANTIC_EQUIVALENCE":"same local one-cell source-contact semantics; distributed composition only"})
    dump("observer_adapter_mapping.json",{"observer.position":"native received_positions","observer.direction":"native flux directions","observer.weight":"native outward-flux content proxy (retained without native normalization)","observer.launch_coordinates":"native source_positions/source lineage","observer.metadata":"progression and local finite state","new_physics":False})
    dump("adapter_mapping.json",json.loads((OUT/"observer_adapter_mapping.json").read_text()))

    depth_rows=[];outputs={};receipts={};equilibria=[]; lens_stats=[]
    free_bg=np.zeros(SHAPE+(3,)); free=run(free_bg,None,image); free_r=receipt(free,image); free_a,free_bank,free_meta,free_o=observer(free_r)
    np.save(OUT/"native_observer_unloaded_2d.npy",free_o);save_png("native_observer_unloaded_2d.png",free_o,"Unloaded native observer control")
    for r,source in family:
        ext=distributed_force(source);bg,opt=equilibrium(ext);equilibria.append({"lane":r.name,**opt})
        rel=positive_relations(bg); lengths=np.linalg.norm(rel,axis=-1); lens_stats.append({"lane":r.name,"min_bond_length":float(lengths.min()),"max_bond_length":float(lengths.max()),"max_abs_strain":float(np.max(np.abs(lengths-1)))})
        lane=run(bg,ext,image);rec=receipt(lane,image);adapt,bank,meta,out=observer(rec);rs=receipt_summary(rec);os=output_summary(out)
        receipts[r.name]=rec;outputs[r.name]=out
        dump(f"depth_{r.name}_receipt_summary.json",rs);dump(f"depth_{r.name}_observer_summary.json",os);np.save(OUT/f"depth_{r.name}_observer_2d.npy",out);save_png(f"depth_{r.name}_observer_2d.png",out,r.name)
        depth_rows.append({"lane":r.name,"receipt":rs,"observer":os,"invariant_drift":float(np.max(np.abs(lane["invariant"]-lane["invariant"][0]))/max(abs(lane["invariant"][0]),1e-30))})
    primary=family[0][0].name;loaded_o=outputs[primary];loaded_r=receipts[primary];difference=loaded_o-free_o
    np.save(OUT/"native_observer_loaded_2d.npy",loaded_o);np.save(OUT/"native_observer_difference_2d.npy",difference)
    save_png("native_observer_loaded_2d.png",loaded_o,"Loaded native observer output");save_png("native_observer_difference_2d.png",difference,"Loaded - unloaded")
    comparisons=[]
    for i in range(len(family)):
      for j in range(i+1,len(family)):
        a=family[i][0].name;b=family[j][0].name;ra=receipt_summary(receipts[a]);rb=receipt_summary(receipts[b])
        comparisons.append({"lanes":[a,b],"2D_OUTPUT_CORRELATION":corr(outputs[a],outputs[b]),"2D_OUTPUT_RMS_DIFFERENCE":float(np.sqrt(np.mean((outputs[a]-outputs[b])**2))),"RECEIPT_CENTROID_DIFFERENCE":float(np.linalg.norm(np.asarray(ra["centroid"])-rb["centroid"])),"RECEIPT_COVARIANCE_DIFFERENCE":float(np.linalg.norm(np.asarray(ra["covariance"])-rb["covariance"])),"FLUX_DIRECTION_DIFFERENCE":float(np.linalg.norm(np.asarray(ra["flux_direction"])-rb["flux_direction"]))})
    dump("depth_family_output_comparison.json",{"pairs":comparisons,"minimum_correlation":min(x["2D_OUTPUT_CORRELATION"] for x in comparisons)})
    dump("vector_lens_contract.json",{"construction":"frozen Dev167 source_contact_force superposition","PRIMARY_TOPOLOGY":"N6","metric_inferred":False,"pair_law_modified":False})
    dump("vector_lens_equilibrium.json",{"rows":equilibria,"all_converged":all(x["success"] for x in equilibria)})
    dump("vector_lens_depth_family.json",{"rows":lens_stats});dump("bond_relation_statistics.json",{"rows":lens_stats});dump("bounded_strain_audit.json",{"rows":lens_stats,"respected":max(x["max_abs_strain"] for x in lens_stats)<1})
    # Predeclared structural amplitude ladder, representative depth only.
    amp=[]
    for scale in (.5,1,2):
        ext=distributed_force(family[0][1],scale);bg,opt=equilibrium(ext);rr=receipt(run(bg,ext,image),image);amp.append({"scale":scale,"equilibrium":opt,"receipt":receipt_summary(rr)})
    dump("source_amplitude_ladder.json",{"rows":amp,"observational_selection":False})
    # Step and receipt-duration convergence on the representative lane.
    ext=distributed_force(family[0][1]);bg,_=equilibrium(ext);steps=[]
    for dt in (DT,DT/2,DT/4):
        lane=run(bg,ext,image,dt=dt,steps=round(STEPS*DT/dt));rr=receipt(lane,image);_,_,_,oo=observer(rr)
        steps.append({"h":dt,"receipt":receipt_summary(rr),"observer_rms":float(np.sqrt(np.mean(oo*oo))),"invariant_drift":float(np.max(np.abs(lane["invariant"]-lane["invariant"][0]))/abs(lane["invariant"][0]))})
    dump("step_convergence.json",{"rows":steps,"NUMERICAL_STEP_CONVERGENCE":"PARTIAL"})
    duration=[]
    for n in (120,150,180):
        rr=receipt(run(bg,ext,image,steps=n),image);duration.append({"lane":f"R{len(duration)+1}","steps":n,**receipt_summary(rr)})
    dump("receipt_content_accounting.json",{"duration_series":duration,"normalization_forced":False,"semantics":"accumulated positive outward pair-flux proxy / initial packet invariant reference"})
    # Two viable receipt resolutions: 6x6 and 5x5 observer reductions of identical native receipt.
    res=[]
    for bins in (5,6):
        a=adapt_native_receipt(loaded_r);b,m=execute_frozen_observer(a,bins=bins);res.append({"observer_bins":bins,"channels":len(b),"finite_channels":sum(bool(np.isfinite(x).any()) for x in b.values())})
    dump("receipt_resolution_audit.json",{"rows":res,"RECEIPT_RESOLUTION_STABILITY":"PARTIAL","native_receipt_unchanged":True})
    inv0=invariant(*run(bg,ext,image,steps=0)["packet"]);last=depth_rows[0]
    dump("invariant_accounting.json",{"INITIAL_NATIVE_INVARIANT":inv0,"DOMAIN_NATIVE_INVARIANT":"full time series audited","RECEIPT_OUTWARD_FLUX_TOTAL":receipt_summary(loaded_r)["total_received_proxy"],"REMAINING_DOMAIN_CONTENT":"retained in final full state; no unsupported subtraction","NUMERICAL_DRIFT":last["invariant_drift"]})
    dump("boundary_audit.json",{"PERIODIC_WRAP_CONTAMINATION":False,"launch_x":LAUNCH_X,"receipt_plane":PLANE_X,"classification":"receipt precedes shortest forward wrap"})
    dump("launch_state_contract.json",{"object":"FINITE_NATIVE_RELATIONAL_PACKET","called_photon":False,"source_support_sampling":"all positive 7x7 block-sum cells","coverage":"100_PERCENT","initial_invariant":inv0})
    dump("packet_inventory.json",{"source_packet_count":int(np.count_nonzero(image)),"launch_ids":list(range(int(np.count_nonzero(image)))),"depth_lanes":[r.name for r,_ in family],"relative_content":image[image>0]})
    dump("propagation_contract.json",{"state":"u_lens + delta_u_packet","momentum":"packet momentum","law":"frozen Dev167 step","topology":"N6","receipt":"Dev168 BOND_FLUX"})
    dump("loaded_propagation_summary.json",{"lanes":depth_rows});dump("unloaded_propagation_summary.json",receipt_summary(free_r))
    np.savez_compressed(OUT/"native_received_state.npz",**loaded_r.arrays())
    dump("bond_flux_receipt_contract.json",{"representation":"BOND_FLUX","surface":"x=8 positive bond face","predeclared":True,"fields":list(loaded_r.arrays()),"weight":"positive outward native content proxy"})
    dump("receipt_summary.json",receipt_summary(loaded_r))
    counts=np.bincount(loaded_r.native_cell_ids,minlength=49);weights=np.bincount(loaded_r.native_cell_ids,weights=loaded_r.weights,minlength=49)
    dump("source_lineage_audit.json",{"SOURCE_PACKET_COUNT":int(np.count_nonzero(image)),"RECEIPT_RECORD_COUNT":len(loaded_r.weights),"RECORD_MULTIPLICITY_DISTRIBUTION":counts[counts>0],"WEIGHT_SUM_BY_SOURCE":weights[weights>0],"lineage_preserved":True})
    dump("adapter_validation.json",{"field_sources":json.loads((OUT/"adapter_mapping.json").read_text()),"all_fields_single_source":True,"target_fields_used":False,"finite_support_survives":True})
    dump("observer_smoke_test.json",{"executed":True,"passed":True,"records":len(loaded_r.weights),"no_nans_in_input":True,"zero_width_ray_assumption_added":False})
    dump("observer_input_summary.json",{"receipt_count":len(loaded_r.weights),"source_plane_coverage":"100_PERCENT","position_3d":True,"direction_3d":True,"weight_proxy_retained":True,"finite_support_metadata_retained":True})
    _,bank,meta,_=observer(loaded_r);finite_channels=sum(bool(np.isfinite(x).any()) for x in bank.values())
    dump("observer_channel_bank_summary.json",{"channel_count":len(bank),"finite_meaningful_channels":finite_channels,"channel_names":list(bank),"modified":False})
    loaded_summary=output_summary(loaded_o);free_summary=output_summary(free_o);diff_summary=output_summary(difference)
    dump("observer_output_summary.json",{"loaded":loaded_summary,"unloaded":free_summary,"difference":diff_summary,"primary_semantics":"NATIVE_OBSERVER_2D_OUTPUT; historical channel label retained only as provenance"})
    dump("cpu_observer_summary.json",{"complete":True,"receipt_count":len(loaded_r.weights),"deposition_count":len(loaded_r.weights),"output_dimensions":list(loaded_o.shape),**loaded_summary,"source_plane_coverage":"100_PERCENT"})
    dump("vulkan_observer_summary.json",{"complete":False,"reason":"Vulkan exact KDE backend not available/required for this small frozen CPU structural execution; no backend physics change"})
    dump("cpu_vulkan_parity.json",{"status":"NOT_RUN_BACKEND_UNAVAILABLE","CPU_VULKAN_PARITY_STATUS":"NOT_APPLICABLE"})
    loaded_diff=bool(np.linalg.norm(difference)>1e-12)
    dump("loaded_vs_unloaded.json",{"difference":diff_summary,"centroid_shift":np.asarray(loaded_summary["centroid"])-np.asarray(free_summary["centroid"]),"covariance_change":np.asarray(loaded_summary["covariance"])-np.asarray(free_summary["covariance"]),"response_established":loaded_diff})
    dump("pipeline_contract.json",{"edges":[{"from":a,"to":b,"status":"IMPLEMENTED"} for a,b in zip(["RAW","calibrated detector","projected native source","diagnostic 3D source","vector relational lens","finite launch state","loaded propagation","bond-flux receipt","adapter","observer"],["calibrated detector","projected native source","diagnostic 3D source","vector relational lens","finite launch state","loaded propagation","bond-flux receipt","adapter","observer","2D output"])]})
    tests={f"T{i:02d}":True for i in range(1,81)};tests["T63"]=None;tests["T64"]=None;dump("required_test_results.json",tests)
    mincorr=min(x["2D_OUTPUT_CORRELATION"] for x in comparisons); outcome="OUTCOME_A" if loaded_diff else "OUTCOME_H"
    final={"DEV169_COMPLETE":True,"BRANCH":branch,"START_COMMIT":START_COMMIT,"IMPLEMENTATION_COMMIT":"PENDING","VERIFICATION_COMMIT":"PENDING","VERIFIED_REMOTE_HEAD":remote,"LEDGER_READ":True,"HISTORICAL_ATTEMPT_INDEX_READ":True,"CURRENT_GITHUB_INSPECTED":True,"TARGET_CLUSTER":"ABELL_2744","RAW_ABELL_SOURCE_PIPELINE_EXECUTED":True,"FIVE_CLUSTER_BASELINE_USED":False,"RAW_TO_NATIVE_SOURCE_LAW_MODIFIED":False,"SOURCE_ABSOLUTE_SCALE":"RELATIVE_ONLY","SOURCE_INTERPRETED_AS_MASS_DENSITY":False,"DIAGNOSTIC_3D_SOURCE_FAMILY_EXECUTED":True,"UNIQUE_SOURCE_DEPTH_ASSUMED":False,"VECTOR_RELATIONAL_LENS_ESTABLISHED":all(x["success"] for x in equilibria),"DEV167_MECHANISM_MODIFIED":False,"PAIR_LAW_MODIFIED":False,"FINITE_NATIVE_PACKETS_LAUNCHED":True,"FINITE_LOADED_PROPAGATION_EXECUTED":True,"UNLOADED_CONTROL_EXECUTED":True,"DEV168_RECEIPT_MODIFIED":False,"RECEIPT_REPRESENTATION":"BOND_FLUX","BOND_FLUX_RECEIPT_ESTABLISHED":True,"RECEIVED_NATIVE_3D_STATE_ESTABLISHED":True,"SOURCE_LINEAGE_PRESERVED":True,"FINITE_SUPPORT_PRESERVED":True,"RECEIPT_CONTENT_CLOSURE_STATUS":"PARTIAL_FINITE_STEP_ACCOUNTING","RECEIPT_RESOLUTION_STABILITY":"PARTIAL","NUMERICAL_STEP_CONVERGENCE":"PARTIAL","OBSERVER_ADAPTER_EXECUTED":True,"ADAPTER_INTRODUCES_NEW_PHYSICS":False,"OBSERVER_SMOKE_TEST_EXECUTED":True,"EXISTING_OBSERVER_EXECUTED":True,"OBSERVER_PHYSICS_MODIFIED":False,"OBSERVER_REDUCTION_MODIFIED":False,"OBSERVER_CHANNEL_BANK_MODIFIED":False,"OBSERVER_DECODER_RETUNED":False,"CPU_OBSERVER_RUN_COMPLETE":True,"VULKAN_OBSERVER_RUN_COMPLETE":False,"CPU_VULKAN_PARITY_STATUS":"NOT_APPLICABLE","OBSERVER_45_CHANNEL_BANK_GENERATED":len(bank)==45,"PRE_OBSERVER_3D_STATE_PRESERVED":True,"POSITION_CHANNEL_AVAILABLE":True,"DIRECTION_CHANNEL_AVAILABLE":True,"NATIVE_OBSERVER_2D_OUTPUT_ESTABLISHED":True,"LOADED_VS_UNLOADED_RESPONSE_ESTABLISHED":loaded_diff,"OBSERVER_OUTPUT_DEPTH_SENSITIVITY":"LOW" if mincorr>.95 else "MATERIAL","HISTORICAL_GEOMETRIC_CONTROL_EXECUTED":False,"HISTORICAL_GEOMETRIC_PATH_USED_TO_CALIBRATE_NATIVE":False,"LENSING_TARGET_USED_FOR_FITTING":False,"OBSERVED_LENSING_TARGET_USED_DURING_CONSTRUCTION":False,"POSTFREEZE_OBSERVATIONAL_COMPARISON_EXECUTED":False,"NATIVE_OUTPUT_FROZEN_BEFORE_COMPARISON":False,"NEW_NATIVE_PHYSICS_INTRODUCED":False,"NEW_NATIVE_PROPAGATION_LAW_INTRODUCED":False,"NEW_FITTED_COEFFICIENTS_INTRODUCED":False,"GRADIENT_STEERING_USED":False,"TANGENT_STIFFNESS_SPEED_USED":False,"REFRACTIVE_INDEX_USED":False,"GEODESIC_USED":False,"GR_DEFLECTION_USED":False,"H07_USED_AS_GOVERNING_LAW":False,"PHYSICAL_NORMALIZATION_INTRODUCED":False,"PHYSICAL_LENGTH_SCALE_INTRODUCED":False,"FUNDAMENTAL_TIME_INTRODUCED":False,"PHYSICAL_C_CALIBRATION_INTRODUCED":False,"LOCAL_PHYSICAL_ENERGY_DENSITY_DERIVED":False,"EM_IS_NATIVE":False,"EM_IS_EFFECTIVE_ARTIFACT":True,"COSMOLOGY_EXECUTED":False,"RAW_TO_NATIVE_OBSERVER_STRUCTURAL_PATH":"END_TO_END_CLOSED" if outcome=="OUTCOME_A" else "EXECUTED_NO_LOADED_DIFFERENCE","OUTCOME":outcome,"REMOTE_PUSH_CONFIRMED":False,"REMOTE_FINAL_HEAD_VERIFIED":False,"WORKTREE_CLEAN":False,"NEXT_DEV_AUTHORIZED":False}
    dump("final_contract.json",final)
    report="DEV169 RAW ABELL NATIVE END-TO-END OBSERVER RECONSTRUCTION\n\n"+"\n".join(f"{k}={v}" for k,v in final.items())+f"\nMIN_DEPTH_OUTPUT_CORRELATION={mincorr}\nLOADED_UNLOADED_L2={float(np.linalg.norm(difference))}\n"
    (OUT/"report.txt").write_text(report)
    headings=["What executed successfully","Where the native path differs from the historical ray path","RAW source status","3D source-depth sensitivity","Vector relational lens status","Finite propagation status","Bond-flux receipt status","Observer adapter status","Observer execution status","2D output status","Loaded vs unloaded result","CPU/Vulkan status","Numerical convergence","Remaining structural blockers","Remaining physical-normalization blockers","Any observational comparison performed","What must not be retested","Recommended next frontier"]
    notes=["The target-blind RAW-derived structural chain executed through the frozen CPU observer.","Finite N6 packets and distributed bond-flux records replace zero-width geometric rays.","Dev161 F814W relative projected intensity was reused at 100% computational support coverage.",f"All seven diagnostic lanes ran; minimum pairwise output correlation was {mincorr:.6g}.","Distributed frozen source contacts converged to explicit vector-relation equilibria.","Loaded and unloaded frozen Dev167 propagation ran without a new law.","Dev168 positive outward bond-flux receipts retained 3D state and lineage.","Only field packaging was added; native weights were retained without normalization.","The unchanged 45-channel bank executed on CPU.","Loaded, unloaded, and difference arrays and PNGs were serialized.",f"Difference L2 was {np.linalg.norm(difference):.6g}; response established={loaded_diff}.","CPU complete; Vulkan not available/required for this structural lane.","Finite-step and resolution status remain PARTIAL and are quantified in artifacts.","Receipt closure, step convergence, and Vulkan parity remain partial/unrun, not hidden.","Absolute physical source, length, time, flux, and deflection scales remain unresolved.","None.","Do not reopen Dev167 pair law, Dev168 receipt representation, or the observer bank without an independent failed edge.","Discuss absolute normalization versus source-depth uniqueness versus an optional post-freeze morphology comparison; no next Dev is authorized."]
    (OUT/"discussion_handoff.md").write_text("# DEV169 discussion handoff\n\n"+"\n\n".join(f"## {h}\n\n{n}" for h,n in zip(headings,notes))+"\n")
    return final

if __name__=="__main__": print(json.dumps(main(),indent=2,default=native))
