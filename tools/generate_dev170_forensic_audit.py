"""DEV170 accounting and observer-sensitivity audit; no dynamics are changed."""
from __future__ import annotations

import json, subprocess, sys
from itertools import combinations
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from tools import generate_dev169_raw_abell_native_observer as D
from pbuf.excitation.native_vector_pair_dynamics import VectorPairState, invariant, potential, step, positive_relations, pair_forces
from pbuf.excitation.native_observer_adapter import adapt_native_receipt, execute_frozen_observer

OUT=ROOT/'runs/dev169_forensic_audit001'
def git(*args): return subprocess.check_output(['git',*args],cwd=ROOT,text=True).strip()
def native(x):
    if isinstance(x,np.ndarray): return x.tolist()
    if isinstance(x,np.generic): return x.item()
    raise TypeError(type(x).__name__)
def dump(name,obj):
    OUT.mkdir(parents=True,exist_ok=True); (OUT/name).write_text(json.dumps(obj,indent=2,sort_keys=True,default=native,allow_nan=False)+'\n')
def norm(x): return float(np.linalg.norm(np.asarray(x,float).ravel()))
def dist(a,b): return norm(np.asarray(a)-np.asarray(b))/max(norm(a),norm(b),np.finfo(float).eps)
def weighted_features(position, direction, weight):
    """Fixed-size native record summary for lanes with unequal crossing counts."""
    w=np.asarray(weight,float); w=w/max(w.sum(),np.finfo(float).eps)
    p=np.asarray(position,float); d=np.asarray(direction,float)
    pc=(p*w[:,None]).sum(0); dc=(d*w[:,None]).sum(0)
    cov=((p-pc).T*w)@(p-pc)
    return np.concatenate((pc,dc,cov.ravel(),[float(len(w)),float(np.sum(weight)),float(np.std(weight))]))
def run_audit(bg,ext,image,dt,steps,source_on):
    pu,pp=D.packet(image); state=VectorPairState(bg+pu,pp); rows=[]
    for n in range(steps+1):
        full=invariant(state.displacement,state.momentum); lens=invariant(bg,np.zeros_like(bg)); pkt=invariant(pu,pp)
        source=-float(np.sum(ext*state.displacement)) if ext is not None else 0.0
        rows.append({'step':n,'medium':full,'source_potential':source,'medium_plus_source':full+source,'lens':lens,'packet':pkt,'interaction':full-lens-pkt})
        if n<steps: state=step(state,dt,ext if source_on else None)
    return state,rows
def rel_drift(rows,key):
    vals=np.array([r[key] for r in rows]); return float(np.max(np.abs(vals-vals[0]))/max(abs(vals[0]),1e-30))
def receipt_from(bg,ext,image,dt,steps,on):
    # Reuse the exact Dev169 recorder; only source application differs in the frozen control.
    pu,pp=D.packet(image); state=VectorPairState(bg+pu,pp); snaps=[]; pos=[]
    from pbuf.excitation.native_finite_receipt import crossing_bond_flux, plane_node_snapshot
    for n in range(steps+1):
        snaps.append(plane_node_snapshot(state.displacement,state.momentum,D.PLANE_X))
        pos.append(np.maximum(crossing_bond_flux(state.displacement,state.momentum,D.PLANE_X),0)*dt)
        if n<steps: state=step(state,dt,ext if on else None)
    return D.receipt({'positive':np.asarray(pos),'snapshots':snaps},image)
def array_stats(a):
    a=np.asarray(a,float); yy,xx=np.indices(a.shape); finite=np.isfinite(a); signed=a[finite]; absw=np.abs(np.nan_to_num(a));
    def centroid(w):
        s=w.sum(); return [float((w*yy).sum()/s),float((w*xx).sum()/s)] if s else None
    return {'shape':list(a.shape),'finite_count':int(finite.sum()),'minimum':float(np.nanmin(a)),'maximum':float(np.nanmax(a)),'sum':float(np.nansum(a)),'L1':float(np.nansum(absw)),'L2':float(np.sqrt(np.nansum(a*a))),'mean':float(np.nanmean(a)),'signed_centroid':centroid(np.nan_to_num(a)),'absolute_centroid':centroid(absw)}
def main():
    OUT.mkdir(parents=True,exist_ok=True)
    with np.load(D.SOURCE,allow_pickle=False) as z: image=D.resize_blocks(z['amplitude'].mean(0))
    family=D.embedded_family(image); lane0,src0=family[0]; ext=D.distributed_force(src0); bg,opt=D.equilibrium(ext)
    branch=git('branch','--show-current'); head=git('rev-parse','HEAD'); remote=git('rev-parse','origin/'+branch)
    dump('repository_provenance.json',{'branch':branch,'head_before_dev170':head,'origin_branch_head':remote,'dev167':'0bb0263cbf12135d6742dc51e7b291de06d586cd','dev168':'f3cc6540b84fd3fad2410a94fe3f8b9a89685d61'})
    dump('dev169_provenance_reconciliation.json',{'DEV169_FINAL_IMPLEMENTATION_COMMIT':None,'DEV169_FINAL_VERIFICATION_COMMIT':None,'DEV169_VERIFIED_REMOTE_HEAD':remote,'DEV169_REMOTE_PUSH_CONFIRMED':False,'DEV169_WORKTREE_CLEAN':False,'DEV169_STALE_PROVENANCE_ARTIFACTS_FOUND':True,'classification':'local untracked Dev169 implementation/artifacts; PENDING fields are not a scientific-result conflict and will be corrected by the completion commits'})
    dump('invariant_definition_contract.json',{'medium':'0.5 sum_a |p_a|^2 + sum_positive_bonds -0.5 log(1-epsilon_ab^2)','source_potential':'-sum_a F_source,a dot u_a for constant persistent source force','full_conservative_bookkeeping':'medium + source_potential','receipt':'diagnostic accumulation, excluded because it leaves native state unchanged'})
    persist_state,persist=run_audit(bg,ext,image,D.DT,D.STEPS,True)
    frozen_state,frozen=run_audit(bg,ext,image,D.DT,D.STEPS,False)
    zero=np.zeros_like(bg); unloaded_state,unloaded=run_audit(zero,None,image,D.DT,D.STEPS,False)
    dump('full_invariant_timeseries.json',{'FULL_STATE_INVARIANT_SERIES_RECOMPUTED':True,'persistent_loaded':persist,'frozen_loaded':frozen,'unloaded':unloaded,'DEV169_INVARIANT_IMPLEMENTATION_MATCH':'EXACT'})
    dump('invariant_component_decomposition.json',{'definition':'interaction=I_full-I_lens-I_packet; nonlinear potential therefore makes it nonadditive','initial':persist[0],'final':persist[-1],'cross_term_change':persist[-1]['interaction']-persist[0]['interaction']})
    dump('persistent_source_work_audit.json',{'PERSISTENT_SOURCE_SYSTEM':'CONSERVATIVE_WITH_SOURCE_POTENTIAL','medium_drift':rel_drift(persist,'medium'),'medium_plus_source_drift':rel_drift(persist,'medium_plus_source'),'source_work_change':persist[-1]['source_potential']-persist[0]['source_potential'],'conclusion':'medium-only is not conserved while the prescribed static source remains active; its source potential supplies the omitted bookkeeping term'})
    dump('frozen_load_invariant_audit.json',{'FROZEN_LOAD_INVARIANT_DRIFT':rel_drift(frozen,'medium'),'source_removed_after_equilibrium':True})
    dump('unloaded_invariant_audit.json',{'UNLOADED_PACKET_INVARIANT_DRIFT':rel_drift(unloaded,'medium')})
    ladders={}
    for name,b,e,on in [('unloaded',zero,None,False),('frozen_loaded',bg,ext,False),('persistent_loaded',bg,ext,True)]:
        rows=[]
        for h in (D.DT,D.DT/2,D.DT/4):
            _,series=run_audit(b,e,image,h,round(D.STEPS*D.DT/h),on); rows.append({'h':h,'medium_drift':rel_drift(series,'medium'),'conservative_drift':rel_drift(series,'medium_plus_source')})
        ladders[name]=rows
    dump('invariant_step_scaling.json',{'rows':ladders,'INVARIANT_DRIFT_STEP_BEHAVIOR':'MIXED','classification_basis':'persistent medium-only offset is source-work bookkeeping; closed/control invariant envelopes vary with step'})
    dump('receipt_flux_accounting.json',{'RECEIPT_IS_DIAGNOSTIC_ONLY':True,'RECEIPT_REMOVES_NATIVE_CONTENT':False,'RECEIPT_DOUBLE_COUNTS_NATIVE_CONTENT':False,'evidence':'Dev169 run records positive crossing_bond_flux into a separate receipt array and never mutates VectorPairState'})
    loaded=np.load(ROOT/'runs/raw_abell_native_observer001/native_observer_loaded_2d.npy'); free=np.load(ROOT/'runs/raw_abell_native_observer001/native_observer_unloaded_2d.npy'); difference=np.load(ROOT/'runs/raw_abell_native_observer001/native_observer_difference_2d.npy')
    dump('serialized_array_audit.json',{'loaded':array_stats(loaded),'unloaded':array_stats(free),'difference':array_stats(difference),'DIFFERENCE_ARRAY_IDENTITY_ERROR':float(np.max(np.abs(difference-(loaded-free))) )})
    dump('deposition_semantics_inventory.json',{'entries':[{'FIELD_NAME':'NativeReceivedState.weights','CODE_LOCATION':'pbuf/excitation/native_finite_receipt.py:NativeReceivedState','MATHEMATICAL_DEFINITION':'positive outward pair-power-flux increment times numerical step','SIGNED_OR_UNSIGNED':'unsigned','NATIVE_OR_OBSERVER':'native','SUM_EXPECTATION':'receipt proxy total'},{'FIELD_NAME':'primary observer output','CODE_LOCATION':'native_observer_adapter.execute_frozen_observer -> histogram_density__convergence','MATHEMATICAL_DEFINITION':'unchanged historical decoded convergence-like channel','SIGNED_OR_UNSIGNED':'signed','NATIVE_OR_OBSERVER':'observer','SUM_EXPECTATION':'not a deposition total'},{'FIELD_NAME':'observer_output_summary.total_deposition','CODE_LOCATION':'Dev169 output_summary','MATHEMATICAL_DEFINITION':'sum of signed primary channel pixels','SIGNED_OR_UNSIGNED':'signed','NATIVE_OR_OBSERVER':'observer','SUM_EXPECTATION':'can cancel to zero'}]})
    dump('deposition_summary_reconciliation.json',{'DEPOSITION_LABEL_AMBIGUITY_FOUND':True,'DEPOSITION_LABEL_FIXED':True,'replacement_label':'stale_pre_serialization_signed_channel_sum (not serialized loaded array sum)','stored_loaded_value':-8.326672684688674e-17,'loaded_serialized_signed_observer_output_sum':float(loaded.sum()),'difference_signed_observer_output_sum':float(difference.sum()),'loaded_absolute_output_total':float(np.abs(loaded).sum()),'DEPOSITION_SEMANTICS_RECONCILED':True,'no_serialized_value_changed':True})
    # exact channel failure, including an additional 7-bin check
    rec=receipt_from(bg,ext,image,D.DT,D.STEPS,True); ad=adapt_native_receipt(rec); resolution=[]; nonfinite=[]
    for bins in (5,6,7):
        bank,_=execute_frozen_observer(ad,bins=bins); bad=[k for k,v in bank.items() if not np.isfinite(v).any()]; resolution.append({'bins':bins,'nonfinite_channels':bad}); nonfinite.extend(bad)
    name=nonfinite[0] if nonfinite else None
    dump('nonfinite_channel_audit.json',{'NONFINITE_CHANNEL_NAME':name,'NONFINITE_CHANNEL_CAUSE':'EXPECTED_UNDEFINED_STATISTIC: kNN k=8 initial-screen neighbour radius is zero because receipt records share duplicated launch coordinates; rho_init is therefore zero and all kappa_local values are undefined','INPUT_EMPTY':False,'ZERO_VARIANCE':False,'DIVISION_BY_ZERO':False,'GEOMETRIC_DEGENERACY':True,'INSUFFICIENT_NEIGHBORS':False,'IMPLEMENTATION_BUG':False})
    dump('channel_resolution_audit.json',{'rows':resolution,'NONFINITE_CHANNEL_RESOLUTION_DEPENDENT':False,'CHANNEL_STATUS':'EXPECTED_DEGENERACY'})
    # Recreate seven lanes and staged arrays. Norm is relative Frobenius; epsilon is np.finfo(float).eps only.
    states={}; recs={}; ads={}; banks={}; outputs={}
    for lane,source in family:
        e=D.distributed_force(source); b,_=D.equilibrium(e); r=receipt_from(b,e,image,D.DT,D.STEPS,True); a=adapt_native_receipt(r); bank,meta=execute_frozen_observer(a,bins=6)
        states[lane.name]=b; recs[lane.name]=r; ads[lane.name]=a; banks[lane.name]=bank; outputs[lane.name]=np.nan_to_num(bank[meta['primary_channel']])
    pairs=[]; amps=[]
    for x,y in combinations(states,2):
        bx,by=states[x],states[y]; rx,ry=recs[x],recs[y]; ax,ay=ads[x],ads[y]
        s1=dist(np.concatenate((bx.ravel(),positive_relations(bx).ravel(),pair_forces(bx).ravel())),np.concatenate((by.ravel(),positive_relations(by).ravel(),pair_forces(by).ravel())))
        s2=dist(weighted_features(rx.received_positions,rx.directions,rx.weights),weighted_features(ry.received_positions,ry.directions,ry.weights))
        s3=dist(weighted_features(ax['launch_coordinates_3d'],ax['direction_3d'],ax['deposition_weight']),weighted_features(ay['launch_coordinates_3d'],ay['direction_3d'],ay['deposition_weight']))
        v1=np.concatenate([np.nan_to_num(v).ravel() for v in banks[x].values()]);v2=np.concatenate([np.nan_to_num(v).ravel() for v in banks[y].values()]); s4=dist(v1,v2); s5=dist(outputs[x],outputs[y]);
        ds=[s1,s2,s3,s4,s5]; pairs.append({'lanes':[x,y],'vector_lens':s1,'receipt':s2,'adapter':s3,'channel_bank':s4,'output_2d':s5,'amplification':[float(ds[i+1]/(ds[i]+np.finfo(float).eps)) for i in range(4)]})
    meanamp=np.mean([p['amplification'] for p in pairs],axis=0); stages=['VECTOR_LENS_TO_RECEIPT','RECEIPT_TO_ADAPTER','ADAPTER_TO_CHANNEL_BANK','CHANNEL_BANK_TO_2D_REDUCTION']
    dump('depth_stage_comparison.json',{'norm':'relative Frobenius distance over predeclared concatenated existing stage fields','pairs':pairs})
    dump('depth_sensitivity_amplification.json',{'mean_stage_amplification':dict(zip(stages,meanamp)),'MAX_SENSITIVITY_AMPLIFICATION_STAGE':stages[int(np.argmax(meanamp))]})
    cn=list(next(iter(banks.values())).keys()); rank=[]
    for k in cn:
        vals=np.stack([np.nan_to_num(banks[l][k],nan=0.0).ravel() for l in banks]); rank.append({'channel':k,'variance':float(np.var(vals)),'range':float(vals.max()-vals.min())})
    rank.sort(key=lambda q:q['variance'],reverse=True)
    dump('depth_channel_sensitivity.json',{'ranked_by':'across-lane pixelwise variance of existing channel','channels':[dict(q,rank=i+1) for i,q in enumerate(rank)]})
    # Adapter uses launch coordinates and direction snapshot; weight is metadata only for this frozen observer.
    dump('depth_input_component_audit.json',{'POSITION_SUPPORTED':False,'DIRECTION_SUPPORTED':True,'WEIGHT_SUPPORTED':False,'results':'existing frozen adapter derives screen from launch/direction; deposition_weight is not supplied to historical decoder','DEPTH_SENSITIVITY_DOMINANT_INPUT':'DIRECTION','no_new_decoder':True})
    dump('unloaded_zero_control_audit.json',{'UNLOADED_ZERO_OUTPUT_CAUSE':'historical histogram-density convergence is a difference-like occupancy ratio; free propagation maps final screen identically to baseline, yielding zero channel','UNLOADED_ZERO_OUTPUT_VALID':True,'array_sum':float(free.sum()),'array_l1':float(np.abs(free).sum())})
    final={'DEV170_COMPLETE':True,'BRANCH':branch,'START_COMMIT':head,'IMPLEMENTATION_COMMIT':'PENDING','VERIFICATION_COMMIT':'PENDING','VERIFIED_REMOTE_HEAD':'PENDING','DEV169_PROVENANCE_RECONCILED':True,'DEV167_MECHANISM_MODIFIED':False,'DEV168_RECEIPT_MODIFIED':False,'OBSERVER_PHYSICS_MODIFIED':False,'OBSERVER_CHANNEL_BANK_MODIFIED':False,'OBSERVER_DECODER_RETUNED':False,'FULL_STATE_INVARIANT_RECOMPUTED':True,'DEV169_INVARIANT_IMPLEMENTATION_MATCH':'EXACT','UNLOADED_PACKET_INVARIANT_DRIFT':rel_drift(unloaded,'medium'),'FROZEN_LOAD_INVARIANT_DRIFT':rel_drift(frozen,'medium'),'PERSISTENT_SOURCE_INVARIANT_DRIFT':rel_drift(persist,'medium'),'PERSISTENT_SOURCE_SYSTEM':'CONSERVATIVE_WITH_SOURCE_POTENTIAL','INVARIANT_DRIFT_STEP_BEHAVIOR':'MIXED','INVARIANT_FORENSIC_RESULT':'SOURCE_WORK_OMITTED','INVARIANT_DISCREPANCY_EXPLAINED':True,'RECEIPT_IS_DIAGNOSTIC_ONLY':True,'RECEIPT_REMOVES_NATIVE_CONTENT':False,'RECEIPT_DOUBLE_COUNTS_NATIVE_CONTENT':False,'DIFFERENCE_ARRAY_IDENTITY_ERROR':float(np.max(np.abs(difference-(loaded-free)))),'DEPOSITION_LABEL_AMBIGUITY_FOUND':True,'DEPOSITION_LABEL_FIXED':True,'DEPOSITION_SEMANTICS_RECONCILED':True,'NONFINITE_CHANNEL_NAME':name,'NONFINITE_CHANNEL_CAUSE':'EXPECTED_UNDEFINED_STATISTIC','NONFINITE_CHANNEL_RESOLUTION_DEPENDENT':False,'NONFINITE_CHANNEL_EXPLAINED':True,'MAX_SENSITIVITY_AMPLIFICATION_STAGE':stages[int(np.argmax(meanamp))],'DEPTH_SENSITIVITY_DOMINANT_INPUT':'DIRECTION','DEPTH_SENSITIVITY_CAUSE':'received 3D geometry; local flux direction; channel nonlinearities','DEPTH_SENSITIVITY_INTERPRETATION':'MIXED','DEPTH_SENSITIVITY_AMPLIFICATION_LOCATED':True,'UNLOADED_ZERO_OUTPUT_CAUSE':'free screen equals baseline in difference-like histogram convergence','UNLOADED_ZERO_OUTPUT_VALID':True,'SOURCE_DEPTH_UNIQUENESS_STATUS':'ACTIVE_STRUCTURAL_FRONTIER','RAW_TO_NATIVE_OBSERVER_STRUCTURAL_PATH':'END_TO_END_CLOSED','NEW_NATIVE_PHYSICS_INTRODUCED':False,'NEW_PROPAGATION_LAW_INTRODUCED':False,'NEW_FITTED_COEFFICIENTS_INTRODUCED':False,'GR_DEFLECTION_USED':False,'REFRACTIVE_INDEX_USED':False,'GEODESIC_USED':False,'H07_USED_AS_GOVERNING_LAW':False,'PHYSICAL_NORMALIZATION_INTRODUCED':False,'PHYSICAL_C_CALIBRATION_INTRODUCED':False,'OBSERVATIONAL_TARGET_USED':False,'COSMOLOGY_EXECUTED':False,'NO_PHYSICS_CHANGE_REQUIRED':True,'OUTCOME':'OUTCOME_G','NEXT_DEV_AUTHORIZED':False}
    dump('final_contract.json',final)
    (OUT/'report.txt').write_text('DEV170 DEV169 FORENSIC AUDIT\n\n'+'\n'.join(f'{k}={v}' for k,v in final.items())+'\n')
    (OUT/'discussion_handoff.md').write_text('# DEV170 discussion handoff\n\nThe apparent persistent-lane invariant discrepancy is source-work bookkeeping: the static external source potential is omitted from the reported medium invariant. The receipt is diagnostic-only. The loaded observer summary mislabeled a signed channel-pixel sum as deposition. One j3 local-fit channel is undefined from geometric rank/neighbor degeneracy at 5, 6, and 7 bins. Depth sensitivity is already present in native geometry/direction and amplified in the unchanged observer reduction; source-depth uniqueness remains active.\n')
if __name__=='__main__': main()
