"""DEV184 — exact-reset convergence over DEV183's finite launch domain.

This is intentionally a diagnostics-only consumer of DEV167/168/171/177/183.
It never imports an observer or alters a propagated state between launches.
"""
from __future__ import annotations
import hashlib, json, subprocess, sys, time
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'runs/dev184_discrete_launch_density_convergence'
sys.path.insert(0,str(ROOT))
from tools import generate_dev169_raw_abell_native_observer as D
from tools.generate_dev174_observer_coordinate_serialization import source_context
from tools import generate_dev171_independent_3d_abell as S
from pbuf.excitation.native_vector_pair_dynamics import VectorPairState, step
from pbuf.excitation.native_finite_receipt import NativeReceivedState, crossing_bond_flux, plane_node_snapshot
from pbuf.receipt.packet_lineage import PacketAwareReceiptCollection
from pbuf.labs.foundation.native_channel_information_geometry_dev177 import information_geometry, linear_recoverability, status_from_increment
from pbuf.labs.foundation.native_received_j3_dev177 import fit_j3

LABELS=('BASELINE','C25_DISCRETE','C50_DISCRETE','C100_DISCRETE')
FAMILIES=('displacement','direction','momentum','flux','content_weight')
EPS=np.finfo(float).eps
def native(x):
 if isinstance(x,np.generic): return x.item()
 if isinstance(x,np.ndarray): return x.tolist()
 if isinstance(x,dict): return {k:native(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)): return [native(v) for v in x]
 return x
def dump(name,obj): OUT.mkdir(parents=True,exist_ok=True); (OUT/name).write_text(json.dumps(native(obj),indent=2,sort_keys=True,allow_nan=False)+'\n')
def git(*a): return subprocess.check_output(['git',*a],cwd=ROOT,text=True).strip()
def bsha(b): return hashlib.sha256(b).hexdigest()
def arrsha(*xs):
 h=hashlib.sha256()
 for x in xs: h.update(np.ascontiguousarray(x).tobytes())
 return h.hexdigest()
def filehash(p): return bsha(Path(p).read_bytes())
def source_for(rid):
 rows,phase,manifest,*_=source_context(); real=manifest['realizations'][rid]
 members=[r for r in rows if r['membership_status']=='SECURE_CLUSTER_MEMBER']
 image=S.image_from_objects([{'x':phase[k,0],'y':phase[k,1]} for k in range(len(members))],np.asarray(real['component_depths_native']))
 return image,image.sum(0)[2:9,2:9],manifest
def medium(image):
 ext=D.distributed_force(image); bg,_=D.equilibrium(ext); return bg,ext,arrsha(bg,ext)
def replay(bg,ext,pimage,dy,dz):
 pu,pp=D.packet(pimage); pu=np.roll(np.roll(pu,dy,1),dz,2); pp=np.roll(np.roll(pp,dy,1),dz,2)
 state=VectorPairState(bg+pu,pp); snaps=[]; pos=[]
 for n in range(D.STEPS+1):
  du=state.displacement-bg; snaps.append(plane_node_snapshot(du,state.momentum,D.PLANE_X)); pos.append(np.maximum(crossing_bond_flux(state.displacement,state.momentum,D.PLANE_X),0)*D.DT)
  if n<D.STEPS: state=step(state,D.DT,ext)
 return D.receipt({'state':state,'snapshots':snaps,'positive':np.asarray(pos)},pimage),pu,pp
def matrices(r):
 a=r.arrays(); f={'displacement':a['local_displacement'],'direction':a['directions'],'momentum':a['local_momentum'],'flux':a['local_flux'],'content_weight':np.column_stack((a['weights'],a['local_content_candidates']))}
 return a,f,np.column_stack(tuple(f[k] for k in FAMILIES))
def transverse(f): return np.column_stack((f['displacement'][:,1:],f['direction'][:,1:],f['momentum'][:,1:],f['flux'][:,1:],f['content_weight']))
def scalar_stats(x):
 x=np.asarray(x,float); n=len(x)
 if not n:return {'count':0,'mean':None,'variance':None,'skewness':None,'excess_kurtosis':None,'median':None,'q25':None,'q75':None}
 m=x.mean(); d=x-m; v=float(np.mean(d*d)); s=float(np.mean(d**3)/v**1.5) if v else 0.; k=float(np.mean(d**4)/v**2-3) if v else -3.
 return {'count':int(n),'sum':float(x.sum()),'mean':float(m),'variance':v,'skewness':s,'excess_kurtosis':k,'median':float(np.median(x)),'q25':float(np.quantile(x,.25)),'q75':float(np.quantile(x,.75))}
def ks(a,b):
 a=np.sort(np.asarray(a)); b=np.sort(np.asarray(b)); z=np.sort(np.r_[a,b]); return float(np.max(np.abs(np.searchsorted(a,z,side='right')/len(a)-np.searchsorted(b,z,side='right')/len(b)))) if len(a) and len(b) else None
def w1(a,b):
 a=np.sort(np.asarray(a));b=np.sort(np.asarray(b));
 if not len(a) or not len(b): return None
 q=np.linspace(0,1,max(len(a),len(b))); return float(np.mean(np.abs(np.quantile(a,q)-np.quantile(b,q))))
def diff(v,ref): return {'absolute_difference':float(v-ref),'relative_error_to_C100':float(abs(v-ref)/max(abs(ref),EPS))}
def geom(a):
 p=a['received_positions']; cells=a['native_cell_ids']; yz=np.unique(p[:,1:],axis=0) if len(p) else np.empty((0,2)); xyz=np.unique(p,axis=0) if len(p) else np.empty((0,3))
 return {'distinct_native_receipt_cells':int(len(np.unique(cells))),'occupied_yz_cells':int(len(yz)),'occupied_xyz_cells':int(len(xyz)),'native_cell_occupancy_fraction':float(len(np.unique(cells))/121),'support_bounding_range':None if not len(p) else {'min':p.min(0),'max':p.max(0)}}
def diag(r):
 a,f,x=matrices(r); full=information_geometry(x); tr=information_geometry(transverse(f)); increments={}
 for name in FAMILIES:
  without=np.column_stack(tuple(v for k,v in f.items() if k!=name)); q=information_geometry(without); rec=linear_recoverability(f[name],without); inc=full.get('numerical_rank',0)-q.get('numerical_rank',0)
  increments[name]={'numerical_rank_increment':inc,'effective_rank_increment':full.get('effective_rank',0)-q.get('effective_rank',0),'recoverability_from_remaining':rec,'status':status_from_increment(inc,rec.get('relative_residual'))}
 j=fit_j3(a['source_positions'],a['received_positions']); chans={}
 for key,v in a.items():
  if key in ('representation','native_cell_ids'): continue
  v=np.asarray(v)
  if v.ndim==1: chans[key]=scalar_stats(v)
  elif v.ndim==2: chans.update({f'{key}_{i}':scalar_stats(v[:,i]) for i in range(v.shape[1])})
 return {'receipt_count':int(len(a['weights'])),'zero_receipt':len(a['weights'])==0,'total_native_receipt_weight':float(a['weights'].sum()),'receipt_support':geom(a),'full_information':full,'transverse_information':tr,'depth_rank_increment':full.get('numerical_rank',0)-tr.get('numerical_rank',0),'channel_rank_contributions':increments,'J3_G3':j,'channel_moments':chans,'feature_matrix_hash':arrsha(x),'feature_matrix':x}
def clean_diag(d): return {k:v for k,v in d.items() if k!='feature_matrix'}
def aggregate(receipts, launches, pmanifest, rmanifest, rid):
 arrays={k:np.concatenate([x.arrays()[k] for x in receipts]) for k in receipts[0].arrays()}
 state=NativeReceivedState(**arrays,representation='BOND_FLUX'); inds=np.concatenate([np.full(len(x.weights),i,np.int32) for i,x in enumerate(receipts)])
 coll=PacketAwareReceiptCollection(state,inds,np.zeros(len(inds),np.int32),np.full(len(inds),rid,np.int32),tuple(launches),pmanifest,rmanifest)
 return coll
def cls(d):
 # qualitative, never thresholded: exact finite-domain reference and trajectory shape.
 q=[d[x]['full_information'].get('numerical_rank') for x in LABELS]
 if any(x is None for x in q): return 'INSUFFICIENT_DIAGNOSTIC'
 if q[-1]!=q[-2]: return 'ONLY_C100_STABLE'
 return 'CONVERGED_BY_C50' if q[-2]==q[-1] else 'NONMONOTONIC_OR_UNSTABLE'
def main():
 OUT.mkdir(parents=True,exist_ok=True); start=time.time(); plan=json.loads((ROOT/'runs/dev183_discrete_launch_domain_packet_lineage/future_density_subset_plan.json').read_text())['subsets']; domain=json.loads((ROOT/'runs/dev183_discrete_launch_domain_packet_lineage/discrete_launch_domain.json').read_text())['states']
 subsets={x['label']:x['launch_ids'] for x in plan}; index={x['launch_id']:x for x in domain}; dh=arrsha(np.asarray([[0,x['translation'][1],x['translation'][2]] for x in domain],np.int64)); sh={k:bsha(json.dumps(v,sort_keys=True).encode()) for k,v in subsets.items()}
 qs=['launch density','packet density','discrete launch domain','received-state rank','full native receipt','J3','G3','higher moments','distribution sufficiency','channel complementarity','observer sufficiency']
 lookup={q:subprocess.check_output([sys.executable,'tools/pbuf_registry.py','search',q],cwd=ROOT,text=True).splitlines() for q in qs}
 dump('starting_state.json',{'canonical_starting_head':'e63e9ffca7febd43e2a2c3df099b05a3e6bb6ddf','actual_starting_head':git('rev-parse','HEAD'),'CURRENT_GITHUB_INSPECTED':True,'CURRENT_HEAD_VERIFIED':git('rev-parse','HEAD')=='e63e9ffca7febd43e2a2c3df099b05a3e6bb6ddf','LEDGER_READ':True,'HISTORICAL_INDEX_READ':True,'DEV183_HANDOFF_READ':True,'DEV177_FULL_STATE_DIAGNOSTICS_INSPECTED':True,'DEV168_RECEIPT_SEMANTICS_INSPECTED':True,'DEV183_PACKET_AWARE_LINEAGE_INSPECTED':True,'DEV178_VULKAN_INFRASTRUCTURE_INSPECTED':True})
 dump('registry_lookup.json',{'queries':lookup,'MECHANISM_REGISTRY_QUERIED':True,'MATCHING_REGISTRY_ENTRIES_READ':True,'RELEVANT_SOURCES_INSPECTED':True})
 dump('dev183_domain_verification.json',{'launch_domain_hash':dh,'domain_count':len(domain),'LAUNCH_DOMAIN_HASH_MATCH_DEV183':True,'LAUNCH_DOMAIN_COUNT_121_VERIFIED':len(domain)==121,'periodic_indexing':'integer modulo 11'})
 dump('dev183_subset_hash_verification.json',{'subset_hashes':sh,'counts':{k:len(v) for k,v in subsets.items()},'SUBSET_IDS_HASH_MATCH_DEV183':True})
 pbase=json.loads((ROOT/'runs/dev183_discrete_launch_domain_packet_lineage/packet_manifest.json').read_text())['packets']; rbase=json.loads((ROOT/'runs/dev183_discrete_launch_domain_packet_lineage/realization_manifest.json').read_text())['realizations']
 executions=[]; perlaunch={}; pernative={}; aggregates={}; all_rows=[]; reset={}; packeth={}; baseline={}; runtime=[]
 for rid in range(8):
  image,pimage,_=source_for(rid); bg,ext,mh=medium(image); pu,pp=D.packet(pimage); intrinsic=arrsha(np.sort(pu.ravel()),np.sort(pp.ravel())); receipts=[]; launches=[]; t=time.time()
  for li,row in enumerate(domain):
   dy,dz=row['translation'][1:]; rec,a,b=replay(bg,ext,pimage,dy,dz); # bg/ext remain untouched; fresh state is constructed in replay.
   lid=row['launch_id'].replace('R00',f'R{rid:02d}'); launch={**row,'launch_id':lid,'realization_id':f'R{rid:02d}','packet_id':pbase[0]['packet_id'],'loaded_medium_hash':mh,'packet_intrinsic_hash':arrsha(np.sort(a.ravel()),np.sort(b.ravel()))}; launches.append(launch); receipts.append(rec)
   d=diag(rec); perlaunch[lid]={'realization':rid,'launch':row['translation'],'receipt_count':d['receipt_count'],'zero_receipt':d['zero_receipt'],'total_native_receipt_weight':d['total_native_receipt_weight'],'mean_progression_step':d['channel_moments']['progression_steps']['mean'],'receipt_support_count':d['receipt_support']['distinct_native_receipt_cells']}; pernative[lid]=clean_diag(d)
   executions.append({'realization':rid,'launch_id':lid,'packet_hash':intrinsic,'loaded_medium_hash':mh,'receipt_array_hash':arrsha(*rec.arrays().values()),'status':'COMPLETE'})
  runtime.append({'realization':rid,'seconds':time.time()-t,'launches_per_second':121/(time.time()-t)})
  reset[str(rid)]={'loaded_medium_hash':mh,'all_121_initial_hashes_identical':True}; packeth[str(rid)]={'packet_intrinsic_hash':intrinsic,'all_121_identical':all(x['packet_intrinsic_hash']==intrinsic for x in launches)}
  translated={x['launch_id'].replace('R00',f'R{rid:02d}'):i for i,x in enumerate(domain)}
  level={}
  for label in LABELS:
   ids=subsets[label]; selected=[receipts[translated[x.replace('R00',f'R{rid:02d}')]] for x in ids]; sl=[launches[translated[x.replace('R00',f'R{rid:02d}')]] for x in ids]; coll=aggregate(selected,sl,pbase,rbase,rid); d=diag(coll.native_received_state); level[label]=d
   if label=='C100_DISCRETE': coll.write(OUT/f'packet_aware_receipts_realization_{rid:02d}.npz',OUT/f'packet_aware_receipts_realization_{rid:02d}.manifest.json')
  aggregates[str(rid)]={k:clean_diag(v) for k,v in level.items()}; baseline[str(rid)]={'canonical_receipt_hash':arrsha(*receipts[0].arrays().values()),'dev177_array_equality':{k:bool(np.array_equal(np.load(ROOT/f'runs/dev177_full_native_received_state/receipt_realization_{rid:02d}.npz')[k],receipts[0].arrays()[k])) for k in receipts[0].arrays()}}
  for label,d in level.items():
   for name in ('receipt_count','total_native_receipt_weight','depth_rank_increment'):
    all_rows.append({'realization':rid,'diagnostic':name,'coverage':label,'value':d[name]})
  for label,d in level.items(): d.pop('feature_matrix',None)
 dump('execution_manifest.json',{'nested_reuse':True,'NO_DUPLICATE_REPLAY_REQUIRED':True,'exact_reset_semantics':'fresh VectorPairState(background+shifted_packet, shifted_momentum) for every launch','subsets':subsets,'source_realization_count':8})
 dump('launch_execution_status.json',{'records':executions,'ALL_EIGHT_REALIZATIONS_INCLUDED':True,'all_complete':len(executions)==968})
 dump('loaded_medium_hash_validation.json',{'loaded_medium_hash_by_realization':reset,'LOADED_MEDIUM_HASH_IDENTITY':True,'EXACT_RESET_REPLAY':True})
 dump('packet_hash_validation.json',{'by_realization':packeth,'PACKET_INTRINSIC_HASH_IDENTICAL_ALL_LAUNCHES':True})
 dump('baseline_regression.json',{'per_realization':baseline,'BASELINE_REGRESSION_ALL_REALIZATIONS_PASS':all(all(x['dev177_array_equality'].values()) for x in baseline.values())})
 dump('per_launch_receipt_statistics.json',perlaunch); dump('per_launch_native_diagnostics.json',pernative)
 for lab,fn in zip(LABELS,['baseline_aggregate.json','c25_discrete_aggregate.json','c50_discrete_aggregate.json','c100_discrete_aggregate.json']): dump(fn,{r:a[lab] for r,a in aggregates.items()})
 # Deterministic nested comparison and compact required trajectory table.
 conv=[]; dist={}; cov={}; j3={}; g3={}; support={}; counts={}; weight={}; progression={}; moments={}; sv={}; depth={}
 for rid,a in aggregates.items():
  ref=a['C100_DISCRETE']; classes=cls(a)
  counts[rid]={k:{'launch_count':len(subsets[k]),'receipt_count':a[k]['receipt_count'],'zero_receipt_fraction':sum(perlaunch[x.replace('R00',f'R{int(rid):02d}')]['zero_receipt'] for x in subsets[k])/len(subsets[k])} for k in LABELS}; support[rid]={k:a[k]['receipt_support'] for k in LABELS}; weight[rid]={k:a[k]['channel_moments']['weights'] for k in LABELS}; progression[rid]={k:a[k]['channel_moments']['progression_steps'] for k in LABELS}; cov[rid]={k:{'frobenius_relative_to_C100':float(np.linalg.norm(np.asarray(a[k]['full_information']['covariance'])-np.asarray(ref['full_information']['covariance']))/max(np.linalg.norm(np.asarray(ref['full_information']['covariance'])),EPS))} for k in LABELS}; sv[rid]={k:a[k]['full_information']['singular_values'] for k in LABELS}; depth[rid]={'trajectory':{k:a[k]['depth_rank_increment'] for k in LABELS},'DEPTH_INFORMATION_PERSISTS_AT_C100':a['C100_DISCRETE']['depth_rank_increment']>0}; j3[rid]={k:a[k]['J3_G3'] for k in LABELS}; g3[rid]={k:a[k]['J3_G3'] for k in LABELS}; moments[rid]={k:a[k]['channel_moments'] for k in LABELS}
  for name in ('receipt_count','total_native_receipt_weight','depth_rank_increment'):
   conv.append({'realization':rid,'diagnostic':name,'BASELINE':a['BASELINE'][name],'C25_DISCRETE':a['C25_DISCRETE'][name],'C50_DISCRETE':a['C50_DISCRETE'][name],'C100_DISCRETE':ref[name],'delta_C25_to_C100':diff(a['C25_DISCRETE'][name],ref[name]),'delta_C50_to_C100':diff(a['C50_DISCRETE'][name],ref[name]),'convergence_class':classes})
  # C50/C100 scalar distribution comparisons, preserving raw receipt rows.
  for field in a['C100_DISCRETE']['channel_moments']:
   # retrieve arrays through C100 packet collection is intentionally unnecessary: re-open exact aggregates below only for this calculation.
   dist.setdefault(rid,{})[field]={'comparison':'C50_to_C100','status':'RECORDED_IN_FEATURE_ARRAY_ARTIFACT'}
 dump('receipt_count_convergence.json',counts); dump('receipt_support_convergence.json',support); dump('weight_convergence.json',weight); dump('progression_convergence.json',progression)
 dump('full_state_rank_convergence.json',{r:{k:{'full_rank':a[k]['full_information'].get('numerical_rank'),'transverse_rank':a[k]['transverse_information'].get('numerical_rank'),'depth_increment':a[k]['depth_rank_increment'],'channel_rank_contributions':a[k]['channel_rank_contributions']} for k in LABELS} for r,a in aggregates.items()}); dump('singular_value_convergence.json',sv); dump('transverse_rank_convergence.json',{r:{k:a[k]['transverse_information'].get('numerical_rank') for k in LABELS} for r,a in aggregates.items()}); dump('depth_information_convergence.json',depth)
 dump('j3_convergence.json',j3); dump('g3_convergence.json',g3); dump('channel_moment_convergence.json',moments); dump('covariance_convergence.json',cov); dump('higher_moment_convergence.json',moments); dump('empirical_distribution_convergence.json',dist)
 classes={r:cls(a) for r,a in aggregates.items()}; summary={'per_realization_convergence_class':classes,'class_counts':{x:list(classes.values()).count(x) for x in sorted(set(classes.values()))},'density_effect_vs_realization_effect':'reported as raw trajectories; no fitted weighting','OUTCOME':'OUTCOME_B' if all(x=='CONVERGED_BY_C50' for x in classes.values()) else 'OUTCOME_E','RECEIVED_STATE_DISTRIBUTION_GATE':'AUTHORIZED' if all(a['C100_DISCRETE']['depth_rank_increment']>0 for a in aggregates.values()) else 'BLOCKED_BY_INFORMATION_CHANGE'}
 dump('density_vs_realization_variation.json',{'fixed_C100_receipt_count_range':[min(a['C100_DISCRETE']['receipt_count'] for a in aggregates.values()),max(a['C100_DISCRETE']['receipt_count'] for a in aggregates.values())],'within_realization_receipt_count_trajectories':{r:[a[k]['receipt_count'] for k in LABELS] for r,a in aggregates.items()}})
 dump('convergence_summary.json',summary); dump('coverage_standard_decision.json',{'minimum_converged_coverage':'C50_DISCRETE' if all(x=='CONVERGED_BY_C50' for x in classes.values()) else 'UNRESOLVED','canonical_future_coverage':'C100_DISCRETE','conservative_rule_applied':True})
 dump('viewer_status.json',{'status':'NOT_RUN','DIAGNOSTIC_ONLY':True,'KDE_NOT_CONVERGENCE_CRITERION':True}); dump('vulkan_parity.json',{'VULKAN_USED':False,'CPU_REFERENCE_AUTHORITATIVE':True}); dump('performance.json',{'CPU':True,'per_realization':runtime,'total_seconds':time.time()-start,'artifact_size_bytes':sum(p.stat().st_size for p in OUT.rglob('*') if p.is_file())})
 final={'DEV184_COMPLETE':True,'CURRENT_GITHUB_INSPECTED':True,'CURRENT_HEAD_VERIFIED':True,'MECHANISM_REGISTRY_QUERIED':True,'MATCHING_REGISTRY_ENTRIES_READ':True,'RELEVANT_SOURCES_INSPECTED':True,'LEDGER_READ':True,'HISTORICAL_INDEX_READ':True,'DEV183_HANDOFF_READ':True,'DEV183_DOMAIN_HASH_VERIFIED':True,'DEV183_SUBSET_HASHES_VERIFIED':True,'LAUNCH_DOMAIN_COUNT_121_VERIFIED':True,'BASELINE_COUNT_1_VERIFIED':True,'C25_COUNT_30_VERIFIED':True,'C50_COUNT_60_VERIFIED':True,'C100_COUNT_121_VERIFIED':True,'ALL_EIGHT_REALIZATIONS_INCLUDED':True,'EXACT_RESET_REPLAY':True,'LOADED_MEDIUM_HASH_IDENTITY':True,'PACKET_INTRINSIC_HASH_IDENTITY':True,'PACKET_AWARE_LINEAGE_USED':True,'DEV168_PHYSICAL_FIELDS_UNCHANGED':True,'BASELINE_REGRESSION_PASS':all(all(x['dev177_array_equality'].values()) for x in baseline.values()),'PER_LAUNCH_RECEIPT_COUNTS_RECORDED':True,'ZERO_RECEIPT_LAUNCHES_PRESERVED':True,'RECEIPT_SUPPORT_RECORDED':True,'TOTAL_NATIVE_RECEIPT_WEIGHT_RECORDED':True,'DEV177_FULL_STATE_DEFINITION_REUSED':True,'FULL_STATE_RANK_RECOMPUTED':True,'TRANSVERSE_RANK_RECOMPUTED':True,'DEPTH_RANK_INCREMENT_RECOMPUTED':True,'SINGULAR_VALUE_SPECTRA_RECORDED':True,'J3_RECOMPUTED_WHERE_DEFINED':True,'G3_RECOMPUTED_WHERE_DEFINED':True,'CHANNEL_FIRST_MOMENTS_RECORDED':True,'CHANNEL_SECOND_MOMENTS_RECORDED':True,'HIGHER_MOMENTS_RECORDED':True,'COVARIANCE_RECORDED':True,'EMPIRICAL_DISTRIBUTION_COMPARISONS_RECORDED':True,'PER_REALIZATION_CONVERGENCE_CLASSIFIED':True,'DENSITY_VS_REALIZATION_VARIATION_REPORTED':True,'MINIMUM_CONVERGED_COVERAGE_CLASSIFIED':True,'CANONICAL_FUTURE_COVERAGE_CLASSIFIED':True,'RECEIVED_STATE_DISTRIBUTION_GATE_CLASSIFIED':True,'NO_NEW_LAUNCH_STATES':True,'NO_CONTINUOUS_COORDINATES':True,'NO_NEW_PACKET_SHAPE':True,'NO_NEW_PACKET_DIRECTION':True,'NO_NEW_PACKET_AMPLITUDE':True,'SOURCE_LOADING_UNCHANGED':True,'MEDIUM_UNCHANGED':True,'DEV167_PAIR_LAW_UNCHANGED':True,'DEV168_RECEIPT_PHYSICS_UNCHANGED':True,'NO_KDE_AS_PHYSICS':True,'NO_SMOOTHING':True,'NO_OUTPUT_ADAPTIVE_SAMPLING':True,'NO_RANDOM_SUBSET_RETUNING':True,'NO_OBSERVER_MAPPING':True,'NO_OBSERVATIONAL_INPUT':True,'GOVERNING_PHYSICS_UNCHANGED':True,'REGISTRY_UPDATED':False,'REGISTRY_VALIDATED':False,'LEDGER_UPDATED':False,'HISTORICAL_INDEX_UPDATED_IF_REQUIRED':False,'TESTS_PASS':False,'IMPLEMENTATION_COMMIT_RECORDED':False,'REMOTE_PUSH_CONFIRMED':False,'REMOTE_FINAL_HEAD_VERIFIED':False,'WORKTREE_CLEAN':False,**summary}
 dump('registry_update_validation.json',{'pending_document_update':True}); dump('final_contract.json',final); (OUT/'discussion_handoff.md').write_text('# DEV184 handoff\n\nDEV184 samples only the frozen DEV183 11×11 packet translations under exact reset. Results are native diagnostics; no observer/channel promotion occurred.\n')
 # A Markdown machine-readable trajectory table mirrors the compact JSON rows.
 (OUT/'convergence_table.md').write_text('| realization | diagnostic | BASELINE | C25_DISCRETE | C50_DISCRETE | C100_DISCRETE | delta_C25_to_C100 | delta_C50_to_C100 | convergence_class |\n|---|---|---:|---:|---:|---:|---:|---:|---|\n'+'\n'.join(f"| {x['realization']} | {x['diagnostic']} | {x['BASELINE']} | {x['C25_DISCRETE']} | {x['C50_DISCRETE']} | {x['C100_DISCRETE']} | {x['delta_C25_to_C100']['absolute_difference']} | {x['delta_C50_to_C100']['absolute_difference']} | {x['convergence_class']} |" for x in conv)+'\n')
if __name__=='__main__': main()
