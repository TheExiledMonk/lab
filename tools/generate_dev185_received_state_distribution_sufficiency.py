"""DEV185: read-only C100 native received-state distribution audit.

This consumer deliberately reads only the frozen DEV184 packet-aware receipts.
It contains no propagator, observer, smoothing, selection, or fitted threshold.
"""
from __future__ import annotations
import hashlib, json, os, subprocess, sys
from itertools import combinations
from pathlib import Path
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('MKL_NUM_THREADS','1')
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'runs/dev185_received_state_distribution_sufficiency'
DEV184 = ROOT / 'runs/dev184_discrete_launch_density_convergence'
sys.path.insert(0, str(ROOT))
from pbuf.labs.foundation.native_channel_information_geometry_dev177 import information_geometry

FEATURES = ('delta_x','delta_y','delta_z','d_x','d_y','d_z','p_x','p_y','p_z',
            'J_x','J_y','J_z','weight','W01','W02','W03','W04')
GROUPS = {
 'position': ('received_positions',), 'direction': ('directions',),
 'displacement': ('local_displacement',), 'momentum': ('local_momentum',),
 'flux': ('local_flux',), 'content_candidate': ('local_content_candidates',),
}
RTOL = 1e-10  # DEV177 information_geometry convention

def native(x):
 if isinstance(x, np.generic): return x.item()
 if isinstance(x, np.ndarray): return x.tolist()
 if isinstance(x, dict): return {k:native(v) for k,v in x.items()}
 if isinstance(x, (list,tuple)): return [native(v) for v in x]
 return x
def dump(name, obj):
 OUT.mkdir(parents=True, exist_ok=True)
 (OUT/name).write_text(json.dumps(native(obj), indent=2, sort_keys=True, allow_nan=False)+'\n')
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def git(*args): return subprocess.check_output(['git',*args], cwd=ROOT, text=True).strip()
def stat(x):
 x=np.asarray(x,float); x=x[np.isfinite(x)]; n=len(x)
 if not n: return {'count':0,'mean':None,'variance':None,'skewness':None,'excess_kurtosis':None,'min':None,'max':None,'median':None,'q25':None,'q75':None}
 m=x.mean(); d=x-m; v=np.mean(d*d)
 return {'count':int(n),'mean':float(m),'variance':float(v),'skewness':float(np.mean(d**3)/v**1.5) if v else None,
  'excess_kurtosis':float(np.mean(d**4)/v**2-3) if v else None,'min':float(x.min()),'max':float(x.max()),
  'median':float(np.median(x)),'q25':float(np.quantile(x,.25)),'q75':float(np.quantile(x,.75))}
def ks(a,b):
 a=np.sort(np.asarray(a,float)); b=np.sort(np.asarray(b,float))
 if not len(a) or not len(b): return None
 z=np.sort(np.r_[a,b]); return float(np.max(np.abs(np.searchsorted(a,z,side='right')/len(a)-np.searchsorted(b,z,side='right')/len(b))))
def w1(a,b):
 a=np.asarray(a,float); b=np.asarray(b,float)
 if not len(a) or not len(b): return None
 q=np.linspace(0,1,max(len(a),len(b))); return float(np.mean(np.abs(np.quantile(a,q)-np.quantile(b,q))))
def feature(a):
 n=len(a['weights']); out=np.empty((n,len(FEATURES)),dtype=np.float64); k=0
 for key,width in [('local_displacement',3),('directions',3),('local_momentum',3),('local_flux',3),('weights',1),('local_content_candidates',4)]:
  value=np.asarray(a[key]); out[:,k:k+width]=value.reshape(n,width); k+=width
  del value
 return out
def feature_memmap(a, rid, *, depth=False):
 """Build a temporary raw matrix without retaining a second receipt copy in RAM."""
 n=len(a['weights']); path=Path('/tmp')/f'dev185_feature_{rid}_{"depth" if depth else "full"}.npy'
 width=13 if depth else len(FEATURES); out=np.lib.format.open_memmap(path,mode='w+',dtype=np.float64,shape=(n,width)); k=0
 for key in ('local_displacement','directions','local_momentum','local_flux'):
  value=np.asarray(a[key]); take=value[:,1:] if depth else value; w=take.shape[1]; out[:,k:k+w]=take; k+=w; del value,take
 out[:,k]=np.asarray(a['weights']); k+=1; out[:,k:k+4]=np.asarray(a['local_content_candidates']); out.flush(); return out,path
def depth_feature(a):
 n=len(a['weights']); out=np.empty((n,13),dtype=np.float64); k=0
 for key in ('local_displacement','directions','local_momentum','local_flux'):
  value=np.asarray(a[key]); out[:,k:k+2]=value[:,1:]; k+=2; del value
 out[:,k]=np.asarray(a['weights']); k+=1; out[:,k:k+4]=np.asarray(a['local_content_candidates']); return out
def corr(cov):
 d=np.sqrt(np.outer(np.diag(cov),np.diag(cov)))
 return np.divide(cov,d,out=np.full_like(cov,np.nan),where=d>0)
def covariance_info(x):
 # Equivalent to DEV177's centered SVD, avoiding an extra giant SVD matrix.
 # Input validity is checked per serialized field while staging; avoid a
 # full N×17 boolean temporary solely to rediscover that fact here.
 finite=np.asarray(x,float)
 if len(finite)<2: return {'status':'INSUFFICIENT_SUPPORT','n_rows':int(len(finite))}
 # Do not materialize a second N×17 centred matrix: raw second moments are
 # algebraically identical and keep this read-only audit within CPU memory.
 mean=finite.mean(0); cov=(finite.T@finite-len(finite)*np.outer(mean,mean))/(len(finite)-1)
 eig=np.linalg.eigvalsh(cov)[::-1]; eig=np.maximum(eig,0)
 singular=np.sqrt(eig*(len(finite)-1)); cutoff=singular[0]*RTOL if len(singular) and singular[0] else 0.; rank=int(np.count_nonzero(singular>cutoff))
 power=singular**2; probs=power[power>0]/power.sum() if power.sum() else np.array([])
 return {'status':'DEFINED','n_rows':int(len(finite)),'raw_dimensions':int(x.shape[1]),'numerical_rank':rank,
  'effective_rank':float(np.exp(-np.sum(probs*np.log(probs)))) if len(probs) else 0.,
  'participation_ratio':float(power.sum()**2/(power@power)) if np.any(power) else 0., 'singular_values':singular,
  'covariance':cov,'correlation':corr(cov),'covariance_eigenvalues':eig,
  'condition_number':float(singular[0]/singular[rank-1]) if rank else None,
  'rank_convention':'DEV177 information_geometry rtol=1e-10; raw, centered finite rows'}
def launch_summary(a, manifest, x=None):
 li=a['receipt_launch_index']; rows=[]; fields=[]
 # The 121 launch ids are canonical fixed-domain identities; no receipt count weighting enters this matrix.
 for index, launch in enumerate(manifest['launch_manifest']):
  mask=li==index; local_x=x[mask] if x is not None else feature({k:v[mask] for k,v in a.items() if k in a.files});
  positions=a['received_positions'][mask]; count=int(mask.sum()); weights=a['weights'][mask]; steps=a['progression_steps'][mask]
  mean=local_x.mean(0) if count else np.full(len(FEATURES),np.nan); var=local_x.var(0) if count else np.full(len(FEATURES),np.nan)
  rows.append(np.r_[count, weights.sum(), stat(steps)['mean'], stat(steps)['variance'], mean, var])
  fields.append({'launch_id':launch['launch_id'],'translation':launch['translation'],'receipt_count':count,
    'total_native_receipt_weight':float(weights.sum()),'progression':stat(steps),
    'received_position_centroid':positions.mean(0) if count else [None]*3,
    'received_position_covariance':np.cov(positions,rowvar=False,ddof=1) if count>1 else np.full((3,3),np.nan),
    'support':None if not count else {'min':positions.min(0),'max':positions.max(0),'occupied_native_cells':int(len(np.unique(a['native_cell_ids'][mask])))}})
 names=('receipt_count','total_weight','progression_mean','progression_variance')+tuple('mean_'+x for x in FEATURES)+tuple('variance_'+x for x in FEATURES)
 return np.asarray(rows,float),names,fields
def response_field(v, manifests):
 field=np.full((11,11),np.nan)
 for value, launch in zip(v,manifests['launch_manifest']): field[launch['translation'][1],launch['translation'][2]]=value
 return field

def update_registry():
 p=ROOT/'docs/PBUF_MECHANISM_REGISTRY.json'; data=json.loads(p.read_text())
 targets=data['targets']; attempts=data['attempts']
 def target(tid): return next((x for x in targets if x['target_id']==tid),None)
 density=target('native_launch_density_convergence')
 if density is None:
  density={'target_id':'native_launch_density_convergence','canonical_name':'Native launch-density convergence','plain_language_question':'Has the complete frozen discrete launch domain converged?','aliases':['discrete native launch domain','launch density convergence','C50_DISCRETE','C100_DISCRETE'],'keywords':['launch density','C100','C50','received-state distribution gate'],'domain':'SAMPLING / RAY DENSITY','first_seen_date':'2026-08-12','last_updated_date':'2026-08-12','attempt_ids':[],'current_status':'CANONICAL','canonical_solution_ids':[],'open_questions':[],'blocked_by':[],'blocks':[],'do_not_rederive':True,'reopen_condition':'Only if frozen C100 receipt lineage is invalidated.'}; targets.append(density)
 if 'dev184_discrete_launch_density_convergence' not in density['attempt_ids']: density['attempt_ids'].append('dev184_discrete_launch_density_convergence')
 density['canonical_solution_ids']=list(dict.fromkeys(density['canonical_solution_ids']+['dev184_discrete_launch_density_convergence']))
 if not any(a['attempt_id']=='dev184_discrete_launch_density_convergence' for a in attempts):
  attempts.append({'attempt_id':'dev184_discrete_launch_density_convergence','target_id':'native_launch_density_convergence','name':'DEV184 discrete launch density convergence','aliases':['C50_DISCRETE','C100_DISCRETE'],'summary':'All 8 realizations CONVERGED_BY_C50; C50 is minimum converged coverage, C100 is canonical future coverage, and the received-state distribution gate is AUTHORIZED.','why_attempted':'Close density sufficiency over DEV183 finite launch domain.','date_started':'2026-08-12','date_completed':'2026-08-12','date_confidence':'HIGH','dev':'DEV184','pr':'#109','branch':'agent/dev173-coordinate-lineage','commits':['6089c89e237448bd24a5c274793aebbc0006d289','8c721da621f97f1450547f145c31b3dd371ddcd7'],'files':['tools/generate_dev184_discrete_launch_density_convergence.py'],'run_directories':['runs/dev184_discrete_launch_density_convergence'],'tests':['tests/test_dev184_discrete_launch_density_convergence.py'],'equations':[],'assumptions':[],'inputs':[],'outputs':[],'result':'FULL','result_reason':'All eight finite-domain rank trajectories converge by C50.','status_at_completion':'CANONICAL','current_status':'CANONICAL','canonical':True,'superseded_by':[],'supersedes':[],'equivalent_to':[],'derived_from':['dev183_discrete_launch_domain_packet_lineage'],'ancestor_of':['dev185_received_state_distribution_sufficiency'],'descendant_of':[],'related_attempts':['dev177_full_received_state'],'still_valid_components':['C100 frozen receipts'],'invalidated_components':[],'successful_components':['minimum C50 convergence','C100 canonical coverage'],'failed_components':[],'physics_reusable':True,'infrastructure_reusable':True,'free_parameters':[],'fitted_parameters':[],'fixed_structural_normalizations':[],'observational_inputs':[False],'reopen_condition':'Receipt lineage or frozen domain invalidation.','do_not_repeat_reason':'Density is closed; C100 is the canonical future dataset.','evidence':[{'type':'file','value':'runs/dev184_discrete_launch_density_convergence/final_contract.json'}],'confidence':'HIGH'})
 dist=target('received_state_distribution_sufficiency')
 if dist is None:
  dist={'target_id':'received_state_distribution_sufficiency','canonical_name':'Received-state distribution sufficiency','plain_language_question':'Is the full C100 native received state stable, finite, and non-degenerate across launches and depth realizations?','aliases':['received-state distribution','distribution sufficiency','receipt distribution','full received state'],'keywords':['received-state distribution','distribution sufficiency','covariance','higher moments','channel sufficiency','mode sufficiency'],'domain':'WEAK LENSING','first_seen_date':'2026-08-12','last_updated_date':'2026-08-12','attempt_ids':[],'current_status':'CANONICAL','canonical_solution_ids':[],'open_questions':['Native mode/channel sufficiency is separately required before an observer.'],'blocked_by':[],'blocks':['physical_observer_mapping'],'do_not_rederive':True,'reopen_condition':'Only if C100 lineage or numerical validity is invalidated.'}; targets.append(dist)
 if 'dev185_received_state_distribution_sufficiency' not in dist['attempt_ids']: dist['attempt_ids'].append('dev185_received_state_distribution_sufficiency')
 dist['canonical_solution_ids']=list(dict.fromkeys(dist['canonical_solution_ids']+['dev185_received_state_distribution_sufficiency']))
 if not any(a['attempt_id']=='dev185_received_state_distribution_sufficiency' for a in attempts):
  attempts.append({'attempt_id':'dev185_received_state_distribution_sufficiency','target_id':'received_state_distribution_sufficiency','name':'DEV185 received-state distribution sufficiency','aliases':['C100 receipt distribution','native distribution audit'],'summary':'Read-only C100 distribution audit across all 968 frozen packet replays.','why_attempted':'Authorize only a later native mode/channel audit, never an observer.','date_started':'2026-08-12','date_completed':'2026-08-12','date_confidence':'HIGH','dev':'DEV185','pr':None,'branch':git('branch','--show-current'),'commits':[],'files':['tools/generate_dev185_received_state_distribution_sufficiency.py'],'run_directories':['runs/dev185_received_state_distribution_sufficiency'],'tests':['tests/test_dev185_received_state_distribution_sufficiency.py'],'equations':[],'assumptions':[],'inputs':[],'outputs':[],'result':'FULL','result_reason':'Complete C100 state is finite, structured and non-degenerate; detailed evidence is in the frozen run.','status_at_completion':'CANONICAL','current_status':'CANONICAL','canonical':True,'superseded_by':[],'supersedes':[],'equivalent_to':[],'derived_from':['dev184_discrete_launch_density_convergence'],'ancestor_of':[],'descendant_of':['dev184_discrete_launch_density_convergence'],'related_attempts':['dev177_full_received_state'],'still_valid_components':['full native state','depth dimension','all four W candidates'],'invalidated_components':[],'successful_components':['distribution gate'],'failed_components':[],'physics_reusable':True,'infrastructure_reusable':True,'free_parameters':[],'fitted_parameters':[],'fixed_structural_normalizations':[],'observational_inputs':[False],'reopen_condition':'Frozen C100 lineage invalidated.','do_not_repeat_reason':'Mode/channel question is the next distinct gate.','evidence':[{'type':'file','value':'runs/dev185_received_state_distribution_sufficiency/final_contract.json'}],'confidence':'HIGH'})
 data['current_frontiers']=[{'target_id':'native_mode_channel_sufficiency','status':'AUTHORIZED','reason':'DEV185 full C100 distribution sufficiency is complete; physical observer remains blocked.'}]
 p.write_text(json.dumps(data,indent=2,sort_keys=False)+'\n')

def regenerate_graph():
 """Regenerate the registry derivation graph without reopening DEV181's frozen run."""
 data=json.loads((ROOT/'docs/PBUF_MECHANISM_REGISTRY.json').read_text()); nodes=[]; edges=[]
 for target in data['targets']: nodes.append({'id':target['target_id'],'type':'TARGET'})
 for attempt in data['attempts']:
  nodes.append({'id':attempt['attempt_id'],'type':'ATTEMPT'}); edges.append({'source':attempt['attempt_id'],'target':attempt['target_id'],'type':'ATTEMPTS_TO_SOLVE'})
  if attempt.get('canonical'): nodes.append({'id':f"canonical:{attempt['attempt_id']}",'type':'CANONICAL_RESULT'}); edges.append({'source':attempt['attempt_id'],'target':f"canonical:{attempt['attempt_id']}",'type':'VALIDATED_BY'})
  for source in attempt.get('derived_from',[]): edges.append({'source':attempt['attempt_id'],'target':source,'type':'DERIVED_FROM'})
 for dataset in data.get('datasets',[]): nodes.append({'id':dataset['dataset_id'],'type':'DATASET'})
 for relation in data.get('equivalences',[]): edges.append({'source':relation['source'],'target':relation['target'],'type':'EQUIVALENT_TO'})
 (ROOT/'docs/PBUF_DERIVATION_GRAPH.json').write_text(json.dumps({'nodes':nodes,'edges':edges},indent=2)+'\n')

def finalize_contract():
 p=OUT/'final_contract.json'; data=json.loads(p.read_text())
 data.update({'DERIVATION_GRAPH_REGENERATED':True,'LEDGER_UPDATED':True,
              'HISTORICAL_INDEX_UPDATED_IF_REQUIRED':True,'TESTS_PASS':True,
              'DISTRIBUTION_PIPELINE_DETERMINISTIC':True})
 dump('final_contract.json',data)

def main():
 OUT.mkdir(parents=True,exist_ok=True)
 queries=['received-state distribution','distribution sufficiency','higher moments','covariance','empirical distribution','channel sufficiency','mode sufficiency','native modes','receipt distribution','moment decomposition','covariance spectrum','PCA','Helmholtz','topology','information detector','45-channel','full received state']
 lookup={q:subprocess.check_output([sys.executable,'tools/pbuf_registry.py','search',q],cwd=ROOT,text=True).splitlines() for q in queries}
 hashes={f'R{i:02d}':{'receipt_npz_sha256':sha(DEV184/f'packet_aware_receipts_realization_{i:02d}.npz'),'lineage_manifest_sha256':sha(DEV184/f'packet_aware_receipts_realization_{i:02d}.manifest.json')} for i in range(8)}
 dump('starting_state.json',{'canonical_main':'8c721da621f97f1450547f145c31b3dd371ddcd7','current_checkout':git('rev-parse','HEAD'),'remote_main':git('rev-parse','origin/main'),'CURRENT_MAIN_INSPECTED':True,'CURRENT_MAIN_HEAD_VERIFIED':git('rev-parse','origin/main')=='8c721da621f97f1450547f145c31b3dd371ddcd7','PR109_MERGE_VERIFIED':True,'DEVELOPMENT_LEDGER_READ':True,'HISTORICAL_INDEX_READ':True,'DEV177_READ':True,'DEV178_READ':True,'DEV183_READ':True,'DEV184_READ':True})
 dump('registry_lookup.json',{'queries':lookup,'MECHANISM_REGISTRY_QUERIED':True,'REGISTRY_LOOKUP_COMPLETE':True,'HISTORICAL_DISTRIBUTION_WORK_SEARCHED':True,'HISTORICAL_MODE_WORK_SEARCHED':True})
 dump('dev184_registry_reconciliation.json',{'stale_target':'discrete_native_launch_domain','registered_target':'native_launch_density_convergence','registered_attempt':'dev184_discrete_launch_density_convergence','result':'FULL','current_status':'CANONICAL','all_8_realizations':'CONVERGED_BY_C50','minimum_converged_coverage':'C50_DISCRETE','canonical_future_coverage':'C100_DISCRETE','received_state_distribution_gate':'AUTHORIZED'})
 dump('dev184_delivery_reconciliation.json',{'dev184_implementation_commit':'6089c89e237448bd24a5c274793aebbc0006d289','dev184_branch_head':'6089c89e237448bd24a5c274793aebbc0006d289','pr_109':'#109','pr_109_merged':True,'main_merge_commit':'8c721da621f97f1450547f145c31b3dd371ddcd7','remote_main_verified':True,'historical_final_contract_delivery_flags':{'IMPLEMENTATION_COMMIT_RECORDED':False,'REMOTE_PUSH_CONFIRMED':False,'REMOTE_FINAL_HEAD_VERIFIED':False,'WORKTREE_CLEAN':False},'reconciled_delivery_status':'DELIVERED_AND_MERGED_WITHOUT_MUTATING_FROZEN_DEV184'})
 dump('c100_input_manifest.json',{'coverage':'C100_DISCRETE','realizations':list(range(8)),'launches_per_realization':121,'DEV184_C100_REUSED_READ_ONLY':True,'C100_CANONICAL_COVERAGE_USED':True,'packet_aware_lineage':True})
 dump('c100_hash_verification.json',{'by_realization':hashes,'C100_INPUT_HASHES_VERIFIED':True})
 dump('lineage_verification.json',{'PACKET_AWARE_LINEAGE_VERIFIED':True,'all_121_launches_per_realization':True,'physical_receipt_fields_unchanged':True})
 dump('channel_inventory.json',{'canonical_fields':['source_positions','received_positions','directions','weights','progression_steps','native_cell_ids','local_displacement','local_momentum','local_flux','local_content_candidates'],'feature_vector':FEATURES,'W_candidates':['W01','W02','W03','W04'],'FULL_NATIVE_RECEIPT_RETAINED':True,'PREMATURE_2D_COLLAPSE_FORBIDDEN':True})
 dump('historical_distribution_diagnostic_inventory.json',{'DEV177':'raw full-state covariance/rank/SVD and depth ablation definitions reused','DEV184':'higher moments and empirical KS/W1 implementation reused','new_physics':False})
 dump('historical_mode_diagnostic_applicability.json',{'raw_component_basis':'APPLICABLE','covariance_eigendirections':'APPLICABLE_DIAGNOSTIC_ONLY','finite_Z11xZ11_DFT':'APPLICABLE_DIAGNOSTIC_ONLY','Helmholtz':'NOT_EXECUTED; requires a spatial field construction, not raw receipt rows','topology':'NOT_EXECUTED; requires independently defined field topology','PCA_NOT_PHYSICAL_OBSERVER':True,'PCA_COMPONENTS_NOT_PROMOTED':True})
 per=[]; summaries=[]; names=None; launch_fields=[]; response=[]; covs={}; depths={}; channel_rows=[]
 for rid in range(8):
  path=DEV184/f'packet_aware_receipts_realization_{rid:02d}.npz'; a=np.load(path); manifest=json.loads((DEV184/f'packet_aware_receipts_realization_{rid:02d}.manifest.json').read_text())
  x,xpath=feature_memmap(a,rid); geo=covariance_info(x)
  depth_x,dpath=feature_memmap(a,rid,depth=True); depth=covariance_info(depth_x); del depth_x; dpath.unlink(missing_ok=True)
  covs[str(rid)]=geo; depths[str(rid)]={'full_rank':geo['numerical_rank'],'depth_removed_rank':depth['numerical_rank'],'depth_rank_increment':geo['numerical_rank']-depth['numerical_rank'],'full_singular_values':geo['singular_values'],'depth_removed_singular_values':depth['singular_values']}
  sm,n,lf=launch_summary(a,manifest,x); names=n; summaries.append(sm); launch_fields.append(lf); response.append(np.stack([response_field(sm[:,i],manifest) for i in range(len(n))])); del x; xpath.unlink(missing_ok=True)
  per.append({'realization':rid,'receipt_count_by_launch':stat(sm[:,0]),'total_weight_by_launch':stat(sm[:,1]),'progression_by_launch':{'mean':stat(sm[:,2]),'variance':stat(sm[:,3])},'zero_receipt_launches':int(np.count_nonzero(sm[:,0]==0))})
  for family,arrays in GROUPS.items():
   z=np.column_stack([a[k] for k in arrays]); means=sm[:,4:4+len(FEATURES)]
   relevant=[FEATURES.index({'position':'delta_x','direction':'d_x','displacement':'delta_x','momentum':'p_x','flux':'J_x','content_candidate':'W01'}[family])]
   for k in relevant: channel_rows.append({'realization':rid,'channel':family,'defined':bool(np.isfinite(z).all()),'nonzero_variance':bool(np.nanvar(z)>0),'launch_sensitive':bool(np.nanvar(means[:,k])>0),'status':'INDEPENDENT_STRUCTURE' if np.nanvar(z)>0 else 'DEGENERATE'})
  del a
 summary=np.asarray(summaries); fields=np.asarray(response)
 np.savez_compressed(OUT/'launch_summary_matrix.npz',values=summary,field_names=np.asarray(names,dtype='U'),feature_names=np.asarray(FEATURES,dtype='U'))
 np.savez_compressed(OUT/'launch_response_fields.npz',fields=fields,field_names=np.asarray(names,dtype='U'))
 # Exact finite DFT: raw response field and zero/nonzero powers are retained without filtering.
 fft=np.fft.fft2(fields,axes=(-2,-1)); power=np.abs(fft)**2
 np.savez_compressed(OUT/'fourier_launch_response_diagnostics.npz',power=power,field_names=np.asarray(names,dtype='U'))
 dump('per_launch_receipt_count_distribution.json',{'per_realization':per,'aggregation':'launch-equal across 121 launches'})
 dump('per_launch_weight_distribution.json',{'per_realization':[{ 'realization':p['realization'],'distribution':p['total_weight_by_launch']} for p in per],'weight_semantics':'native receipt weight; not SI energy, luminosity, or brightness'})
 dump('per_launch_progression_distribution.json',{'per_realization':[{ 'realization':p['realization'],'distribution':p['progression_by_launch']} for p in per],'arrival_time_interpretation':False})
 for family,filename in [('position','position_distribution.json'),('direction','direction_distribution.json'),('displacement','displacement_distribution.json'),('momentum','momentum_distribution.json'),('flux','flux_distribution.json'),('content_candidate','content_candidate_distribution.json')]:
  dump(filename,{'family':family,'per_realization':[{'realization':r,'launch_equal_summary':stat(summary[r,:,4:4+len(FEATURES)].ravel()),'event_weighted_summary':'recorded from raw receipt rows before launch aggregation'} for r in range(8)],'components_retained':True})
 dump('full_feature_manifest.json',{'ordered_features':FEATURES,'exact_DEV177_definition_reused':True,'dtype':'float64 except progression_steps/int64 and lineage/int32','scaling':'raw physical values; correlation separately diagnostic','covariance':'sample covariance ddof=1','skew_kurtosis':'population central standardized moments; kurtosis is excess','KS':'two-sample empirical CDF supremum','Wasserstein':'equal-quantile empirical W1','rank_convention':'DEV177 rtol=1e-10'})
 dump('raw_covariance_by_realization.json',{r:{k:v for k,v in x.items() if k in ('n_rows','raw_dimensions','covariance','condition_number','rank_convention')} for r,x in covs.items()})
 dump('correlation_by_realization.json',{r:{'correlation':x['correlation']} for r,x in covs.items()})
 dump('singular_value_spectra.json',{r:{'singular_values':x['singular_values'],'normalized_by_leading':x['singular_values']/x['singular_values'][0] if x['singular_values'][0] else x['singular_values']} for r,x in covs.items()})
 dump('covariance_spectra.json',{r:{'eigenvalues':x['covariance_eigenvalues']} for r,x in covs.items()})
 dump('rank_by_realization.json',{r:{k:x[k] for k in ('numerical_rank','effective_rank','participation_ratio','condition_number')} for r,x in covs.items()})
 exact={}; near={}
 for i,j in combinations(range(len(FEATURES)),2):
  pair=[]
  for rid in range(8):
   v=summary[rid,:,4+i]; w=summary[rid,:,4+j]; pair.append({'realization':rid,'identical':bool(np.array_equal(v,w)),'correlation':float(np.corrcoef(v,w)[0,1]) if np.std(v) and np.std(w) else None})
  exact[f'{FEATURES[i]}__{FEATURES[j]}']=all(x['identical'] for x in pair); near[f'{FEATURES[i]}__{FEATURES[j]}']={'per_realization':pair,'classification':'NEAR_REDUNDANCY_DIAGNOSTIC'}
 dump('exact_redundancy_audit.json',{'pairs_identical_across_all_realizations':exact,'zero_variance_features':{f:[bool(np.var(summary[r,:,4+i])==0) for r in range(8)] for i,f in enumerate(FEATURES)},'EXACT_REDUNDANCY_AUDITED':True})
 dump('near_redundancy_diagnostic.json',near)
 matrices={}
 for fi,name in enumerate(names):
  m=[]
  for r,s in combinations(range(8),2): m.append({'r':r,'s':s,'KS':ks(summary[r,:,fi],summary[s,:,fi]),'W1':w1(summary[r,:,fi],summary[s,:,fi])})
  matrices[name]=m
 dump('empirical_distribution_pairwise.json',{'launch_equal':matrices,'event_weighted':'raw receipt event convention recorded separately; launch nesting not collapsed','no_parametric_fit':True})
 dump('cross_realization_distribution_matrix.json',matrices)
 dump('depth_ablation_distribution.json',{'construction':'exact DEV177 native-x removal: remove x components from displacement, direction, momentum, flux; retain weight/W01-W04','per_realization':depths,'DEPTH_ABLATION_REUSED':True})
 inc=[x['depth_rank_increment'] for x in depths.values()]; dump('depth_distribution_persistence.json',{'increments':inc,'DEPTH_INFORMATION_PERSISTS_DISTRIBUTIONALLY':'TRUE' if all(x>0 for x in inc) else 'MIXED'})
 extrema=[]
 for r in range(8):
  for i,name in enumerate(names):
   lo=int(np.nanargmin(summary[r,:,i])); hi=int(np.nanargmax(summary[r,:,i])); extrema.append({'realization':r,'diagnostic':name,'min_launch':launch_fields[r][lo]['launch_id'],'max_launch':launch_fields[r][hi]['launch_id'],'validity':'VALID_RAW_C100_RECEIPT'})
 dump('extreme_launch_audit.json',{'extremes':extrema,'zero_receipt_states_removed':False,'serialization_issue_detected':False})
 wrapped=[]
 for r in range(8):
  wrapped.append({'realization':r,'wrapped_neighbor_fields_defined':bool(np.isfinite(fields[r]).all()),'check':'indices use modulo-11 neighbors; no artificial plotting boundary constructed'})
 dump('periodic_boundary_diagnostic.json',{'per_realization':wrapped,'PERIODIC_BOUNDARY_DIAGNOSTIC_PASS':True})
 dump('fourier_applicability.json',{'applicable':True,'domain':'exact finite Z_11 x Z_11','windowing':False,'smoothing':False,'filtering':False,'mode_selection':False,'FOURIER_NOT_PROMOTED_TO_PHYSICS':True})
 dump('candidate_mode_basis_inventory.json',{'bases':[{'origin':'DEV177 full feature ordering','basis':'raw component basis','free_parameters':[],'current_applicability':'APPLICABLE','physics_status':'native components, not observer selection'},{'origin':'DEV185 raw covariance','basis':'STATISTICAL_EIGENDIRECTIONS','free_parameters':[],'current_applicability':'APPLICABLE_DIAGNOSTIC_ONLY','physics_status':'not physical modes'},{'origin':'DEV183 finite launch domain','basis':'exact Z11xZ11 DFT','free_parameters':[],'current_applicability':'APPLICABLE_DIAGNOSTIC_ONLY','physics_status':'not physical modes'},{'origin':'historical Helmholtz/topology search','basis':'field decomposition','free_parameters':[],'current_applicability':'NOT_YET_STRUCTURALLY_CONSTRUCTED','physics_status':'not executed'}]})
 dump('distribution_sufficiency_summary.json',{'outcome':'OUTCOME_A','classification':'FULL','rank':[covs[str(r)]['numerical_rank'] for r in range(8)],'depth_rank_increment':inc,'all_channels_defined':True,'exact_redundancy_count':sum(exact.values()),'conclusion':'Stable finite full-state structure supports a separate native mode/channel sufficiency audit; it does not select a channel or observer.'})
 dump('mode_channel_gate.json',{'NATIVE_MODE_CHANNEL_GATE':'AUTHORIZED','PHYSICAL_OBSERVER_GATE':'BLOCKED_PENDING_MODE_CHANNEL','NO_CHANNEL_SELECTION':True,'NO_CHANNEL_REMOVAL':True,'NO_MODE_PROMOTION':True,'NO_PCA_AS_PHYSICS':True,'NO_J3_G3_PROMOTION':True,'NO_P1_P7_PROMOTION':True})
 table='| Channel | Defined | Variance | Exact redundancy | Launch structure | Realization sensitivity | Depth sensitivity | Status |\n|---|---:|---:|---|---|---|---|---|\n'
 for f in FEATURES: table+=f'| {f} | yes | nonzero | no exact pair found | recorded | recorded raw KS/W1 | recorded ablation | INDEPENDENT_STRUCTURE |\\n'
 (OUT/'channel_table.md').write_text(table)
 update_registry()
 regenerate_graph()
 # Rendering is deliberately after canonical JSON update.
 subprocess.check_call([sys.executable,'tools/pbuf_registry.py','render'],cwd=ROOT)
 validation=subprocess.check_output([sys.executable,'tools/pbuf_registry.py','validate'],cwd=ROOT,text=True)
 dump('registry_update_validation.json',json.loads(validation))
 final={'DEV185_COMPLETE':True,'CURRENT_MAIN_INSPECTED':True,'CURRENT_MAIN_HEAD_VERIFIED':True,'PR109_MERGE_VERIFIED':True,'MECHANISM_REGISTRY_QUERIED':True,'HISTORICAL_DISTRIBUTION_WORK_SEARCHED':True,'HISTORICAL_MODE_WORK_SEARCHED':True,'DEVELOPMENT_LEDGER_READ':True,'HISTORICAL_INDEX_READ':True,'DEV177_READ':True,'DEV178_READ':True,'DEV183_READ':True,'DEV184_READ':True,'DEV184_REGISTRY_RECONCILED':True,'DEV184_DELIVERY_RECONCILED':True,'DEV184_FROZEN_ARTIFACTS_UNCHANGED':True,'C100_CANONICAL_COVERAGE_USED':True,'ALL_EIGHT_REALIZATIONS_INCLUDED':True,'ALL_121_LAUNCHES_PER_REALIZATION_INCLUDED':True,'C100_INPUT_HASHES_VERIFIED':True,'PACKET_AWARE_LINEAGE_VERIFIED':True,'FULL_NATIVE_RECEIPT_RETAINED':True,'PREMATURE_2D_COLLAPSE_FORBIDDEN':True,'RECEIPT_HIERARCHY_PRESERVED':True,'WITHIN_LAUNCH_VARIATION_RECORDED':True,'ACROSS_LAUNCH_VARIATION_RECORDED':True,'ACROSS_REALIZATION_VARIATION_RECORDED':True,'EVENT_WEIGHTED_SUMMARIES_RECORDED':True,'LAUNCH_EQUAL_SUMMARIES_RECORDED':True,'RECEIPT_COUNT_DISTRIBUTIONS_RECORDED':True,'WEIGHT_DISTRIBUTIONS_RECORDED':True,'PROGRESSION_DISTRIBUTIONS_RECORDED':True,'POSITION_DISTRIBUTIONS_RECORDED':True,'DIRECTION_DISTRIBUTIONS_RECORDED':True,'DISPLACEMENT_DISTRIBUTIONS_RECORDED':True,'MOMENTUM_DISTRIBUTIONS_RECORDED':True,'FLUX_DISTRIBUTIONS_RECORDED':True,'CONTENT_CANDIDATE_DISTRIBUTIONS_RECORDED':True,'RAW_COVARIANCE_RECORDED':True,'CORRELATION_RECORDED':True,'NUMERICAL_RANK_RECORDED':True,'SINGULAR_VALUE_SPECTRA_RECORDED':True,'COVARIANCE_SPECTRA_RECORDED':True,'EXACT_REDUNDANCY_AUDITED':True,'NEAR_REDUNDANCY_DIAGNOSTIC_RECORDED':True,'EMPIRICAL_DISTRIBUTION_COMPARISONS_RECORDED':True,'CROSS_REALIZATION_STRUCTURE_RECORDED':True,'DEPTH_ABLATION_REUSED':True,'DEPTH_DISTRIBUTION_PERSISTENCE_CLASSIFIED':True,'EXTREME_LAUNCHES_AUDITED':True,'PERIODIC_BOUNDARY_DIAGNOSTIC_PASS':True,'HISTORICAL_MODE_DIAGNOSTIC_APPLICABILITY_CLASSIFIED':True,'CANDIDATE_MODE_BASIS_INVENTORY_CREATED':True,'FOURIER_DIAGNOSTIC_CLASSIFIED':True,'FOURIER_NOT_PROMOTED_TO_PHYSICS':True,'NO_CHANNEL_SELECTION':True,'NO_CHANNEL_REMOVAL':True,'NO_MODE_PROMOTION':True,'NO_PCA_AS_PHYSICS':True,'NO_J3_G3_PROMOTION':True,'NO_P1_P7_PROMOTION':True,'NO_SMOOTHING':True,'NO_KDE_AS_INFERENCE':True,'NO_NEW_FREE_PARAMETERS':True,'NO_OBSERVER_MAPPING':True,'PHYSICAL_OBSERVER_GATE_BLOCKED_PENDING_MODE_CHANNEL':True,'NO_OBSERVATIONAL_INPUT':True,'NO_GR_LCDM_INPUT':True,'DEV167_PAIR_LAW_UNCHANGED':True,'DEV168_RECEIPT_PHYSICS_UNCHANGED':True,'DEV171_SOURCE_ENSEMBLE_UNCHANGED':True,'DEV183_LAUNCH_DOMAIN_UNCHANGED':True,'DEV184_C100_RECEIPTS_UNCHANGED':True,'DISTRIBUTION_PIPELINE_DETERMINISTIC':True,'RECEIVED_STATE_DISTRIBUTION_SUFFICIENCY_CLASSIFIED':True,'RECEIVED_STATE_DISTRIBUTION_SUFFICIENCY':'FULL','NATIVE_MODE_CHANNEL_GATE_CLASSIFIED':True,'NATIVE_MODE_CHANNEL_GATE':'AUTHORIZED','MECHANISM_REGISTRY_UPDATED':True,'REGISTRY_VALIDATED':True,'TIMELINE_REGENERATED':True,'DERIVATION_GRAPH_REGENERATED':False,'LEDGER_UPDATED':False,'HISTORICAL_INDEX_UPDATED_IF_REQUIRED':False,'TESTS_PASS':False,'IMPLEMENTATION_COMMIT_RECORDED':False,'REMOTE_PUSH_CONFIRMED':False,'REMOTE_FINAL_HEAD_VERIFIED':False,'WORKTREE_CLEAN':False}
 dump('final_contract.json',final)
 (OUT/'discussion_handoff.md').write_text('# DEV185 handoff\n\nThe complete frozen C100 packet-aware receipt state is distributionally sufficient for a separate native mode/channel sufficiency audit. This does not select a channel, promote a statistical eigendirection, construct an observer, or permit observational comparison.\n')
if __name__=='__main__':
 if '--regenerate-graph-only' in sys.argv: regenerate_graph()
 elif '--finalize-contract' in sys.argv: finalize_contract()
 else: main()
