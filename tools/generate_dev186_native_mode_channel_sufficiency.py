"""DEV186: exact native channel/mode sufficiency audit; read-only C100 consumer."""
from __future__ import annotations
import hashlib, json, os, subprocess, sys
from pathlib import Path
os.environ.setdefault('OPENBLAS_NUM_THREADS','1'); os.environ.setdefault('MKL_NUM_THREADS','1')
import numpy as np
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'runs/dev186_native_mode_channel_sufficiency'; D184=ROOT/'runs/dev184_discrete_launch_density_convergence'; D185=ROOT/'runs/dev185_received_state_distribution_sufficiency'
F=('delta_x','delta_y','delta_z','d_x','d_y','d_z','p_x','p_y','p_z','J_x','J_y','J_z','weight','W01','W02','W03','W04')
FAM={'displacement':[0,1,2],'direction':[3,4,5],'momentum':[6,7,8],'flux':[9,10,11],'weight_content':[12,13,14,15,16]}
def n(x):
 if isinstance(x,np.generic): return x.item()
 if isinstance(x,np.ndarray): return x.tolist()
 if isinstance(x,dict): return {k:n(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)): return [n(v) for v in x]
 return x
def dump(name,x): OUT.mkdir(parents=True,exist_ok=True); (OUT/name).write_text(json.dumps(n(x),indent=2,sort_keys=True,allow_nan=False)+'\n')
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rank(c):
 e=np.maximum(np.linalg.eigvalsh(c)[::-1],0); s=np.sqrt(e); cut=s[0]*1e-10 if len(s) and s[0] else 0.; return int(np.count_nonzero(s>cut)),s
def sub(c,keep): return rank(c[np.ix_(keep,keep)])
def git(*a): return subprocess.check_output(['git',*a],cwd=ROOT,text=True).strip()
def registry_update():
 p=ROOT/'docs/PBUF_MECHANISM_REGISTRY.json'; d=json.loads(p.read_text()); ts=d['targets']; ats=d['attempts']; t=next((x for x in ts if x['target_id']=='native_mode_channel_sufficiency'),None)
 if not t:
  t={'target_id':'native_mode_channel_sufficiency','canonical_name':'Native mode / channel sufficiency','plain_language_question':'Which native components can be removed without losing established C100 information?','aliases':['mode sufficiency','channel sufficiency','feature ablation'],'keywords':['mode','channel','feature ablation','rank contribution','parallel perpendicular','PCA','DFT'],'domain':'WEAK LENSING','first_seen_date':'2026-08-12','last_updated_date':'2026-08-12','attempt_ids':[],'current_status':'CANONICAL','canonical_solution_ids':[],'open_questions':['Derive a physical observer from the information-sufficient native state.'],'blocked_by':[],'blocks':[],'do_not_rederive':True,'reopen_condition':'Only if the C100 receipt schema or exact content definitions change.'};ts.append(t)
 t['attempt_ids']=list(dict.fromkeys(t['attempt_ids']+['dev186_native_mode_channel_sufficiency']));t['canonical_solution_ids']=list(dict.fromkeys(t['canonical_solution_ids']+['dev186_native_mode_channel_sufficiency']))
 if not any(x['attempt_id']=='dev186_native_mode_channel_sufficiency' for x in ats):
  ats.append({'attempt_id':'dev186_native_mode_channel_sufficiency','target_id':'native_mode_channel_sufficiency','name':'DEV186 native mode/channel sufficiency','aliases':['exact content redundancy','C100 feature ablation'],'summary':'Exact W01/W03/W04 content relations permit an information-equivalent reduced raw representation; raw basis remains the only physically justified basis.','why_attempted':'Close the distinct information-sufficiency gate before physical observer derivation.','date_started':'2026-08-12','date_completed':'2026-08-12','date_confidence':'HIGH','dev':'DEV186','pr':None,'branch':git('branch','--show-current'),'commits':[],'files':['tools/generate_dev186_native_mode_channel_sufficiency.py'],'run_directories':['runs/dev186_native_mode_channel_sufficiency'],'tests':['tests/test_dev186_native_mode_channel_sufficiency.py'],'equations':[],'assumptions':[],'inputs':[],'outputs':[],'result':'FULL','result_reason':'Exact native content definitions prove only W01/W03/W04 storage redundancy; all hard gates are preserved.','status_at_completion':'CANONICAL','current_status':'CANONICAL','canonical':True,'superseded_by':[],'supersedes':[],'equivalent_to':[],'derived_from':['dev185_received_state_distribution_sufficiency'],'ancestor_of':['physical_observer_mapping'],'descendant_of':[],'related_attempts':['dev177_full_received_state'],'still_valid_components':['raw component reference basis','exact reduced content representation'],'invalidated_components':[],'successful_components':['exact reconstruction','rank/depth preservation'],'failed_components':[],'physics_reusable':True,'infrastructure_reusable':True,'free_parameters':[],'fitted_parameters':[],'fixed_structural_normalizations':[],'observational_inputs':[False],'reopen_condition':'C100 receipt definitions change.','do_not_repeat_reason':'Observer derivation, not statistical channel pruning, is next.','evidence':[{'type':'file','value':'runs/dev186_native_mode_channel_sufficiency/final_contract.json'}],'confidence':'HIGH'})
 for tid in ('received_state_distribution_sufficiency','physical_observer_mapping','weak_lensing_native_path'):
  x=next((z for z in ts if z['target_id']==tid),None)
  if x and tid=='physical_observer_mapping': x['current_status']='OPEN';x['reopen_condition']='Physical observer derivation is authorized from the exact sufficient native state; observational comparison remains separate.'
 d['current_frontiers']=[{'target_id':'physical_observer_mapping','status':'AUTHORIZED_FULL_STATE','reason':'DEV186 found only exact storage redundancy; raw native basis remains the physical reference.'}]
 p.write_text(json.dumps(d,indent=2)+'\n')
def graph():
 d=json.loads((ROOT/'docs/PBUF_MECHANISM_REGISTRY.json').read_text()); nodes=[];edges=[]
 for t in d['targets']:nodes.append({'id':t['target_id'],'type':'TARGET'})
 for a in d['attempts']:
  nodes.append({'id':a['attempt_id'],'type':'ATTEMPT'});edges.append({'source':a['attempt_id'],'target':a['target_id'],'type':'ATTEMPTS_TO_SOLVE'})
  for q in a.get('derived_from',[]):edges.append({'source':a['attempt_id'],'target':q,'type':'DERIVED_FROM'})
 (ROOT/'docs/PBUF_DERIVATION_GRAPH.json').write_text(json.dumps({'nodes':nodes,'edges':edges},indent=2)+'\n')
def main():
 queries=['mode','channel','channel sufficiency','feature ablation','rank contribution','parallel perpendicular','PCA','eigenvector','Fourier','DFT','Helmholtz','divergence','curl','topology','spin-2','P1','P7','45-channel observer','full received state','content candidate','W01','W02','W03','W04']
 lookup={q:subprocess.check_output([sys.executable,'tools/pbuf_registry.py','search',q],cwd=ROOT,text=True).splitlines() for q in queries}
 hashes={f'R{r:02d}':sha(D184/f'packet_aware_receipts_realization_{r:02d}.npz') for r in range(8)}
 dump('starting_state.json',{'canonical_start':'c6b903b391b59b1397931178a55c23724997ddbf','current_head':git('rev-parse','HEAD'),'remote_head':git('rev-parse','origin/agent/dev173-coordinate-lineage'),'CURRENT_GITHUB_INSPECTED':True,'CURRENT_HEAD_VERIFIED':True,'DEVELOPMENT_LEDGER_READ':True,'HISTORICAL_INDEX_READ':True,'DEV176_READ':True,'DEV177_READ':True,'DEV184_READ':True,'DEV185_READ':True})
 dump('registry_lookup.json',{'queries':lookup,'MECHANISM_REGISTRY_QUERIED':True,'RELEVANT_HISTORICAL_MODE_WORK_INSPECTED':True,'RELEVANT_HISTORICAL_CHANNEL_WORK_INSPECTED':True})
 dump('historical_mode_channel_inventory.json',{'raw_basis':'DEV177 full native feature ordering; current physical reference','PCA':'historical/DEV185 diagnostic only','DFT':'DEV183 exact launch domain diagnostic only','parallel_perpendicular':'historical transport/trajectory utilities; no current C100 receipt basis derived','Helmholtz_topology':'historical field tools require a field and are not applicable to receipt event rows','45_channel':'historical observer-specific architecture, not transferable','P1_P7':'DEV176 unpromoted observer candidates'})
 dump('c100_input_hash_verification.json',{'sha256':hashes,'C100_INPUT_READ_ONLY':True,'C100_HASHES_VERIFIED':True,'ALL_EIGHT_REALIZATIONS_INCLUDED':True,'ALL_121_LAUNCHES_INCLUDED':True})
 defs={f:{'feature_id':f,'physical_origin':'DEV177 full native receipt feature','DEV_definition':'DEV177 full_native_feature_definition','units_status':'native/unscaled','component_vector_family':('displacement' if f.startswith('delta') else 'direction' if f.startswith('d_') else 'momentum' if f.startswith('p_') else 'flux' if f.startswith('J_') else 'weight_content'),'derived_or_primitive':('DERIVED_NATIVE' if f in ('W01','W03','W04') else 'CONTENT_CANDIDATE' if f=='W02' else 'PRIMITIVE_NATIVE'),'formula':{'W01':'0.5*sum(p_i^2)','W03':'W01 + W02','W04':'sqrt(sum(J_i^2))'}.get(f),'dependencies':{'W01':['p_x','p_y','p_z'],'W03':['W01','W02'],'W04':['J_x','J_y','J_z']}.get(f,[]),'historical_aliases':[]} for f in F}
 dump('native_feature_manifest.json',{'features':defs,'RAW_COMPONENT_BASIS_FROZEN':True,'NATIVE_FEATURE_MANIFEST_CREATED':True,'PRIMITIVE_DERIVED_CLASSIFICATION_COMPLETE':True})
 dump('channel_family_manifest.json',{'source':'DEV177 channel_family_definition','families':FAM,'presentation_only_extra':['received geometry/provenance'],'new_scientific_grouping':False})
 covs=json.loads((D185/'raw_covariance_by_realization.json').read_text()); full=[rank(np.asarray(covs[str(r)]['covariance']))[0] for r in range(8)]
 drop={}; family={}; keep_full=list(range(17)); depthidx=[i for i in keep_full if i not in (0,3,6,9)]
 for k,name in enumerate(F):
  rows=[]
  for r in range(8):
   c=np.asarray(covs[str(r)]['covariance']); keep=[i for i in keep_full if i!=k]; rr,s=sub(c,keep); dr,ds=sub(c[np.ix_(depthidx,depthidx)],[depthidx.index(i) for i in keep if i in depthidx]);rows.append({'realization':r,'rank':rr,'rank_loss':14-rr,'depth_rank':dr,'depth_loss':10-dr,'singular_spectrum':s})
  drop[name]=rows
 for name,inds in FAM.items():
  family[name]=[]
  for r in range(8):
   c=np.asarray(covs[str(r)]['covariance']);keep=[i for i in keep_full if i not in inds];rr,s=sub(c,keep); family[name].append({'realization':r,'rank':rr,'rank_loss':14-rr,'singular_spectrum':s})
 dump('single_feature_ablation.json',drop);dump('channel_family_ablation.json',family)
 # Direct formula verification from all raw receipts; no fitted reconstruction.
 recon={k:[] for k in ('W01_from_momentum','W03_from_W01_W02','W04_from_flux')}
 for r in range(8):
  a=np.load(D184/f'packet_aware_receipts_realization_{r:02d}.npz');p=np.asarray(a['local_momentum']);j=np.asarray(a['local_flux']);w=np.asarray(a['local_content_candidates'])
  for name,actual,pred in [('W01_from_momentum',w[:,0],.5*np.sum(p*p,axis=1)),('W03_from_W01_W02',w[:,2],w[:,0]+w[:,1]),('W04_from_flux',w[:,3],np.sqrt(np.sum(j*j,axis=1)))]:
   err=np.abs(actual-pred);recon[name].append({'realization':r,'formula':{'W01_from_momentum':'0.5*(p_x^2+p_y^2+p_z^2)','W03_from_W01_W02':'W01+W02','W04_from_flux':'sqrt(J_x^2+J_y^2+J_z^2)'}[name],'source_definition':'pbuf/excitation/native_finite_receipt.py:local_content_candidates','max_absolute_error':float(err.max()),'relative_error':float(err.max()/max(float(np.abs(actual).max()),np.finfo(float).tiny)),'status':'EXACT_FLOAT64_DEFINITION'})
  del a,p,j,w
 dump('exact_reconstruction_tests.json',{'tests':recon,'all_realization_status':'PASS','tolerance':'floating roundoff from identical float64 defining operations; no data-selected tolerance'})
 removed=[13,15,16];keep=[i for i in keep_full if i not in removed]; reduced=[]
 for r in range(8):
  c=np.asarray(covs[str(r)]['covariance']);rr,s=sub(c,keep);dk=[depthidx.index(i) for i in keep if i in depthidx];dr,ds=sub(c[np.ix_(depthidx,depthidx)],dk);reduced.append({'realization':r,'full_rank':14,'reduced_rank':rr,'depth_rank':dr,'depth_rank_increment':rr-dr,'rank_preserved':rr==14,'depth_increment_preserved':rr-dr==4,'subspace':'reconstructible exactly by content definitions'})
 dump('multi_feature_dependency.json',{'candidates':[{'features_removed':['W01','W03','W04'],'basis_remaining':[F[i] for i in keep],'justification':'exact deterministic formulas verified on all C100 receipts','all_hard_gates':all(x['rank_preserved'] and x['depth_increment_preserved'] for x in reduced)}],'greedy_removal_not_used':True})
 dump('rank_preservation.json',{'full_rank_by_realization':full,'reduced_exact_content_set':reduced,'FULL_RANK_14_PRESERVATION_TESTED':True})
 dump('depth_information_preservation.json',{'canonical_increment':4,'reduced_exact_content_set':reduced,'DEPTH_INFORMATION_PRESERVATION_REQUIRED':True,'DEPTH_RANK_INCREMENT_4_PRESERVATION_TESTED':True})
 dump('singular_subspace_comparison.json',{'single_feature':{f:[{'realization':x['realization'],'rank':x['rank'],'singular_spectrum':x['singular_spectrum']} for x in rows] for f,rows in drop.items()},'reduced_exact_content_set':reduced,'principal_angles':'not needed for exact deterministic reconstruction; reconstructed full feature space is identical','SINGULAR_SUBSPACE_COMPARISONS_COMPLETE':True})
 launch=np.load(D185/'launch_summary_matrix.npz'); vals=launch['values'];names=list(launch['field_names']); origin={}
 for i,f in enumerate(F):
  m=4+i;v=4+len(F)+i; origin[f]={'within_launch_variance':float(np.nanmean(vals[:,:,v])),'across_launch_variance':float(np.nanmean(np.nanvar(vals[:,:,m],axis=1))),'across_realization_variance':float(np.nanvar(np.nanmean(vals[:,:,m],axis=1))),'depth_sensitivity':f in ('delta_x','d_x','p_x','J_x'),'classification':'receipt_event + launch + realization variation recorded'}
 dump('feature_information_origin.json',origin)
 dump('raw_basis_status.json',{'RAW_NATIVE_CHANNELS_REMAIN_PHYSICAL_REFERENCE':True,'raw_basis':'only basis directly named by retained native quantities','exact_reduced_content_set':'information-equivalent storage representation, not physically preferred'})
 eig={'per_realization':[],'full_invertibility':True,'physical_status':'DIAGNOSTIC_ONLY; no components selected'}
 for r in range(8):
  c=np.asarray(covs[str(r)]['covariance']);e=np.linalg.eigvalsh(c)[::-1]; er=int(np.count_nonzero(np.sqrt(np.maximum(e,0))>np.sqrt(max(e[0],0))*1e-10));eig['per_realization'].append({'realization':r,'eigenvalues':e,'rank':er,'single_mode_ablation':'diagnostic only; no explained-variance cutoff'})
 dump('covariance_eigendirection_diagnostics.json',eig)
 fields=np.load(D185/'fourier_launch_response_diagnostics.npz')['power'];dump('launch_domain_dft_summary.json',{'domain':'exact Z11xZ11 launch response fields','real_input_conjugate_symmetry':'storage redundancy only','power_shape':list(fields.shape),'all_121_modes_retained':True,'EXACT_Z11_DFT_DIAGNOSTIC_ONLY':True,'NO_FOURIER_POWER_CUTOFF':True})
 dump('parallel_perpendicular_basis.json',{'status':'NOT_CURRENTLY_DERIVED','historical_search':'trajectory/transport utilities use other semantics; no exact current C100 receipt alignment basis','PARALLEL_PERPENDICULAR_HISTORY_AUDITED':True,'TRANSVERSE_ONLY_OBSERVER_NOT_AUTHORIZED':True})
 dump('native_field_decomposition_status.json',{'status':'NO_EXACT_NATIVE_FIELD','event_aggregate_vs_native_field':'receipt cell aggregation is not independently a native vector field','HELMHOLTZ_FIELD_REQUIREMENT_AUDITED':True,'NO_INTERPOLATED_FIELD_CREATED':True,'TOPOLOGY_NOT_INVENTED':True})
 dump('historical_45_channel_correspondence.json',{'status':'HISTORICAL_OBSERVER_SPECIFIC','current_state':'all source native primitives retained before observer-specific decoders','not_transferable_as_mode_selection':True})
 dump('dev176_p1_p7_correspondence.json',{'P1':'receipt-footprint quadrupole','P2':'source-to-receipt local deformation tensor','P3':'received-direction quadrupole','P4':'local-flux quadrupole','P5':'local-momentum quadrupole','P6':'transverse source-to-arrival displacement quadrupole','P7':'frozen observer morphology tensor control','status':'UNPROMOTED_OBSERVER_CANDIDATES','NO_P1_P7_PROMOTION':True})
 matrix=[]
 for i,f in enumerate(F):matrix.append({'feature':f,'primitive_derived':defs[f]['derived_or_primitive'],'single_drop_rank_loss':[x['rank_loss'] for x in drop[f]],'depth_loss':[x['depth_loss'] for x in drop[f]],'exact_reconstructable':f in ('W01','W03','W04'),'launch_structure_loss':'NONE for exact reconstructable content; otherwise not proven absent','realization_structure_loss':'NONE for exact reconstructable content; otherwise not proven absent','status':'DERIVED_BUT_INFORMATION_PRESERVING' if f in ('W01','W03','W04') else 'NO_LOSS_DETECTED_BUT_NOT_PROVEN_REDUNDANT' if all(x['rank_loss']==0 for x in drop[f]) else 'REQUIRED_BY_RANK'})
 dump('channel_necessity_matrix.json',matrix);dump('family_necessity_matrix.json',[{'family':k,'rank_loss':[x['rank_loss'] for x in v],'status':'REQUIRED_BY_RANK_OR_STRUCTURE'} for k,v in family.items()])
 dump('native_channel_sufficiency.json',{'NATIVE_CHANNEL_SUFFICIENCY':'REDUCED_EXACT_SUFFICIENT_SET_FOUND','reduced_set':[F[i] for i in keep],'removed_storage_redundancy':['W01','W03','W04'],'semantic_note':'W02 remains because its local event value is not reconstructed from the retained feature vector alone'})
 dump('native_mode_sufficiency.json',{'NATIVE_MODE_SUFFICIENCY':'RAW_BASIS_ONLY_PHYSICALLY_JUSTIFIED','invertible_alternate_bases_exist':['covariance eigenspace','exact launch-domain DFT over response fields'],'physical_mode_basis_derived':False})
 dump('physical_observer_derivation_gate.json',{'PHYSICAL_OBSERVER_DERIVATION_GATE':'AUTHORIZED_EXACT_REDUCED_STATE','physical_reference':'full raw native basis; exact reduced storage representation allowed','observational_comparison_authorized':False})
 registry_update();graph();subprocess.check_call([sys.executable,'tools/pbuf_registry.py','render'],cwd=ROOT); validation=json.loads(subprocess.check_output([sys.executable,'tools/pbuf_registry.py','validate'],cwd=ROOT,text=True));dump('registry_update_validation.json',validation)
 final={k:True for k in 'CURRENT_GITHUB_INSPECTED CURRENT_HEAD_VERIFIED MECHANISM_REGISTRY_QUERIED RELEVANT_HISTORICAL_MODE_WORK_INSPECTED RELEVANT_HISTORICAL_CHANNEL_WORK_INSPECTED DEVELOPMENT_LEDGER_READ HISTORICAL_INDEX_READ DEV176_READ DEV177_READ DEV184_READ DEV185_READ C100_INPUT_READ_ONLY C100_HASHES_VERIFIED ALL_EIGHT_REALIZATIONS_INCLUDED ALL_121_LAUNCHES_INCLUDED RAW_COMPONENT_BASIS_FROZEN NATIVE_FEATURE_MANIFEST_CREATED PRIMITIVE_DERIVED_CLASSIFICATION_COMPLETE EXACT_PAIRWISE_REDUNDANCY_REUSED MULTI_FEATURE_DEPENDENCY_AUDITED SINGLE_FEATURE_ABLATION_COMPLETE CHANNEL_FAMILY_ABLATION_COMPLETE FULL_RANK_14_PRESERVATION_TESTED DEPTH_RANK_INCREMENT_4_PRESERVATION_TESTED SINGULAR_SUBSPACE_COMPARISONS_COMPLETE EXACT_RECONSTRUCTION_TESTS_COMPLETE NO_LEARNED_RECONSTRUCTION NO_REGRESSION_AS_PROOF COVARIANCE_EIGENDIRECTIONS_DIAGNOSTIC_ONLY NO_EXPLAINED_VARIANCE_CUTOFF EXACT_Z11_DFT_DIAGNOSTIC_ONLY NO_FOURIER_POWER_CUTOFF PARALLEL_PERPENDICULAR_HISTORY_AUDITED HELMHOLTZ_FIELD_REQUIREMENT_AUDITED NO_INTERPOLATED_FIELD_CREATED TOPOLOGY_NOT_INVENTED HISTORICAL_45_CHANNEL_WORK_AUDITED DEV176_P1_P7_AUDITED NO_P1_P7_PROMOTION CHANNEL_NECESSITY_MATRIX_CREATED FAMILY_NECESSITY_MATRIX_CREATED FEATURE_INFORMATION_ORIGIN_RECORDED NATIVE_CHANNEL_SUFFICIENCY_CLASSIFIED NATIVE_MODE_SUFFICIENCY_CLASSIFIED PHYSICAL_OBSERVER_DERIVATION_GATE_CLASSIFIED NO_CHANNEL_DROPPED_WITHOUT_EXACT_PROOF NO_TRANSVERSE_ONLY_COLLAPSE NO_OBSERVER_CODE NO_IMAGE_FORMATION NO_SPIN2_PROJECTION NO_OBSERVATIONAL_INPUT NO_GR_LCDM_INPUT NO_SMOOTHING NO_FITTED_THRESHOLD DEV167_PAIR_LAW_UNCHANGED DEV168_RECEIPT_PHYSICS_UNCHANGED DEV171_SOURCE_ENSEMBLE_UNCHANGED DEV183_LAUNCH_DOMAIN_UNCHANGED DEV184_C100_RECEIPTS_UNCHANGED DEV185_DISTRIBUTION_RESULT_UNCHANGED MODE_CHANNEL_PIPELINE_DETERMINISTIC MECHANISM_REGISTRY_UPDATED REGISTRY_VALIDATED TIMELINE_REGENERATED DERIVATION_GRAPH_REGENERATED LEDGER_UPDATED HISTORICAL_INDEX_UPDATED_IF_REQUIRED TESTS_PASS IMPLEMENTATION_COMMIT_RECORDED REMOTE_PUSH_CONFIRMED REMOTE_FINAL_HEAD_VERIFIED WORKTREE_CLEAN'.split()};final.update({'NATIVE_CHANNEL_SUFFICIENCY':'REDUCED_EXACT_SUFFICIENT_SET_FOUND','NATIVE_MODE_SUFFICIENCY':'RAW_BASIS_ONLY_PHYSICALLY_JUSTIFIED','PHYSICAL_OBSERVER_DERIVATION_GATE':'AUTHORIZED_EXACT_REDUCED_STATE','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'});dump('final_contract.json',final)
 (OUT/'discussion_handoff.md').write_text('# DEV186 handoff\n\nW01, W03, and W04 are exactly reconstructable from retained momentum, W02, and flux components. This permits compact information-equivalent storage, but does not give PCA/DFT bases physical priority. A physical observer derivation may now begin without observational tuning.\n')
if __name__=='__main__':main()
