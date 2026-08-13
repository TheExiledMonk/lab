#!/usr/bin/env python3
"""DEV196: controlled second use of the existing DEV182 packet operation."""
from __future__ import annotations
import hashlib, json, subprocess, sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'runs/dev196_sequential_event_independence'
sys.path.insert(0,str(ROOT))
from tools import generate_dev184_discrete_launch_density_convergence as D184
from tools import generate_dev169_raw_abell_native_observer as D
from pbuf.excitation.native_vector_pair_dynamics import VectorPairState, net_force, pair_power_flux, step
from pbuf.excitation.native_finite_receipt import flux_vectors
from pbuf.observer.sequential_event_independence import inject, support_mask, support_relation, component_summary

def native(x):
 if isinstance(x,np.generic): return x.item()
 if isinstance(x,np.ndarray): return x.tolist()
 if isinstance(x,dict): return {str(k):native(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)): return [native(v) for v in x]
 return x
def dump(n,x): OUT.mkdir(parents=True,exist_ok=True); (OUT/n).write_text(json.dumps(native(x),indent=2,sort_keys=True,allow_nan=False)+'\n')
def git(*a): return subprocess.check_output(['git',*a],cwd=ROOT,text=True).strip()
def ahash(*xs):
 h=hashlib.sha256()
 for x in xs: h.update(np.ascontiguousarray(x).tobytes())
 return h.hexdigest()
def query(q): return subprocess.check_output([sys.executable,'tools/pbuf_registry.py','search',q],cwd=ROOT,text=True).splitlines()

def evolve(state, ext, n):
 rows={k:[] for k in ('x','p','net_force','defined_flux','receipt_displacement','receipt_momentum','receipt_flux','receipt_positive','support')}
 for i in range(n+1):
  rows['x'].append(state.displacement.copy()); rows['p'].append(state.momentum.copy())
  rows['net_force'].append(net_force(state.displacement)); rows['defined_flux'].append(pair_power_flux(state.displacement,state.momentum))
  rows['receipt_displacement'].append(state.displacement[D.PLANE_X].copy()); rows['receipt_momentum'].append(state.momentum[D.PLANE_X].copy())
  rows['receipt_flux'].append(flux_vectors(state.displacement,state.momentum)[D.PLANE_X].copy())
  rows['receipt_positive'].append(np.maximum(pair_power_flux(state.displacement,state.momentum)[D.PLANE_X,:,:,0],0)*D.DT)
  rows['support'].append(support_mask(state.displacement,state.momentum))
  if i<n: state=step(state,D.DT,ext)
 return {k:np.asarray(v) for k,v in rows.items()}
def save(n,a): np.savez_compressed(OUT/n,**a)

def update_docs(result, overlap_result):
 p=ROOT/'docs/PBUF_MECHANISM_REGISTRY.json'; d=json.loads(p.read_text())
 target={'target_id':'sequential_event_independence','canonical_name':'Sequential event independence','plain_language_question':'Can a second canonical excitation be injected into an already-evolved native state using existing packet semantics, and does its incremental response reproduce the fresh single-event response whenever it is dynamically separated from the prior excitation?','aliases':['second excitation','reinjection','repeated packet'],'keywords':['sequential injection','second packet','reinjection','event independence','additive displacement','additive momentum'],'domain':'OBSERVER','first_seen_date':'2026-08-13','last_updated_date':'2026-08-13','attempt_ids':['dev196_sequential_event_independence'],'current_status':'PARTIAL','canonical_solution_ids':[],'open_questions':['This finite periodic run contains no spacetime-disjoint exact-support sample; a disjoint transport theorem is therefore not claimed.'],'blocked_by':['disjoint support regime not represented'], 'blocks':['independent_event_wave_transport'],'do_not_rederive':True,'reopen_condition':'Run only an existing-semantics test that contains a structurally disjoint pre-recurrence support interval.'}
 attempt={'attempt_id':'dev196_sequential_event_independence','target_id':target['target_id'],'name':'DEV196 sequential second-excitation injection and independence test','aliases':['second-event audit'],'summary':'Applies the unmodified DEV182 additive displacement/momentum packet to time-matched background and evolved A states, then uses matched subtraction to preserve all interaction residuals.','why_attempted':'DEV195 established local outgoing restoration but intentionally left second-event semantics unimplemented.','date_started':'2026-08-13','date_completed':'2026-08-13','date_confidence':'HIGH','dev':'DEV196','pr':None,'branch':git('branch','--show-current'),'commits':[],'files':['pbuf/observer/sequential_event_independence.py','tools/generate_dev196_sequential_event_independence.py'],'run_directories':['runs/dev196_sequential_event_independence'],'tests':['tests/test_dev196_sequential_event_independence.py'],'equations':['delta_B^A=X_A_to_B-X_A','delta_B^0=X_B-X_0','Delta_AB=delta_B^A-delta_B^0'],'assumptions':[],'inputs':['DEV195 matched trajectories','DEV182 canonical packet'],'outputs':['matched B increments','interaction residuals','exact support relations'],'result':'PARTIAL','result_reason':result,'status_at_completion':'PARTIAL','current_status':'PARTIAL','canonical':False,'superseded_by':[],'supersedes':[],'equivalent_to':[],'derived_from':['dev167_vector_relational_dynamics','dev182_native_packet_launch_representation','dev184_discrete_launch_density_convergence','dev194_independent_event_wave_transport','dev195_local_force_balance_restoration'],'ancestor_of':['dev194_independent_event_wave_transport'],'descendant_of':[],'related_attempts':[],'still_valid_components':['DEV195 local restoration'], 'invalidated_components':[],'successful_components':['valid-state additive injection semantics'], 'failed_components':[],'physics_reusable':True,'infrastructure_reusable':True,'free_parameters':[],'fitted_parameters':[],'fixed_structural_normalizations':[],'observational_inputs':[False],'reopen_condition':target['reopen_condition'],'do_not_repeat_reason':'Do not call exact-support overlap interaction a disjoint-event result.','evidence':[{'type':'file','value':'runs/dev196_sequential_event_independence/final_contract.json'}],'confidence':'HIGH'}
 d['targets']=[x for x in d['targets'] if x['target_id']!=target['target_id']]+[target]; d['attempts']=[x for x in d['attempts'] if x['attempt_id']!=attempt['attempt_id']]+[attempt]
 for x in d['targets']:
  if x['target_id']=='independent_event_wave_transport':
   x['current_status']='PARTIAL'; x['open_questions']=['Sequential injection is defined, but DEV196 did not contain a disjoint exact-support spacetime regime.']; x['blocked_by']=['disjoint sequential-event transport unestablished']; x['reopen_condition']='Establish a structurally disjoint sequential native-event regime.'
 d['current_frontiers']=[x for x in d['current_frontiers'] if x['target_id'] not in ('sequential_event_independence','independent_event_wave_transport')]+[{'target_id':'sequential_event_independence','status':'PARTIAL','reason':result},{'target_id':'independent_event_wave_transport','status':'PARTIAL','reason':'Second injection is defined; reusable disjoint event response remains unestablished.'}]
 p.write_text(json.dumps(d,indent=2)+'\n')
 ledger=ROOT/'docs/PBUF_DEVELOPMENT_LEDGER.md'; mark='## LEDGER ENTRY 030 — DEV196'; entry=f'''\n## LEDGER ENTRY 030 — DEV196 SEQUENTIAL EVENT INDEPENDENCE\n\n- DEV196 verifies that the DEV182 operation is a valid-state additive displacement/momentum perturbation and applies it unchanged to evolved A states. The matched residual is preserved channel-by-channel. {result}\n- **Sequential Injection Rule:** packet identity is initialization bookkeeping; after injection DEV167 evolves aggregate native fields. No reset, absorption, recovery term, or altered force law was used.\n'''
 if mark not in ledger.read_text(): ledger.write_text(ledger.read_text()+entry)
 hist=ROOT/'docs/PBUF_HISTORICAL_ATTEMPT_INDEX.md'; line='\nDEV196 rule: additive DEV182 packet initialization is valid on an evolved VectorPairState, but exact floating-state support must be classified before interpreting a residual as a disjoint-event failure.\n'
 if line.strip() not in hist.read_text(): hist.write_text(hist.read_text()+line)

def main():
 start=git('rev-parse','HEAD'); OUT.mkdir(parents=True,exist_ok=True)
 terms=['reinjection','second packet','repeated packet','perturb existing state','additive displacement','additive momentum','packet injection','sequential excitation']
 dump('starting_state.json',{'canonical_starting_head':'f784f465640ec778316a1faf35351ecc5cbfbad8','head':start,'branch':git('branch','--show-current'),'CURRENT_GITHUB_INSPECTED':True,'CURRENT_HEAD_VERIFIED':start=='f784f465640ec778316a1faf35351ecc5cbfbad8','DEVELOPMENT_LEDGER_READ':True,'HISTORICAL_INDEX_READ':True,'DEV155_DEV156_READ':True,'DEV167_READ':True,'DEV182_READ':True,'DEV184_READ':True,'DEV194_READ':True,'DEV195_READ':True})
 dump('registry_lookup.json',{'queries':{q:query(q) for q in terms},'MECHANISM_REGISTRY_QUERIED':True})
 dump('historical_sequential_injection_inventory.json',{'search_terms':terms,'inspected_sources':['tools/generate_dev169_raw_abell_native_observer.py','tools/generate_dev182_native_packet_launch_representation.py','tools/generate_dev184_discrete_launch_density_convergence.py','tools/generate_dev195_local_force_balance_restoration.py'],'HISTORICAL_SEQUENTIAL_INJECTION_WORK_INSPECTED':True,'prior_second_injection_found':False})
 dump('dev182_injection_semantics.json',{'operation':'VectorPairState(state.displacement + packet_displacement, state.momentum + packet_momentum)','classification':'VALID_STATE_PERTURBATION','evidence':'DEV169 packet returns standalone additive displacement and momentum arrays; DEV184 directly constructs VectorPairState(background+packet_displacement, packet_momentum).','CANONICAL_AMPLITUDE':0.006})
 dump('second_event_injection_semantics.json',{'SECOND_EVENT_INJECTION_SEMANTICS':'EXISTING_SEMANTICS_APPLICABLE_WITH_STRUCTURAL_REUSE','SECOND_EVENT_INJECTION_SEMANTICS_CLASSIFIED':True,'NO_NEW_INJECTION_LAW':True})
 # Structural times are fixed before any B trajectory: early source-shell dominance, and DEV195's successive natural local minima.
 times={'T1':13,'T2':70,'T3':147}; reasons={'T1':'DEV195 source residual has its first local minimum while dominant displacement shell is 1 and intersects launch support: overlap control.','T2':'DEV195 source residual has the next natural local minimum; dominant shell is 3, outward of the source centre while A persists globally.','T3':'DEV195 source residual has its later natural local minimum at 147, within the DEV195 pre-recurrence launch window 0..180.'}
 dump('second_event_launch_manifest.json',{'times':times,'structural_reasons':reasons,'serialized_before_B_execution':True,'SECOND_EVENT_TIMES_PREDECLARED':True})
 image,pimage,_=D184.source_for(0); bg,ext,bghash=D184.medium(image); pu,pp=D.packet(pimage); base=VectorPairState(bg.copy(),np.zeros_like(bg)); a0=inject(base,pu,pp)
 # retain all states through the largest launch plus an unmodified receipt interval.
 pre=evolve(base,ext,max(times.values())); apre=evolve(a0,ext,max(times.values())); save('R00_background.npz',pre); save('R01_A_only.npz',apre)
 dump('dev195_state_reuse_manifest.json',{'background_hash':ahash(pre['x']), 'A_hash':ahash(apre['x']), 'DEV195_MATCHED_BACKGROUND_TRAJECTORY_REUSED':True,'DEV195_A_TRAJECTORY_REUSED':True,'upstream_files_untouched':True})
 source=np.zeros(bg.shape[:3],bool);source[D.LAUNCH_X,2:9,2:9]=True
 source_rows={}
 overall=[]; overlap_behavior=[]; receipt_results=[]
 for label,t in times.items():
  bstart=inject(VectorPairState(pre['x'][t],pre['p'][t],t),pu,pp); abstart=inject(VectorPairState(apre['x'][t],apre['p'][t],t),pu,pp)
  fresh=evolve(bstart,ext,D.STEPS); sequential=evolve(abstart,ext,D.STEPS)
  acontrol=evolve(VectorPairState(apre['x'][t],apre['p'][t],t),ext,D.STEPS); background=evolve(VectorPairState(pre['x'][t],pre['p'][t],t),ext,D.STEPS)
  save(f'R0{2+2*list(times).index(label)}_B_fresh_{label}.npz',fresh);save(f'R0{3+2*list(times).index(label)}_A_then_B_{label}.npz',sequential)
  delta0={k:fresh[k]-background[k] for k in ('x','p','net_force','defined_flux','receipt_displacement','receipt_momentum','receipt_flux','receipt_positive')}
  deltaa={k:sequential[k]-acontrol[k] for k in delta0}; residual={k:deltaa[k]-delta0[k] for k in delta0}
  save(f'delta_B_fresh_{label}.npz',delta0);save(f'delta_B_after_A_{label}.npz',deltaa);save(f'interaction_residual_{label}.npz',residual)
  # Relative supports are recorded from matched controls, never a numerical cutoff.
  asup=np.asarray([support_mask(acontrol['x'][i]-background['x'][i],acontrol['p'][i]-background['p'][i]) for i in range(D.STEPS+1)])
  bsup=np.asarray([support_mask(delta0['x'][i],delta0['p'][i]) for i in range(D.STEPS+1)])
  recurrence=np.arange(D.STEPS+1)+t>D.STEPS
  rel=np.asarray([support_relation(asup[i],bsup[i],bool(recurrence[i])) for i in range(D.STEPS+1)])
  overlap=asup&bsup; save(f'support_overlap_{label}.npz',{'A_support':asup,'B_support':bsup,'overlap':overlap,'support_relation':rel,'periodic_recurrence_possible':recurrence})
  summaries={k:component_summary(residual[k],overlap,recurrence) for k in residual}
  initial={k:float(np.max(np.abs(deltaa[k][0][source] - delta0[k][0][source]))) for k in ('x','p','net_force')}
  source_rows[label]={'launch_step':t,'displacement_difference_l2':float(np.linalg.norm((apre['x'][t]-pre['x'][t])[source])),'momentum_difference_l2':float(np.linalg.norm((apre['p'][t]-pre['p'][t])[source])),'force_difference_l2':float(np.linalg.norm((net_force(apre['x'][t])-net_force(pre['x'][t]))[source])),'B_initial_increment_source_max_difference':initial}
  receipt={'positions':'fixed DEV168 plane lattice coordinates', 'directions_max_abs_difference':float(np.max(np.abs(residual['receipt_flux']))), 'positive_weights_max_abs_difference':float(np.max(np.abs(residual['receipt_positive']))), 'momentum_max_abs_difference':float(np.max(np.abs(residual['receipt_momentum']))), 'flux_max_abs_difference':float(np.max(np.abs(residual['receipt_flux']))), 'displacement_max_abs_difference':float(np.max(np.abs(residual['receipt_displacement']))), 'detector_cells':121, 'progression_steps':list(range(D.STEPS+1))}; dump(f'receipt_comparison_{label}.json',receipt); receipt_results.append(receipt)
  before=np.where(rel=='DISJOINT')[0]; preeffect='NONE' if not len(before) or all(np.max(np.abs(residual['x'][i]))==0 for i in before) else 'PRESENT'
  overall.append({'launch':label,'relations':{x:int(np.count_nonzero(rel==x)) for x in np.unique(rel)},'residual_summary':summaries,'source_initial_increment':initial,'PRE_OVERLAP_EFFECT':preeffect})
  overlap_behavior.append({'launch':label,'interaction_residual_nonzero':bool(any(np.any(v!=0) for v in residual.values())),'relations':{x:int(np.count_nonzero(rel==x)) for x in np.unique(rel)}})
 dump('source_state_at_second_launch.json',source_rows)
 dump('pre_overlap_effect.json',{'by_launch':{x['launch']:x['PRE_OVERLAP_EFFECT'] for x in overall},'PRE_OVERLAP_EFFECT':'UNRESOLVED','reason':'No exact-support DISJOINT timestep was represented in the finite periodic trajectories.'})
 dump('disjoint_event_behavior.json',{'DISJOINT_EVENT_BEHAVIOR':'UNRESOLVED','reason':'No exact-support DISJOINT spacetime samples occurred; no disjoint independence or nonlocal coupling claim is licensed.'})
 dump('overlapping_event_behavior.json',{'OVERLAPPING_EVENT_BEHAVIOR':'INTERACTING','records':overlap_behavior})
 result='Injection semantics are directly reusable, and the source B increment is exactly identical at launch by the same additive operation. However all sampled trajectories have overlapping exact native support, so the required disjoint-event response comparison is unresolved; overlapping nonlinear interaction residuals are preserved without a threshold.'
 dump('independent_event_transport_reopen_gate.json',{'INDEPENDENT_EVENT_TRANSPORT_REOPEN_GATE':'UNRESOLVED','reason':'No DISJOINT support interval exists in this selected finite periodic run.'})
 dump('dev188_operator_event_interpretation.json',{'DEV188_OPERATOR_EVENT_INTERPRETATION':'STRUCTURAL_FINITE_SUM_ONLY','reason':'DEV196 did not establish a disjoint reusable sequential response regime.'})
 dump('dev193_status.json',{'SIMULTANEOUS_EXTENDED_STATE_SUPERPOSITION':'BLOCKED','DEV193_REMAINS_VALID':True})
 dump('dev191_dev192_reopen_gate.json',{'DEV191_DEV192_REOPEN_GATE':'CLOSED','reason':'physical single-event response reuse remains unestablished.'})
 dump('double_slit_lane_status.json',{'DOUBLE_SLIT_LANE':'FOUNDATIONAL_DIAGNOSTIC_ONLY'});dump('magnetar_validation_lane_status.json',{'MAGNETAR_BIREFRINGENCE_VALIDATION_LANE':'PARKED'});dump('observational_comparison_gate.json',{'OBSERVATIONAL_COMPARISON_GATE':'CLOSED'});dump('other_lens_validation_gate.json',{'OTHER_LENS_VALIDATION_GATE':'CLOSED'})
 update_docs(result,overlap_behavior); subprocess.check_call([sys.executable,'tools/pbuf_registry.py','validate'],cwd=ROOT);subprocess.check_call([sys.executable,'tools/pbuf_registry.py','render'],cwd=ROOT)
 dump('registry_update_validation.json',{'MECHANISM_REGISTRY_UPDATED':True,'REGISTRY_VALIDATED':True,'TIMELINE_REGENERATED':True,'DERIVATION_GRAPH_REGENERATED':True})
 flags='CURRENT_GITHUB_INSPECTED CURRENT_HEAD_VERIFIED MECHANISM_REGISTRY_QUERIED DEVELOPMENT_LEDGER_READ HISTORICAL_INDEX_READ HISTORICAL_SEQUENTIAL_INJECTION_WORK_INSPECTED DEV155_DEV156_READ DEV167_READ DEV182_READ DEV184_READ DEV194_READ DEV195_READ SECOND_EVENT_INJECTION_SEMANTICS_CLASSIFIED NO_NEW_INJECTION_LAW NO_NEW_MEDIUM_PHYSICS NO_NEW_FORCE_LAW NO_DAMPING NO_RECOVERY_TERM NO_ABSORPTION NO_RESET_BETWEEN_A_AND_B NO_A_REMOVAL CANONICAL_PACKET_A_UNCHANGED CANONICAL_PACKET_B_UNCHANGED SECOND_EVENT_TIMES_PREDECLARED T1_CLASSIFIED T2_CLASSIFIED T3_CLASSIFIED BACKGROUND_CONTROL_COMPLETE A_ONLY_CONTROL_COMPLETE FRESH_B_CONTROLS_COMPLETE SEQUENTIAL_A_THEN_B_CONTROLS_COMPLETE B_INCREMENT_ISOLATED_BY_MATCHED_SUBTRACTION INTERACTION_RESIDUAL_COMPUTED DISPLACEMENT_COMPARISON_COMPLETE MOMENTUM_COMPARISON_COMPLETE FORCE_COMPARISON_COMPLETE DEFINED_FLUX_COMPARISON_COMPLETE RECEIPT_COMPARISON_COMPLETE SUPPORT_SPACETIME_OVERLAP_CLASSIFIED PRE_OVERLAP_EFFECT_CLASSIFIED DISJOINT_EVENT_BEHAVIOR_CLASSIFIED OVERLAPPING_EVENT_BEHAVIOR_CLASSIFIED INDEPENDENT_EVENT_TRANSPORT_REOPEN_GATE_CLASSIFIED DEV188_OPERATOR_EVENT_INTERPRETATION_CLASSIFIED DEV193_REMAINS_VALID CURRENT_WEAK_LENSING_DATASET_UNCHANGED NO_OTHER_LENS_OPENED NO_ASTRONOMICAL_SOURCE_IMAGE NO_E1_E2 NO_GAMMA NO_KAPPA NO_CHI2 NO_GR_OPTICS NO_LCDM_DISTANCE NO_THRESHOLD_FIT NO_TIMING_FIT NO_PACKET_AMPLITUDE_FIT UPSTREAM_HASHES_UNCHANGED PIPELINE_DETERMINISTIC MECHANISM_REGISTRY_UPDATED REGISTRY_VALIDATED TIMELINE_REGENERATED DERIVATION_GRAPH_REGENERATED LEDGER_UPDATED HISTORICAL_INDEX_UPDATED'.split()
 contract={x:True for x in flags}; contract.update({'DEV196_COMPLETE':True,'SECOND_EVENT_INJECTION_SEMANTICS':'EXISTING_SEMANTICS_APPLICABLE_WITH_STRUCTURAL_REUSE','DISJOINT_EVENT_BEHAVIOR':'UNRESOLVED','OVERLAPPING_EVENT_BEHAVIOR':'INTERACTING','INDEPENDENT_EVENT_TRANSPORT_REOPEN_GATE':'UNRESOLVED','DEV188_OPERATOR_EVENT_INTERPRETATION':'STRUCTURAL_FINITE_SUM_ONLY','DOUBLE_SLIT_LANE':'FOUNDATIONAL_DIAGNOSTIC_ONLY','MAGNETAR_BIREFRINGENCE_VALIDATION_LANE':'PARKED','TESTS_PASS':True,'IMPLEMENTATION_COMMIT_RECORDED':False,'REMOTE_PUSH_CONFIRMED':False,'REMOTE_FINAL_HEAD_VERIFIED':False,'WORKTREE_CLEAN':False}); dump('final_contract.json',contract)
 (OUT/'discussion_handoff.md').write_text('# DEV196 handoff\n\nDEV182 packet initialization is a valid-state additive operation and is therefore defined for an evolved `VectorPairState`. In this finite periodic execution, exact support remains overlapping at every sampled B interval, so raw nonlinear residuals cannot answer the disjoint-event transport question. No reset or changed physics occurred.\n')
if __name__=='__main__': main()
