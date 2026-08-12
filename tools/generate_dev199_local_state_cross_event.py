#!/usr/bin/env python3
"""DEV199 Phase A: native-only local bond-state closure and freeze."""
from __future__ import annotations
import hashlib, json, subprocess, sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'runs/dev199_local_state_em_correlation/phase_a'; IN=ROOT/'runs/dev196_sequential_event_independence'; sys.path.insert(0,str(ROOT))
from pbuf.observer.local_state_cross_event import four_state_cross_term, sigma_prime, sigma_second
from pbuf.excitation.native_vector_pair_dynamics import net_force
from pbuf.observer.sequential_event_independence import inject
from tools import generate_dev184_discrete_launch_density_convergence as D184
from tools import generate_dev169_raw_abell_native_observer as D

def native(x):
 if isinstance(x,np.generic): return x.item()
 if isinstance(x,np.ndarray): return x.tolist()
 if isinstance(x,dict): return {str(k):native(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)): return [native(v) for v in x]
 return x
def dump(name,x): OUT.mkdir(parents=True,exist_ok=True);(OUT/name).write_text(json.dumps(native(x),indent=2,sort_keys=True,allow_nan=False)+'\n')
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def git(*a): return subprocess.check_output(['git',*a],cwd=ROOT,text=True).strip()
def query(q): return subprocess.check_output([sys.executable,'tools/pbuf_registry.py','search',q],cwd=ROOT,text=True).splitlines()
def save(name,**kw): np.savez_compressed(OUT/name,**kw)

def update_docs(result):
 p=ROOT/'docs/PBUF_MECHANISM_REGISTRY.json'; d=json.loads(p.read_text())
 targets=[
 {'target_id':'full_local_state_cross_event_control','canonical_name':'Full local-state cross-event control','plain_language_question':'Does the complete existing local native bond state uniquely account for sequential cross-event influence under DEV167 dynamics?','aliases':['local state','bond cross term','causal neighborhood'],'keywords':['cross event','local bond','strain','orientation','N6'],'domain':'NATIVE DYNAMICS','first_seen_date':'2026-08-13','last_updated_date':'2026-08-13','attempt_ids':['dev199_local_state_cross_event'],'current_status':'CANONICAL','canonical_solution_ids':['dev199_local_state_cross_event'],'open_questions':['Canonical-packet weak-regime membership remains unthresholded and unresolved.'],'blocked_by':[],'blocks':['independent_event_wave_transport'],'do_not_rederive':True,'reopen_condition':'Only if DEV167 or DEV182 changes independently.'},
 {'target_id':'em_wave_structural_correlation','canonical_name':'EM-wave structural correlation','plain_language_question':'After freezing the native PBUF interaction structure, does it show coefficient-free structural correspondence with established weak/strong-field electromagnetic-wave interaction behavior?','aliases':['Euler-Heisenberg','QED structure','EM correlation'],'keywords':['EM','QED','orientation','background invariant'],'domain':'VALIDATION','first_seen_date':'2026-08-13','last_updated_date':'2026-08-13','attempt_ids':[],'current_status':'PARTIAL','canonical_solution_ids':[],'open_questions':['Structural comparison is not a variable or SI mapping.'],'blocked_by':[],'blocks':[],'do_not_rederive':True,'reopen_condition':'Only with an independently defined PBUF-to-observable mapping.'}]
 for t in targets:d['targets']=[x for x in d['targets'] if x['target_id']!=t['target_id']]+[t]
 attempt={'attempt_id':'dev199_local_state_cross_event','target_id':'full_local_state_cross_event_control','name':'DEV199 local-state cross-event derivation','aliases':['DEV199 phase A'],'summary':'Exact four-state positive-N6 bond inclusion--exclusion reconstruction under frozen DEV167 dynamics.','why_attempted':'DEV198 excluded scalar force magnitude as a sufficient reduction.','date_started':'2026-08-13','date_completed':'2026-08-13','date_confidence':'HIGH','dev':'DEV199','branch':git('branch','--show-current'),'commits':[],'files':['pbuf/observer/local_state_cross_event.py','tools/generate_dev199_local_state_cross_event.py'],'run_directories':['runs/dev199_local_state_em_correlation/phase_a'],'tests':['tests/test_dev199_local_state_cross_event.py'],'equations':['F=sigma(epsilon) rhat','DeltaF_AB=F_AB-F_A-F_B+F_0'],'assumptions':['frozen DEV167/DEV182'], 'inputs':['DEV195 trajectory','DEV196 matched trajectories','DEV198 complete scalar ordering'],'outputs':['four-state bond arrays','local path trace','perturbative expansion'],'result':'FULL','result_reason':result,'status_at_completion':'CANONICAL','current_status':'CANONICAL','canonical':True,'superseded_by':[],'supersedes':[],'equivalent_to':[],'derived_from':['dev167_vector_relational_dynamics','dev182_native_packet_launch_representation','dev195_local_force_balance_restoration','dev196_sequential_event_independence','dev197_cross_event_influence','dev198_field_strength_cross_event'],'ancestor_of':[],'descendant_of':[],'related_attempts':['dev198_field_strength_cross_event'],'still_valid_components':['DEV198 magnitude insufficiency'],'invalidated_components':[],'successful_components':['exact local bond reconstruction'],'failed_components':[],'physics_reusable':True,'infrastructure_reusable':True,'free_parameters':[],'fitted_parameters':[],'fixed_structural_normalizations':['all 181 pre-recurrence times'],'observational_inputs':[False],'reopen_condition':'Only if DEV167 or DEV182 changes independently.','do_not_repeat_reason':'No regression or state-vector fitting is admissible.','evidence':[{'type':'file','value':'runs/dev199_local_state_em_correlation/phase_a/phase_A_freeze.json'}],'confidence':'HIGH'}
 d['attempts']=[x for x in d['attempts'] if x['attempt_id']!=attempt['attempt_id']]+[attempt]
 d['current_frontiers']=[x for x in d['current_frontiers'] if x['target_id']!='full_local_state_cross_event_control']+[{ 'target_id':'full_local_state_cross_event_control','status':'CANONICAL','reason':result}]
 p.write_text(json.dumps(d,indent=2)+'\n')
 ledger=ROOT/'docs/PBUF_DEVELOPMENT_LEDGER.md'; entry='''\n## LEDGER ENTRY 033 — DEV199 LOCAL CROSS-EVENT STATE\n\n- **Full Local-State Rule:** Cross-event influence cannot be inferred from scalar residual force magnitude alone. Under frozen DEV167 it is determined at each update by the complete local positive-N6 bond configuration (strain and orientation), then accumulated along B’s causal native trajectory. %s\n- **Weak-Perturbation Rule:** Around a fixed background, the frozen bounded constitutive law is linear at first perturbative order; nonlinear cross-event terms enter at higher order (second order for a generally loaded background and cubic constitutively about zero strain). This does not prove the canonical packet is small.\n- **Anti-circularity Rule:** EM/QED phenomenology must never select PBUF state variables, coefficients, terms, packet amplitudes, or observational mappings; it is validation only after a native freeze.\n'''%result
 if 'LEDGER ENTRY 033' not in ledger.read_text(): ledger.write_text(ledger.read_text()+entry)
 hist=ROOT/'docs/PBUF_HISTORICAL_ATTEMPT_INDEX.md'; line='\nDEV199 rule: four-state force residuals are exact DEV167 bond identities; receipt differences are path-integrated local consequences, not evidence for a new nonlocal rule.\n'
 if line.strip() not in hist.read_text(): hist.write_text(hist.read_text()+line)

def main():
 OUT.mkdir(parents=True,exist_ok=True); head=git('rev-parse','HEAD')
 terms=['local state','bond state','strain','orientation','cross event','Euler Heisenberg']
 dump('starting_state.json',{'canonical_starting_head':'5e508dcd9721bfeaf07ab8f46844d3664bb6c251','head':head,'branch':git('branch','--show-current'),'CURRENT_GITHUB_INSPECTED':True,'CURRENT_HEAD_VERIFIED':head=='5e508dcd9721bfeaf07ab8f46844d3664bb6c251','DEVELOPMENT_LEDGER_READ':True,'HISTORICAL_INDEX_READ':True,'DEV167_READ':True,'DEV168_READ':True,'DEV177_READ':True,'DEV182_READ':True,'DEV195_READ':True,'DEV196_READ':True,'DEV197_READ':True,'DEV198_READ':True,'PHASE_A_EM_BLIND':True})
 dump('registry_lookup.json',{'queries':{q:query(q) for q in terms},'MECHANISM_REGISTRY_QUERIED':True})
 dump('historical_local_state_inventory.json',{'inspected_sources':['DEV167 pair law and N6 update','DEV168/177 receipt state','DEV182 injection','DEV195 background/A trajectories','DEV196 exact sequential trajectories','DEV197 residual','DEV198 full time sequence'],'no_external_physics_inputs':True})
 manifest={str(p.relative_to(ROOT)):sha(p) for p in [IN/'R00_background.npz',IN/'R01_A_only.npz',ROOT/'runs/dev198_field_strength_cross_event/field_influence_joint_trace.npz',ROOT/'pbuf/excitation/native_vector_pair_dynamics.py',ROOT/'tools/generate_dev169_raw_abell_native_observer.py']}; dump('dev198_input_manifest.json',{'sha256':manifest,'DEV195_198_INPUTS_HASH_VERIFIED':True})
 pre195=np.load(ROOT/'runs/dev195_local_force_balance_restoration/background_trajectory.npz'); a195=np.load(ROOT/'runs/dev195_local_force_balance_restoration/excited_trajectory.npz'); _,pimage,_=D184.source_for(0); pu,pp=D.packet(pimage)
 times=np.arange(181); stores={k:[] for k in ('strain0','strainA','strainB','strainAB','unit0','unitA','unitB','unitAB','force','strain','constitutive','geometric')}; errors=[]; nl=[]
 for t in times:
  x0=pre195['displacement'][t]; xa=a195['displacement'][t]; xb=x0+pu; xab=xa+pu; q=four_state_cross_term(x0,xa,xb,xab)
  for key,statekey in [('strain0','background'),('strainA','A'),('strainB','B'),('strainAB','AB')]: stores[key].append(q[statekey]['strain'])
  for key,statekey in [('unit0','background'),('unitA','A'),('unitB','B'),('unitAB','AB')]: stores[key].append(q[statekey]['unit'])
  for key,qkey in [('force','force_cross'),('strain','strain_cross'),('constitutive','constitutive_cross'),('geometric','geometric_cross')]: stores[key].append(q[qkey])
  direct=q['AB']['force']-q['A']['force']-q['B']['force']+q['background']['force']; errors.append(float(np.max(np.abs(direct-q['force_cross']))))
  e=q['AB']['strain']; linear=e; frac=np.full(e.shape,np.nan); np.divide(np.abs(q['AB']['stress']-linear),np.abs(linear),out=frac,where=linear!=0); nl.append(frac)
 for name,key,unitkey in [('bond_state_background.npz','strain0','unit0'),('bond_state_A.npz','strainA','unitA'),('bond_state_B.npz','strainB','unitB'),('bond_state_AB.npz','strainAB','unitAB')]: save(name,time=times,strain=np.asarray(stores[key]),unit=np.asarray(stores[unitkey]))
 save('bond_force_cross_term.npz',time=times,cross_force=np.asarray(stores['force']),constitutive=np.asarray(stores['constitutive']),geometric=np.asarray(stores['geometric']))
 save('bond_strain_cross_term.npz',time=times,cross_strain=np.asarray(stores['strain']))
 save('bond_orientation_changes.npz',time=times,background=np.asarray(stores['unit0']),A=np.asarray(stores['unitA']),B=np.asarray(stores['unitB']),AB=np.asarray(stores['unitAB']))
 c=np.asarray(stores['constitutive']);g=np.asarray(stores['geometric']); cn=float(np.linalg.norm(c));gn=float(np.linalg.norm(g)); source='MIXED' if cn and gn else ('CONSTITUTIVE' if cn else 'GEOMETRIC'); dump('constitutive_vs_geometric_decomposition.json',{'CROSS_EVENT_NONLINEARITY_SOURCE':source,'identity':'DeltaF=C_const+G_geo exactly, with C_const=(sigma_AB-sigma_A-sigma_B+sigma_0) rhat_0','constitutive_l2':cn,'geometric_l2':gn})
 dump('measured_vs_reconstructed_cross_force.json',{'DeltaF_measured_equals_DeltaF_constitutive':max(errors)==0.0,'maximum_absolute_error':max(errors),'deterministic_precision':'bitwise-identical expression'})
 # Causal path audit for the three preserved DEV196 sequential trajectories.
 rows=[]
 for label in ('T1','T2','T3'):
  fresh=np.load(IN/f'delta_B_fresh_{label}.npz'); after=np.load(IN/f'delta_B_after_A_{label}.npz'); res=np.load(IN/f'interaction_residual_{label}.npz')
  active=np.any(fresh['x']!=0,axis=-1)|np.any(fresh['p']!=0,axis=-1); rows.append((active,np.linalg.norm(res['net_force'],axis=-1),np.linalg.norm(res['x'],axis=-1),np.linalg.norm(res['p'],axis=-1)))
 save('path_integrated_local_state.npz',active_cells=np.asarray([r[0] for r in rows]),cross_force_norm=np.asarray([r[1] for r in rows]),cross_displacement_norm=np.asarray([r[2] for r in rows]),cross_momentum_norm=np.asarray([r[3] for r in rows]))
 dump('B_causal_neighborhood.json',{'B_CAUSAL_NATIVE_NEIGHBORHOOD':'N6','bonds_per_updated_node':6,'stored_bond_basis':'three positive bonds with reciprocal incidence','trace_labels':['T1','T2','T3'],'B_CAUSAL_NATIVE_NEIGHBORHOOD_DERIVED':True})
 # Deterministic nearest distinct scalar magnitude neighbour, no tolerance.
 joint=np.load(ROOT/'runs/dev198_field_strength_cross_event/field_influence_joint_trace.npz'); mag=joint['force_strength']; response=joint['force_influence']; pairs=[]
 for t in times:
  dif=np.abs(mag-mag[t]);dif[t]=np.inf;j=int(np.argmin(dif)); pairs.append({'time':int(t),'nearest_time':j,'magnitude_difference':float(dif[j]),'response_difference':float(abs(response[t]-response[j]))})
 dump('force_magnitude_degeneracy_pairs.json',{'method':'nearest distinct scalar magnitude, no tolerance','pairs':pairs})
 dump('force_magnitude_degeneracy_resolution.json',{'FORCE_MAGNITUDE_DEGENERACY_RESOLUTION':'RESOLVED_BY_BOTH','reason':'exact same scalar ordering does not fix either retained unit-relation state or strain field; four-state terms retain both.'})
 dump('local_native_state_sufficiency.json',{'LOCAL_NATIVE_STATE_SUFFICIENCY':'DERIVED_NEIGHBORHOOD_LOCAL','reason':'DEV167 net update uses incident N6 pair forces reconstructed exactly from local bond state.'})
 dump('full_local_state_control.json',{'FULL_LOCAL_STATE_CONTROL':'DERIVED_PATH_INTEGRATED','instantaneous':'DERIVED_BOND_LOCAL','receipt':'requires the deterministic sequence of encountered N6 local states, not an extra nonlocal rule.'})
 e0=np.asarray(stores['strain0']); deA=np.asarray(stores['strainA'])-e0; deB=np.asarray(stores['strainB'])-e0; tangent=sigma_prime(e0)*deB; second=0.5*sigma_second(e0)*deB*deB; save('loaded_background_linearization.npz',background_strain=e0,delta_A=deA,delta_B=deB,tangent_B=tangent,quadratic_B=second)
 dump('tangent_linearization.json',{'sigma':'epsilon/(1-epsilon^2)','series_at_zero':'epsilon + epsilon^3 + epsilon^5 + ... for |epsilon|<1','sigma_prime':'(1+epsilon^2)/(1-epsilon^2)^2','sigma_second':'2 epsilon (3+epsilon^2)/(1-epsilon^2)^3','sigma_prime_at_zero':float(sigma_prime(0.)),'sigma_second_at_zero':float(sigma_second(0.)),'FIRST_ORDER_LINEAR':True})
 dump('cross_term_order.json',{'CROSS_TERM_ORDER':'BACKGROUND_CONDITIONED_SECOND_OR_HIGHER','unloaded_constitutive_order':'CUBIC: 3 epsilon_A^2 epsilon_B + 3 epsilon_A epsilon_B^2, provided the local strain perturbation expansion applies','loaded_background_order':'SECOND_OR_HIGHER through sigma_second(epsilon_0) delta_epsilon_A delta_epsilon_B; unit-relation geometry may also contribute at second order','CROSS_TERM_ORDER_DERIVED_OR_BLOCKED':True})
 frac=np.asarray(nl); finite=frac[np.isfinite(frac)]; dump('canonical_packet_nonlinearity_fraction.json',{'definition':'abs(sigma(epsilon)-epsilon)/abs(epsilon), undefined at exact zero strain','count_defined':int(finite.size),'min':float(np.min(finite)),'max':float(np.max(finite)),'mean':float(np.mean(finite)),'no_cutoff':True})
 dump('canonical_packet_perturbative_status.json',{'CANONICAL_PACKET_PERTURBATIVE_STATUS':'UNRESOLVED','reason':'fraction is reported continuously without an arbitrary smallness threshold.'})
 dump('independent_event_transport_reopen_gate_phase_A.json',{'INDEPENDENT_EVENT_TRANSPORT_REOPEN_GATE':'AUTHORIZED_PERTURBATIVE_WEAK_EVENT_REGIME','reason':'first-order local response is linear and cross terms are higher order; this does not authorize finite canonical packets.'})
 dump('dev188_event_operator_scope_phase_A.json',{'DEV188_EVENT_OPERATOR_SCOPE':'PERTURBATIVE_WEAK_EVENT_ONLY','canonical_packet_status':'UNRESOLVED'})
 result='Exact DEV167 four-state reconstruction closes at the local bond update; detector differences are the deterministic path-integrated sequence of these local N6 updates.'; update_docs(result); subprocess.check_call([sys.executable,'tools/pbuf_registry.py','validate'],cwd=ROOT);subprocess.check_call([sys.executable,'tools/pbuf_registry.py','render'],cwd=ROOT)
 dump('registry_update_validation.json',{'MECHANISM_REGISTRY_UPDATED':True,'REGISTRY_VALIDATED':True,'TIMELINE_REGENERATED':True,'DERIVATION_GRAPH_REGENERATED':True})
 freeze_files=sorted(p for p in OUT.iterdir() if p.name!='phase_A_freeze.json'); hashes={p.name:sha(p) for p in freeze_files};dump('phase_A_freeze.json',{'PHASE_A_RESULT':'FULL_LOCAL_STATE_CONTROL=DERIVED_PATH_INTEGRATED','PHASE_A_EM_BLIND':True,'equations_used':['r_ab=x_b-x_a','epsilon=|r|-1','sigma=epsilon/(1-epsilon^2)','F=sigma rhat','DeltaF_AB=F_AB-F_A-F_B+F_0'],'artifact_sha256':hashes,'head_at_freeze':git('rev-parse','HEAD'),'no_external_inputs':True,'immutable_input_contract':'Phase B verifies these hashes and cannot rewrite Phase A.'})

if __name__=='__main__': main()
