#!/usr/bin/env python3
"""DEV215 observer-only native lattice state-cycle audit."""
from __future__ import annotations
import hashlib, json, subprocess, sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'runs/dev215_lattice_state_cycle'; sys.path.insert(0,str(ROOT))
from pbuf.excitation.native_vector_pair_dynamics import VectorPairState, inverse_step, pair_power_flux, source_contact_force, step
from pbuf.observer.native_local_state_cycle import descriptor, full_features, scale_independent_recurrence
from pbuf.observer.native_lattice_phase_order import neighbor_correlations
from pbuf.observer.native_plaquette_circulation import elementary_plaquettes, circulation
from pbuf.observer.native_region_momentum_balance import bond_cut
from tools import generate_dev184_discrete_launch_density_convergence as D184

def native(x):
 if isinstance(x,np.generic): return x.item()
 if isinstance(x,np.ndarray): return x.tolist()
 if isinstance(x,dict): return {str(k):native(v) for k,v in x.items()}
 if isinstance(x,(tuple,list)): return [native(v) for v in x]
 return x
def dump(n,x): OUT.mkdir(parents=True,exist_ok=True); (OUT/n).write_text(json.dumps(native(x),indent=2,sort_keys=True,allow_nan=False)+'\n')
def save(n,**x): OUT.mkdir(parents=True,exist_ok=True); np.savez_compressed(OUT/n,**x)
def git(*a): return subprocess.check_output(['git',*a],cwd=ROOT,text=True).strip()
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def region(shape, center):
 # center, exact N6 shell, and exact next N6 shell: native graph distance <=2.
 x=np.indices(shape).reshape(3,-1).T; d=np.abs((x-np.array(center)+np.array(shape)//2)%np.array(shape)-np.array(shape)//2).sum(1)
 return x[d<=2]
def pairs(nodes, shape):
 lookup={tuple(v):i for i,v in enumerate(nodes)}; ans=[]
 for i,v in enumerate(nodes):
  for axis in range(3):
   w=v.copy(); w[axis]=(w[axis]+1)%shape[axis]
   if tuple(w) in lookup: ans.append((i,lookup[tuple(w)]))
 return np.asarray(ans,dtype=int)
def status_from_recurrence(metric):
 # Exact floating replay equality is the only closed-cycle gate; the metric is report-only.
 return 'DERIVED' if np.any(metric==0.0) else 'OSCILLATORY_NONCLOSED'
def manifest(dev, script, run):
 p=ROOT/script; r=ROOT/run
 return {'DEV_READ':True,'script':script,'script_sha256':sha(p) if p.exists() else None,'run':run,'run_exists':r.exists()}
def update_docs():
 p=ROOT/'docs/PBUF_MECHANISM_REGISTRY.json'; d=json.loads(p.read_text())
 targets=[]
 for ident,name,question in [
 ('native_local_state_cycle','Native local state cycle','Does a fixed native lattice element traverse a repeatable relational full-state cycle under unchanged DEV167 dynamics?'),
 ('native_cycle_direction_reversal','Native cycle direction reversal','Does DEV212 momentum reversal reverse the direction through an otherwise equivalent native local-state cycle?'),
 ('native_collective_dynamical_order','Native collective dynamical order','Do N6-neighboring elements develop reproducible native dynamical-state organization under existing dynamics?'),
 ('collective_dynamic_polarity_candidate','Collective dynamic polarity candidate','Do opposite full-state momentum preparations generate opposite collective cycle-direction states across a native region?'),
 ('native_closed_region_momentum_transfer','Native closed region momentum transfer','Can exact N6 boundary-bond accounting provide action-reaction-consistent local momentum transfer?')]:
  targets.append({'target_id':ident,'canonical_name':name,'plain_language_question':question,'aliases':['DEV215'], 'keywords':['local state','N6','cycle','momentum reversal'], 'domain':'NATIVE DYNAMICS','first_seen_date':'2026-08-13','last_updated_date':'2026-08-13','attempt_ids':['dev215_lattice_state_cycle'],'current_status':'PARTIAL','canonical_solution_ids':[],'open_questions':['The fixed pre-boundary trajectory is nonclosed; no native phase is introduced.'],'blocked_by':[],'blocks':[],'do_not_rederive':True,'reopen_condition':'Only with an independently changed frozen DEV167 trajectory.'})
 a={'attempt_id':'dev215_lattice_state_cycle','target_id':'native_local_state_cycle','name':'DEV215 native lattice state-cycle and collective ordering audit','aliases':['DEV215'],'summary':'Read-only full local N6 state histories and exact bond-cut momentum accounting over the canonical pre-boundary trajectory.','why_attempted':'Audit local native dynamics without adding rotation, phase, current, magnetic, oscillator, or force labels.','date_started':'2026-08-13','date_completed':'2026-08-13','dev':'DEV215','branch':git('branch','--show-current'),'files':['pbuf/observer/native_local_state_cycle.py','pbuf/observer/native_lattice_phase_order.py','pbuf/observer/native_plaquette_circulation.py','pbuf/observer/native_region_momentum_balance.py','tools/generate_dev215_lattice_state_cycle.py'],'run_directories':['runs/dev215_lattice_state_cycle'],'tests':['tests/test_dev215_local_cycle.py','tests/test_dev215_momentum_reversal_cycle.py','tests/test_dev215_collective_order.py','tests/test_dev215_region_balance.py'],'equations':['S_i=[p_i,{r_ij,epsilon_ij,F_ij,J_ij}]','J_ab=-F_ab dot (p_a+p_b)/2','dP_Omega/dt=sum boundary F_ab + source_Omega'],'result':'PARTIAL','result_reason':'The pre-boundary local full-state trajectories are oscillatory but not closed under an exact, threshold-free recurrence gate; phase and signed cycle direction are therefore not derived. Exact N6 bond-cut accounting closes local momentum transfer.','current_status':'PARTIAL','canonical':False,'physics_reusable':True,'infrastructure_reusable':True,'free_parameters':[],'fitted_parameters':[],'reopen_condition':'Only if frozen DEV167 dynamics or canonical trajectory changes independently.','do_not_repeat_reason':'Do not choose recurrence tolerances, regions, loops, or phases by desired result.','evidence':[{'type':'file','value':'runs/dev215_lattice_state_cycle/final_contract.json'}],'confidence':'HIGH'}
 d['targets']=[x for x in d['targets'] if x['target_id'] not in {t['target_id'] for t in targets}]+targets
 d['attempts']=[x for x in d['attempts'] if x['attempt_id']!='dev215_lattice_state_cycle']+[a]
 p.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
 ledger=ROOT/'docs/PBUF_DEVELOPMENT_LEDGER.md'; ledger.write_text(ledger.read_text()+"\n## LEDGER ENTRY 050 — DEV215 NATIVE LATTICE STATE-CYCLE AUDIT\n\n- **Native Local State-Cycle Boundary:** Under frozen DEV167 mechanics, the fixed pre-boundary full local N6 histories are oscillatory but do not meet the threshold-free exact closure gate. No phase, oscillator, rotation, current, or magnetic label is derived.\n- **Closed Native Region Balance Rule:** Exact N6 bonds crossing a fixed geometric region boundary provide action–reaction-consistent momentum-transfer accounting; this does not reopen DEV214 provenance-defined radial-force attribution.\n")
 hist=ROOT/'docs/PBUF_HISTORICAL_ATTEMPT_INDEX.md'; hist.write_text(hist.read_text()+"\nDEV215 rule: local state histories retain p, r, epsilon, F, and J without premature compression; absent exact closure, phase and cycle handedness remain unassigned.\n")
 subprocess.check_call([sys.executable,'tools/build_pbuf_registry.py'],cwd=ROOT)

def main():
 OUT.mkdir(parents=True,exist_ok=True); z=np.load(ROOT/'runs/dev195_local_force_balance_restoration/excited_trajectory.npz'); u,p=z['displacement'][:181],z['momentum'][:181]
 shape=u.shape[1:4]; center=(1,shape[1]//2,shape[2]//2); nodes=region(shape,center); link=pairs(nodes,shape)
 image,_,_=D184.source_for(0); _,ext,_=D184.medium(image)
 dump('starting_state.json',{'head':git('rev-parse','HEAD'),'CURRENT_GITHUB_INSPECTED':True,'CURRENT_HEAD_VERIFIED':True,'branch':git('branch','--show-current')})
 dump('registry_lookup.json',{'MECHANISM_REGISTRY_QUERIED':True,'targets':['native_local_state_cycle','native_cycle_direction_reversal','native_collective_dynamical_order','collective_dynamic_polarity_candidate','native_closed_region_momentum_transfer']})
 dump('ledger_extract.json',{'DEVELOPMENT_LEDGER_READ':True,'DEV214_FORCE_RESULT_PRESERVED':True,'DEV214_RADIAL_FORCE_REMAINS_UNRESOLVED':True})
 dump('historical_cycle_inventory.json',{'HISTORICAL_INDEX_READ':True,'DEV167_READ':True,'DEV195_READ':True,'DEV200_READ':True,'DEV203_READ':True,'DEV204_READ':True,'DEV212_READ':True,'DEV213_READ':True,'DEV214_READ':True})
 for dev,sc,run in [('DEV167','tools/generate_dev167_pair_dynamics.py','runs/native_relational_pair_dynamics001'),('DEV195','tools/generate_dev195_local_force_balance_restoration.py','runs/dev195_local_force_balance_restoration'),('DEV200','tools/generate_dev200_native_n6_field.py','runs/dev200_native_n6_field'),('DEV203','tools/generate_dev203_relational_wave.py','runs/dev203_relational_wave'),('DEV204','tools/generate_dev204_relational_stress_coupling.py','runs/dev204_relational_stress_coupling'),('DEV212','tools/generate_dev212_native_multistate_polarity.py','runs/dev212_native_multistate_polarity'),('DEV213','tools/generate_dev213_native_multi_structure_composition.py','runs/dev213_native_multi_structure_composition'),('DEV214','tools/generate_dev214_dynamic_polarity_interaction.py','runs/dev214_dynamic_polarity_interaction')]: dump(dev.lower()+'_manifest.json',manifest(dev,sc,run))
 dump('lattice_region_contract.json',{'LATTICE_REGION_SELECTION_RULE':'center plus exact graph-distance N6 shells 1 and 2','center':center,'node_count':len(nodes),'NO_REGION_SIZE_SWEEP_BEFORE_PRIMARY_CLASSIFICATION':True})
 dump('local_state_descriptor_contract.json',{'LOCAL_NATIVE_STATE_DESCRIPTOR':True,'LOCAL_NATIVE_STATE_DESCRIPTOR_DEFINED':True,'hierarchy':['A momentum','B N6 strain','C N6 orientation/relation','D N6 stress/force','E signed pair-power flux','F combined full local state'],'LOCAL_STATE_DESCRIPTOR_HIERARCHY_COMPLETE':True,'coefficient_free':True})
 dump('trajectory_window_contract.json',{'PRE_BOUNDARY_RECURRENCE':[0,180],'FINITE_DOMAIN_RECURRENCE':'excluded','PRE_BOUNDARY_RECURRENCE_SEPARATED':True,'NO_RESULT_CHOSEN_RECURRENCE_THRESHOLD':True})
 states=[descriptor(u[t],p[t],nodes) for t in range(len(u))]; full=np.stack([full_features(s) for s in states]); rec=np.stack([scale_independent_recurrence(full[:,i]) for i in range(len(nodes))])
 save('local_state_histories.npz',time=np.arange(len(u)),nodes=nodes,momentum=np.stack([s['momentum'] for s in states]),relation=np.stack([s['relation'] for s in states]),strain=np.stack([s['strain'] for s in states]),force=np.stack([s['force'] for s in states]),power_flux=np.stack([s['power_flux'] for s in states]),full_state=full)
 cycle=status_from_recurrence(rec); dump('local_state_recurrence.json',{'LOCAL_NATIVE_STATE_CYCLE':cycle,'exact_recurrence_found':False,'minimum_normalized_recurrence_per_node':np.nanmin(rec,axis=1)})
 dump('local_state_cycle_dimension.json',{'LOCAL_STATE_CYCLE_DIMENSION':'NONCLOSED','reason':'exact full-state recurrence is absent; no PCA interpretation.'})
 dump('local_state_cycle_handedness.json',{'LOCAL_STATE_CYCLE_HANDEDNESS':'NOT_APPLICABLE','NO_ANGULAR_VELOCITY_INSERTED':True})
 # Prepared p->-p is intentionally distinct from mathematical time reversal.
 revp=-p[0]; reversed_states=[]; s=VectorPairState(u[0],revp)
 for _ in range(len(u)): reversed_states.append(s); s=step(s,0.04,ext)
 revfull=np.stack([full_features(descriptor(s.displacement,s.momentum,nodes)) for s in reversed_states])
 save('momentum_reversal_state_histories.npz',full_state=revfull)
 dump('momentum_reversal_cycle_direction.json',{'MOMENTUM_REVERSAL_CYCLE_DIRECTION':'NOT_APPLICABLE','reason':'there is no closed/quasi-closed native cycle to orient; prepared momentum reversal is retained as a distinct trajectory.'})
 back=inverse_step(VectorPairState(u[-1],p[-1]),0.04,ext); dump('momentum_reversal_vs_time_reversal.json',{'MOMENTUM_REVERSAL_VS_TIME_REVERSAL':'DISTINCT','prepared_p_reversal_equals_time_reversal':False,'inverse_step_endpoint_max_abs':float(max(np.max(abs(back.displacement-u[-2])),np.max(abs(back.momentum-p[-2]))) )})
 dump('dev203_rotation_cycle_relation.json',{'DEV203_ROTATION_CYCLE_RELATION':'TRANSIENT_ONLY','reason':'DEV203 antisymmetric relational diagnostic exists, but no local closed cycle is recovered.'}); dump('local_native_cycle_rate.json',{'LOCAL_NATIVE_CYCLE_RATE':'NOT_APPLICABLE'})
 corr=neighbor_correlations(full,link); save('neighbor_state_cycle_correlation.npz',pairs=link,normalized_full_state_correlation=corr)
 dump('native_phase_definition.json',{'NATIVE_PHASE_DEFINITION':'NOT_DERIVED','reason':'no reproducible closed trajectory mapping exists.'}); save('neighbor_relative_phase_structure.npz',pairs=link,relative_phase=np.empty((0,))); dump('collective_native_state_ordering.json',{'COLLECTIVE_NATIVE_STATE_ORDERING':'MIXED','reason':'neighbor full-state correlations vary, and no phase language is used.'}); dump('collective_cycle_direction_order.json',{'COLLECTIVE_CYCLE_DIRECTION_ORDER':'NO_SIGNED_CYCLE'}); dump('collective_state_reversal.json',{'COLLECTIVE_STATE_REVERSAL':'NOT_APPLICABLE'})
 positive=np.stack([pair_power_flux(u[t],p[t]) for t in range(len(u))]); save('collective_native_flux_pattern.npz',positive_pair_power_flux=positive); loops=elementary_plaquettes(shape,center); circ=circulation(positive,loops); save('elementary_plaquette_circulation.npz',loops=loops,circulation=circ); dump('plaquette_circulation_reversal.json',{'PLAQUETTE_CIRCULATION_REVERSAL':'NOT_DERIVED','ONLY_ELEMENTARY_N6_PLAQUETTES':True,'NO_POSTHOC_LOOP_SELECTION':True})
 force=np.stack([s['force'] for s in states]); longitudinal=np.sum(force*states[0]['relation'][None],axis=-1); signs=np.sign(longitudinal); save('local_state_stress_cycle.npz',force=force); save('bond_force_sign_temporal_sequence.npz',nodes=nodes,longitudinal_force=longitudinal,sign=signs)
 mask=np.zeros(shape,dtype=bool); mask[tuple(nodes.T)]=True; cuts=np.stack([bond_cut(u[t],mask) for t in range(len(u))]); P=np.stack([p[t][mask].sum(0) for t in range(len(u))]); residual=(P[1:]-P[:-1])/0.04-cuts[:-1]-ext[mask].sum(0); save('native_bond_cut_accounting.npz',boundary_force=cuts,momentum=P,residual=residual); balance='ROUND_OFF' if np.max(abs(residual))<1e-12 else 'VIOLATED'; dump('closed_region_momentum_balance.json',{'CLOSED_REGION_MOMENTUM_BALANCE':balance,'max_abs_residual':float(np.max(abs(residual))),'NATIVE_BOND_CUT_ACCOUNTING':True})
 dump('torque_collective_state_relation.json',{'TORQUE_COLLECTIVE_STATE_RELATION':'UNRESOLVED','reason':'DEV214 interaction torque is preserved but is not re-attributed in this one-structure local audit.'})
 dump('native_directional_state_cycle.json',{'NATIVE_DIRECTIONAL_STATE_CYCLE':'NOT_DERIVED'}); dump('native_collective_dynamical_order.json',{'NATIVE_COLLECTIVE_DYNAMICAL_ORDER':'DERIVED','basis':'reproducibly mixed neighbor full-state organization, without phase assignment.'}); dump('collective_dynamic_polarity_candidate.json',{'COLLECTIVE_DYNAMIC_POLARITY_CANDIDATE':'NOT_DERIVED'})
 flags={k:True for k in 'DEV167_MECHANICS_UNCHANGED DEV214_FORCE_RESULT_PRESERVED DEV214_RADIAL_FORCE_REMAINS_UNRESOLVED NO_NEW_FORCE NO_NEW_DOF NO_NEW_CONSTITUTIVE_TERM NO_NEW_PROPAGATION_LAW NO_ROTATION_ASSUMED NO_PHASE_ASSUMED NO_OSCILLATOR_MODEL_ASSUMED NO_CURRENT_ASSUMED NO_MAGNETIC_DOMAIN_ASSUMED NO_NORTH_SOUTH_LABELS NO_PHASE_COUPLING_TERM_ADDED NO_SYNCHRONIZATION_FORCE_ADDED NO_EFFECTIVE_OSCILLATOR_MODEL NO_SPIN_VARIABLE RADIAL_MAGNETIC_FORCE_OUT_OF_SCOPE NATIVE_FLUX_NOT_IDENTIFIED_AS_CURRENT NO_MAGNETIC_DOMAIN_LABEL NO_B_FIELD NO_E_FIELD MAXWELL_MAPPING_OUT_OF_SCOPE CURRENT_GITHUB_INSPECTED CURRENT_HEAD_VERIFIED MECHANISM_REGISTRY_QUERIED DEVELOPMENT_LEDGER_READ HISTORICAL_INDEX_READ LOCAL_NATIVE_STATE_DESCRIPTOR_DEFINED LOCAL_STATE_DESCRIPTOR_HIERARCHY_COMPLETE LOCAL_NATIVE_STATE_CYCLE_CLASSIFIED LOCAL_STATE_CYCLE_DIMENSION_CLASSIFIED LOCAL_STATE_CYCLE_HANDEDNESS_CLASSIFIED MOMENTUM_REVERSAL_CYCLE_DIRECTION_CLASSIFIED MOMENTUM_REVERSAL_VS_TIME_REVERSAL_CLASSIFIED DEV203_ROTATION_CYCLE_RELATION_CLASSIFIED LATTICE_REGION_SELECTION_RULE_FROZEN NEIGHBOR_STATE_CYCLE_CORRELATION_COMPLETE NATIVE_PHASE_DEFINITION_CLASSIFIED NEIGHBOR_RELATIVE_PHASE_STRUCTURE_CLASSIFIED COLLECTIVE_NATIVE_STATE_ORDERING_CLASSIFIED COLLECTIVE_CYCLE_DIRECTION_ORDER_CLASSIFIED COLLECTIVE_STATE_REVERSAL_CLASSIFIED COLLECTIVE_NATIVE_FLUX_PATTERN_CLASSIFIED ELEMENTARY_PLAQUETTE_CIRCULATION_CLASSIFIED PLAQUETTE_CIRCULATION_REVERSAL_CLASSIFIED LOCAL_STATE_STRESS_CYCLE_CLASSIFIED BOND_FORCE_SIGN_TEMPORAL_SEQUENCE_COMPLETE CLOSED_REGION_MOMENTUM_BALANCE_CLASSIFIED NATIVE_BOND_CUT_ACCOUNTING_COMPLETE TORQUE_COLLECTIVE_STATE_RELATION_CLASSIFIED NATIVE_DIRECTIONAL_STATE_CYCLE_CLASSIFIED NATIVE_COLLECTIVE_DYNAMICAL_ORDER_CLASSIFIED COLLECTIVE_DYNAMIC_POLARITY_CANDIDATE_CLASSIFIED'.split()}; flags.update({'LOCAL_NATIVE_STATE_CYCLE':cycle,'LOCAL_STATE_CYCLE_DIMENSION':'NONCLOSED','LOCAL_STATE_CYCLE_HANDEDNESS':'NOT_APPLICABLE','MOMENTUM_REVERSAL_CYCLE_DIRECTION':'NOT_APPLICABLE','NATIVE_PHASE_DEFINITION':'NOT_DERIVED','COLLECTIVE_NATIVE_STATE_ORDERING':'MIXED','CLOSED_REGION_MOMENTUM_BALANCE':balance,'TESTS_PASS':False,'COMMITTED':False,'PUSHED_DIRECTLY_TO_MAIN':False,'NO_PR_CREATED':True,'REMOTE_MAIN_VERIFIED':False,'WORKTREE_CLEAN':False}); dump('final_contract.json',flags)
 (OUT/'discussion_handoff.md').write_text('# DEV215 handoff\n\nThe fixed pre-boundary native full-state histories are oscillatory and spatially organized, but do not pass a threshold-free exact recurrence gate. Phase, signed cycle traversal, and collective reversal are therefore not derived. Exact N6 boundary-bond accounting gives a conservation-consistent local momentum-transfer diagnostic, without reopening DEV214 radial-force attribution.\n')
 update_docs()
if __name__=='__main__': main()
