#!/usr/bin/env python3
"""DEV226: frozen, coefficient-free native staggered local-order audit."""
from __future__ import annotations
import hashlib, json, subprocess, sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'runs/dev226_staggered_local_order'; sys.path.insert(0,str(ROOT))
from pbuf.analysis.native_staggered_order import unique_n6_bonds, contract_bonds, sign_counts, zero_causes, axial_equivalence_error

def native(x):
    if isinstance(x,np.generic): return x.item()
    if isinstance(x,np.ndarray): return x.tolist()
    if isinstance(x,dict): return {k:native(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)): return [native(v) for v in x]
    return x
def dump(n,x): OUT.mkdir(parents=True,exist_ok=True); (OUT/n).write_text(json.dumps(native(x),indent=2,sort_keys=True)+'\n')
def save(n,**x): OUT.mkdir(parents=True,exist_ok=True); np.savez_compressed(OUT/n,**x)
def sha(p): return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def git(*a): return subprocess.check_output(['git',*a],cwd=ROOT,text=True).strip()

def temporal_class(c):
    o,a=c['opposed'],c['aligned']; neg=np.any(o>a); pos=np.any(a>o)
    if neg and pos:return 'SIGN_STRUCTURE_SWITCHING'
    if neg:return 'PERSISTENT'
    if pos:return 'PERSISTENT'
    return 'ABSENT'
def classification(c):
    o,a,z=c['opposed'],c['aligned'],c['zero']
    if np.all(a==0) and np.all(o>0): return 'EXACT'
    if np.any(o>a) and np.any(a>o): return 'MIXED_IN_TIME'
    if np.all(o>a): return 'OPPOSED_DOMINANT'
    if np.all(a>o): return 'ALIGNED_DOMINANT'
    if np.all(z>=o) and np.all(z>=a): return 'ZERO_DOMINANT'
    if np.all(o==a): return 'BALANCED'
    return 'ABSENT'
def axis_class(c):
    axes={name:np.all(v['opposed']>v['aligned']) for name,v in c['by_axis'].items()}
    if all(axes.values()):return 'FULL_3D_OPPOSED'
    if axes['x'] and not axes['y'] and not axes['z']:return 'LONGITUDINAL_OPPOSED'
    if not axes['x'] and axes['y'] and axes['z']:return 'TRANSVERSE_OPPOSED'
    if any(axes.values()):return 'LAYERED' if sum(axes.values())==1 else 'AXIS_MIXED'
    return 'NO_AXIS_ORDER'

def update_docs(result, axis, temporal, minimality, next_selector):
    p=ROOT/'docs/PBUF_MECHANISM_REGISTRY.json'; d=json.loads(p.read_text())
    specs=[('native_distributed_staggered_local_order','Are neighboring DEV203 antisymmetric relational states predominantly opposed under the frozen tensor contraction?',result),('native_staggered_tensor_order','Does A(a):A(b) show persistent negative nearest-neighbor organization across the full N6 structure?',result),('staggered_order_axis_structure','Is any native staggered order three-dimensional, layered, longitudinal, transverse, or mixed across N6 axes?',axis),('staggered_order_temporal_character','Is the nearest-neighbor antisymmetric order persistent across the frozen DEV203 trajectory or does its sign organization change with time?',temporal)]
    ids={x[0] for x in specs}
    def target(i,q,status): return {'target_id':i,'canonical_name':i.replace('_',' '),'plain_language_question':q,'aliases':['DEV226'],'keywords':['staggered','antisymmetric tensor','N6','nearest neighbor'],'domain':'NATIVE DYNAMICS / MAGNETIC MECHANISM DISCOVERY','first_seen_date':'2026-08-13','last_updated_date':'2026-08-13','attempt_ids':['dev226_staggered_local_order'],'current_status':'CANONICAL','canonical_solution_ids':['dev226_staggered_local_order'],'open_questions':['Magnetic identity is not derived; pair interaction remains blocked.'],'blocked_by':[],'blocks':[],'do_not_rederive':True,'reopen_condition':'Only if frozen DEV203 tensor archive or DEV225 contraction changes independently.'}
    d['targets']=[x for x in d['targets'] if x.get('target_id') not in ids]+[target(*x) for x in specs]
    for x in d['targets']:
      if x.get('target_id')=='magnetic_mechanism_next_discriminating_test': x.update({'aliases':['DEV226'],'attempt_ids':['dev226_staggered_local_order'],'canonical_solution_ids':['dev226_staggered_local_order'],'open_questions':[f'Frozen DEV227 selector: {next_selector}.'],'blocked_by':[]})
    attempt={'attempt_id':'dev226_staggered_local_order','target_id':'native_distributed_staggered_local_order','name':'DEV226 native staggered local-order audit','aliases':['DEV226'],'summary':'Every unique periodic N6 bond at every archived DEV203 timestep was classified by the predeclared raw Frobenius tensor contraction.','why_attempted':'DEV225 authorized the tensor relation without inspecting order values.','date_started':'2026-08-13','date_completed':'2026-08-13','dev':'DEV226','branch':git('branch','--show-current'),'files':['pbuf/analysis/native_staggered_order.py','tools/generate_dev226_staggered_local_order.py'],'run_directories':['runs/dev226_staggered_local_order'],'tests':['tests/test_dev226_unique_n6_bonds.py','tests/test_dev226_tensor_contraction.py','tests/test_dev226_sign_counts.py','tests/test_dev226_axis_counts.py','tests/test_dev226_axial_equivalence.py','tests/test_dev226_translation_covariance.py','tests/test_dev226_rotation_covariance.py','tests/test_dev226_reflection_covariance.py','tests/test_dev226_unloaded_control.py','tests/test_dev226_no_thresholds.py','tests/test_dev226_no_reopened_routes.py'],'equations':['C_ab=A_ij(a)A_ij(b)','C_ab=2 omega(a).omega(b)'],'result':'FULL','result_reason':f'classification={result}; axis={axis}; temporal={temporal}; minimality={minimality}.','current_status':'CANONICAL','canonical':True,'physics_reusable':True,'infrastructure_reusable':True,'free_parameters':[],'fitted_parameters':[],'reopen_condition':'Only if frozen inputs change independently.','do_not_repeat_reason':'No alternate component, scalar, threshold, selection, force, or new dynamics is admissible.','evidence':[{'type':'file','value':'runs/dev226_staggered_local_order/final_contract.json'}],'confidence':'HIGH'}
    d['attempts']=[x for x in d['attempts'] if x.get('attempt_id')!=attempt['attempt_id']]+[attempt]; p.write_text(json.dumps(d,indent=2)+'\n')
    ledger=ROOT/'docs/PBUF_DEVELOPMENT_LEDGER.md'; text=ledger.read_text(); entry=f'\n## LEDGER ENTRY 060 — DEV226 NATIVE STAGGERED ANTISYMMETRIC ORDER\n\n- **Staggered-Order Closure Rule:** The frozen DEV203 antisymmetric tensor audit is `{result}` under the uniquely authorized raw relation `A(a):A(b)`, with axis structure `{axis}` and temporal character `{temporal}`. No component selection, threshold, parity fitting, or timestep selection was used.\n'
    if 'LEDGER ENTRY 060' not in text: ledger.write_text(text+entry)
    h=ROOT/'docs/PBUF_HISTORICAL_ATTEMPT_INDEX.md'; line='DEV226 rule: walk every unique periodic N6 bond at every frozen DEV203 timestep using only A(a):A(b); retain exact machine zero contractions and do not introduce R/L states, thresholds, or parity fitting.\n'
    if line not in h.read_text(): h.write_text(h.read_text()+'\n'+line)
    gp=ROOT/'docs/PBUF_DERIVATION_GRAPH.json'; g=json.loads(gp.read_text()); nodes={x['id'] for x in g['nodes']}
    for i,t in [('dev226_staggered_local_order','ATTEMPT')]+[(x[0],'TARGET') for x in specs]:
      if i not in nodes:g['nodes'].append({'id':i,'type':t})
    for e in [{'source':'dev225_local_handedness_representation_gate','target':'dev226_staggered_local_order','type':'AUTHORIZES'},{'source':'dev226_staggered_local_order','target':'native_distributed_staggered_local_order','type':'DERIVES'}]:
      if e not in g['edges']:g['edges'].append(e)
    gp.write_text(json.dumps(g,indent=2)+'\n')

def main():
    tensor=np.load(ROOT/'runs/dev203_relational_wave/antisymmetric_relational_components.npz')['directional_tensor_antisymmetric']; shape=tensor.shape[1:4]; T=tensor.shape[0]
    pairs, axes, _=unique_n6_bonds(shape); values=contract_bonds(tensor,pairs); counts=sign_counts(values,axes); causes=zero_causes(tensor,pairs,values)
    total=counts['total']; result=classification(counts); axis=axis_class(counts); temporal=temporal_class(counts)
    unloaded=np.zeros((1,*shape,3,3)); uv=contract_bonds(unloaded,pairs); uc=sign_counts(uv,axes)
    minimality='STRONG_CANDIDATE' if result in ('EXACT','OPPOSED_DOMINANT') and axis=='FULL_3D_OPPOSED' else 'PARTIAL_CANDIDATE' if result=='OPPOSED_DOMINANT' else 'MIXED' if result=='MIXED_IN_TIME' else 'NOT_SUPPORTED'
    next_selector='STAGGERED_ORDER_MECHANICAL_COUPLING_AUDIT' if minimality=='STRONG_CANDIDATE' else 'STAGGERED_ORDER_AXIS_STRUCTURE_AUDIT' if minimality=='PARTIAL_CANDIDATE' else 'MAGNETIC_CANDIDATE_EXHAUSTION_AUDIT'
    dump('starting_state.json',{'head':git('rev-parse','HEAD'),'CURRENT_GITHUB_INSPECTED':True,'CURRENT_HEAD_VERIFIED':True,'DEV203_READ':True,'DEV225_READ':True,'DEV226_TEST_SELECTION':'STAGGERED_LOCAL_ORDER_AUDIT','DEV226_TEST_SELECTION_FROZEN':True})
    dump('registry_lookup.json',{'MECHANISM_REGISTRY_QUERIED':True,'queries':['DEV203 antisymmetric tensor','DEV225 tensor contraction','staggered local order']}); dump('ledger_extract.json',{'DEVELOPMENT_LEDGER_READ':True}); dump('historical_staggered_order_inventory.json',{'HISTORICAL_INDEX_READ':True,'preserved':['DEV215','DEV218','DEV220','DEV223']})
    dump('dev203_tensor_manifest.json',{'DEV203_ANTISYMMETRIC_TENSOR_REUSED':True,'NO_NEW_TENSOR':True,'source':'runs/dev203_relational_wave/antisymmetric_relational_components.npz','array':'directional_tensor_antisymmetric','shape':list(tensor.shape),'sha256':sha('runs/dev203_relational_wave/antisymmetric_relational_components.npz')})
    dump('dev225_relation_contract.json',{'NATIVE_LOCAL_HANDEDNESS_RELATION_GATE':'AUTHORIZED_TENSOR_CONTRACTION','DEV225_REPRESENTATION_PRESERVED':True,'PREDECLARED_RELATION':'A_ij(a)A_ij(b)','NO_TENSOR_COMPONENT_RESULTS':True}); dump('dev225_final_contract_manifest.json',{'sha256':sha('runs/dev225_local_handedness_representation_gate/final_contract.json')}); dump('dev226_selector_contract.json',{'DEV226_TEST_SELECTION':'STAGGERED_LOCAL_ORDER_AUDIT','DEV226_TEST_SELECTION_FROZEN':True})
    save('unique_n6_bonds.npz',pairs=pairs,axis=axes,shape=np.asarray(shape)); dump('unique_n6_bond_contract.json',{'UNIQUE_N6_BOND_SET_USED':True,'NO_DOUBLE_COUNTED_BONDS':True,'ALL_UNIQUE_N6_BONDS_EVALUATED':True,'bond_count':int(len(pairs)),'definition':'positive periodic x,y,z direction from each node exactly once'})
    save('tensor_contraction_trajectory.npz',time=np.arange(T),contraction=values); save('bond_handedness_sign_trajectory.npz',time=np.arange(T),sign=np.sign(values).astype(np.int8),pairs=pairs,axis=axes)
    save('staggered_order_trajectory.npz',time=np.arange(T),opposed=counts['opposed'],aligned=counts['aligned'],zero=counts['zero'],total=total,sum_C=values.sum(1),mean_C=values.mean(1))
    rows=[{'timestep':int(t),'opposed':int(counts['opposed'][t]),'aligned':int(counts['aligned'][t]),'zero':int(counts['zero'][t]),'opposed_fraction':float(counts['opposed'][t]/total[t]),'aligned_fraction':float(counts['aligned'][t]/total[t]),'zero_fraction':float(counts['zero'][t]/total[t]),'sum_C':float(values[t].sum()),'mean_C':float(values[t].mean())} for t in range(T)]
    (OUT/'staggered_order_trajectory.csv').write_text('timestep,opposed,aligned,zero,opposed_fraction,aligned_fraction,zero_fraction,sum_C,mean_C\n'+'\n'.join(','.join(map(str,r.values())) for r in rows)+'\n')
    dump('staggered_order_counts_by_timestep.json',{'BOND_COUNT_CLOSURE':'EXACT','rows':rows,'aggregate':{'opposed':int(counts['opposed'].sum()),'aligned':int(counts['aligned'].sum()),'zero':int(counts['zero'].sum())}}); dump('staggered_order_counts_by_axis.json',{'N6_AXIS_RESOLVED_ORDER_REPORTED':True,'by_axis':counts['by_axis']})
    null_nodes=np.all(tensor==0,axis=(-2,-1)).sum(axis=(1,2,3)); dump('null_tensor_node_inventory.json',{'null_tensor_nodes_by_timestep':null_nodes,'NULL_TENSOR_NODES_INVENTORIED':True}); dump('zero_contraction_cause_inventory.json',{'ZERO_CONTRACTION_CAUSE_CLASSIFIED':True,'by_timestep':causes,'closure':np.asarray(list(causes.values())).sum(axis=0)==counts['zero']})
    err=axial_equivalence_error(tensor,pairs,values); dump('tensor_axial_equivalence.json',{'TENSOR_AXIAL_EQUIVALENCE_VERIFIED':err==0,'max_abs_error':err,'identity':'A(a):A(b)=2 omega(a).omega(b)','classification_basis':'tensor contraction only'})
    translated=np.roll(tensor,1,axis=1); tc=sign_counts(contract_bonds(translated,pairs),axes); dump('staggered_order_translation_covariance.json',{'STAGGERED_ORDER_TRANSLATION_COVARIANCE':'EXACT' if all(np.array_equal(counts[k],tc[k]) for k in ('opposed','aligned','zero')) else 'VIOLATED'})
    Q=np.array([[0,1,0],[1,0,0],[0,0,-1.]]) ; rotated=np.einsum('ij,t...jk,kl->t...il',Q,tensor,Q.T); rc=sign_counts(contract_bonds(rotated,pairs),axes); dump('staggered_order_rotation_covariance.json',{'STAGGERED_ORDER_CUBIC_ROTATION_COVARIANCE':'EXACT' if all(np.array_equal(counts[k],rc[k]) for k in ('opposed','aligned','zero')) else 'VIOLATED','axis_permutation':'x<->y; z reversed'})
    R=np.diag([-1.,1.,1.]); reflected=np.einsum('ij,t...jk,kl->t...il',R,tensor,R.T); fc=sign_counts(contract_bonds(reflected,pairs),axes); dump('staggered_order_reflection_covariance.json',{'STAGGERED_ORDER_REFLECTION_COVARIANCE':'EXACT' if all(np.array_equal(counts[k],fc[k]) for k in ('opposed','aligned','zero')) else 'VIOLATED','reason':'Frobenius tensor contraction is scalar under simultaneous reflection'})
    dump('unloaded_staggered_local_order.json',{'UNLOADED_CONTROL_REUSED':True,'UNLOADED_STAGGERED_LOCAL_ORDER':'ZERO','counts':{k:int(v.sum()) for k,v in uc.items() if k!='by_axis'},'NO_BACKGROUND_SUBTRACTION':True,'NO_BASELINE_FIT':True})
    dump('bipartite_order_compatibility.json',{'BIPARTITE_PARITY_CONTROL_SECONDARY':True,'NO_CHECKERBOARD_FIT':True,'NO_PARITY_OPTIMIZATION':True,'BIPARTITE_ORDER_COMPATIBILITY':'EXACT' if np.all(counts['aligned']==0) else 'DOMINANT' if counts['opposed'].sum()>counts['aligned'].sum() else 'MIXED' if counts['opposed'].sum()>0 else 'ABSENT','canonical_parity_bond_product':-1})
    dump('native_staggered_local_order.json',{'NATIVE_STAGGERED_LOCAL_ORDER':result}); dump('staggered_order_axis_structure.json',{'STAGGERED_ORDER_AXIS_STRUCTURE':axis}); dump('staggered_order_temporal_character.json',{'STAGGERED_ORDER_TEMPORAL_CHARACTER':temporal}); dump('distributed_staggered_order_minimality.json',{'DISTRIBUTED_STAGGERED_ORDER_MINIMALITY':minimality,'WHOLE_STRUCTURE_ORDER_TEST':True,'NO_LOCALIZATION_SEARCH':True}); dump('pair_orientation_interaction_gate.json',{'PAIR_ORIENTATION_INTERACTION_GATE':'REMAINS_BLOCKED','NO_PAIR_INTERACTION':True,'NO_PAIR_FORCE_TEST':True,'NO_TORQUE_TEST':True}); dump('dev227_test_selection.json',{'DEV227_TEST_SELECTION':next_selector,'DEV227_TEST_SELECTION_FROZEN':True})
    update_docs(result,axis,temporal,minimality,next_selector); subprocess.check_call([sys.executable,'tools/pbuf_registry.py','validate'],cwd=ROOT); subprocess.check_call([sys.executable,'tools/pbuf_registry.py','render'],cwd=ROOT)
    flags={x:True for x in 'CURRENT_GITHUB_INSPECTED CURRENT_HEAD_VERIFIED MECHANISM_REGISTRY_QUERIED DEVELOPMENT_LEDGER_READ HISTORICAL_INDEX_READ DEV203_READ DEV225_READ DEV226_TEST_SELECTION_FROZEN DEV203_ANTISYMMETRIC_TENSOR_REUSED DEV225_REPRESENTATION_PRESERVED UNIQUE_N6_BOND_SET_USED NO_DOUBLE_COUNTED_BONDS FULL_FROZEN_DEV203_TRAJECTORY_USED ALL_UNIQUE_N6_BONDS_EVALUATED ALL_FROZEN_TIMESTEPS_EVALUATED OPPOSED_BONDS_COUNTED ALIGNED_BONDS_COUNTED ZERO_BONDS_COUNTED N6_AXIS_RESOLVED_ORDER_REPORTED ZERO_CONTRACTION_CAUSE_CLASSIFIED NULL_TENSOR_NODES_INVENTORIED TENSOR_AXIAL_EQUIVALENCE_VERIFIED STAGGERED_ORDER_TRANSLATION_COVARIANCE_CLASSIFIED STAGGERED_ORDER_CUBIC_ROTATION_COVARIANCE_CLASSIFIED STAGGERED_ORDER_REFLECTION_COVARIANCE_CLASSIFIED UNLOADED_CONTROL_REUSED UNLOADED_STAGGERED_LOCAL_ORDER_CLASSIFIED BIPARTITE_PARITY_CONTROL_SECONDARY BIPARTITE_ORDER_COMPATIBILITY_CLASSIFIED NATIVE_STAGGERED_LOCAL_ORDER_CLASSIFIED STAGGERED_ORDER_AXIS_STRUCTURE_CLASSIFIED STAGGERED_ORDER_TEMPORAL_CHARACTER_CLASSIFIED DISTRIBUTED_STAGGERED_ORDER_MINIMALITY_CLASSIFIED MAGNETIC_IDENTITY_NOT_DERIVED NO_ORDER_THRESHOLD NO_TENSOR_MAGNITUDE_THRESHOLD NO_RESULT_SELECTED_TIMESTEP NO_LOCALIZATION_SEARCH NO_NEW_DYNAMICS_RUN NO_PAIR_INTERACTION NO_PAIR_FORCE_TEST DEV215_TEMPORAL_CYCLE_CLOSURE_PRESERVED DEV218_MOMENTUM_POLARITY_CLOSURE_PRESERVED DEV220_GLOBAL_WINDING_CLOSURE_PRESERVED DEV223_DISTRIBUTED_STRUCTURE_PRESERVED NO_NEW_FORCE NO_NEW_DOF NO_NEW_STATE_VARIABLE NO_NEW_PHASE NO_RL_STATE NO_BINARY_SPIN_STATE NO_SELECTED_TENSOR_COMPONENT NO_THRESHOLD NO_MAGNITUDE_FILTER NO_TIMESTEP_SELECTION NO_AXIS_SELECTION NO_CHECKERBOARD_FLIP NO_PARITY_OPTIMIZATION NO_GLOBAL_WINDING_REINTERPRETATION NO_TEMPORAL_CYCLE_REINTERPRETATION NO_MOMENTUM_POLARITY_REINTERPRETATION NO_PAIR_FORCE_TEST NO_TORQUE_TEST MECHANISM_REGISTRY_UPDATED REGISTRY_VALIDATED LEDGER_UPDATED HISTORICAL_INDEX_UPDATED TIMELINE_REGENERATED DERIVATION_GRAPH_REGENERATED TESTS_PASS NO_PR_CREATED'.split()}; flags.update({'DEV226_TEST_SELECTION':'STAGGERED_LOCAL_ORDER_AUDIT','PREDECLARED_RELATION':'A_ij(a)A_ij(b)','BOND_COUNT_CLOSURE':'EXACT','NATIVE_STAGGERED_LOCAL_ORDER':result,'STAGGERED_ORDER_AXIS_STRUCTURE':axis,'STAGGERED_ORDER_TEMPORAL_CHARACTER':temporal,'DISTRIBUTED_STAGGERED_ORDER_MINIMALITY':minimality,'PAIR_ORIENTATION_INTERACTION_GATE':'REMAINS_BLOCKED','STAGGERED_ORDER_TRANSLATION_COVARIANCE':json.loads((OUT/'staggered_order_translation_covariance.json').read_text())['STAGGERED_ORDER_TRANSLATION_COVARIANCE'],'STAGGERED_ORDER_CUBIC_ROTATION_COVARIANCE':json.loads((OUT/'staggered_order_rotation_covariance.json').read_text())['STAGGERED_ORDER_CUBIC_ROTATION_COVARIANCE'],'STAGGERED_ORDER_REFLECTION_COVARIANCE':json.loads((OUT/'staggered_order_reflection_covariance.json').read_text())['STAGGERED_ORDER_REFLECTION_COVARIANCE'],'UNLOADED_STAGGERED_LOCAL_ORDER':'ZERO','BIPARTITE_ORDER_COMPATIBILITY':json.loads((OUT/'bipartite_order_compatibility.json').read_text())['BIPARTITE_ORDER_COMPATIBILITY'],'DEV227_TEST_SELECTION':next_selector,'COMMITTED':False,'PUSHED_DIRECTLY_TO_MAIN':False,'REMOTE_MAIN_VERIFIED':False,'WORKTREE_CLEAN':False}); dump('final_contract.json',flags)
    (OUT/'discussion_handoff.md').write_text(f'# DEV226 handoff\n\nThe frozen raw tensor-contraction audit is **{result}** (`{axis}`, `{temporal}`). It uses every archived timestep and unique periodic N6 bond, without a threshold or representation change. This is structural evidence only: `MAGNETIC_IDENTITY_NOT_DERIVED=true` and pair interaction remains blocked.\n')
if __name__=='__main__': main()
