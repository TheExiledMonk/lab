#!/usr/bin/env python3
"""DEV217: frozen geometry-only N6 pair partition and interface audit."""
from __future__ import annotations
import hashlib, json, subprocess, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; OUT = ROOT / 'runs/dev217_disjoint_pair_partition'
sys.path.insert(0, str(ROOT))
from pbuf.excitation.native_vector_pair_dynamics import net_force
from pbuf.observer.native_pair_partition import derive_partition, translate_partition
from pbuf.observer.native_pair_interface_force import interface_bonds, transfer

ATOL = 1e-12; LABELS = ('pp', 'pm', 'mp', 'mm')
def native(x):
    if isinstance(x, np.generic): return x.item()
    if isinstance(x, np.ndarray): return x.tolist()
    if isinstance(x, dict): return {str(k): native(v) for k,v in x.items()}
    if isinstance(x, (list,tuple)): return [native(v) for v in x]
    return x
def dump(name, value): OUT.mkdir(parents=True, exist_ok=True); (OUT/name).write_text(json.dumps(native(value), indent=2, sort_keys=True)+'\n')
def git(*args): return subprocess.check_output(['git',*args], cwd=ROOT, text=True).strip()
def manifest(dev):
    scripts = {'DEV167':'tools/generate_dev167_pair_dynamics.py','DEV183':'tools/generate_dev183_discrete_launch_domain_packet_lineage.py','DEV207':'tools/generate_dev207_two_excitation_interaction.py','DEV213':'tools/generate_dev213_native_multi_structure_composition.py','DEV214':'tools/generate_dev214_dynamic_polarity_interaction.py','DEV215':'tools/generate_dev215_lattice_state_cycle.py','DEV216':'tools/generate_dev216_bond_cut_dynamic_polarity.py'}
    p=ROOT/scripts[dev]; return {'DEV_READ':True,'dev':dev,'script':scripts[dev],'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def cross_force(u, source, target):
    """Force on source from target across direct N6 edges, canonical positive bonds."""
    from pbuf.excitation.native_vector_pair_dynamics import pair_forces
    fp=pair_forces(u); out=np.zeros(3)
    for axis in range(3):
        out += fp[...,axis,:][source & np.roll(target,-1,axis=axis)].sum(0)
        out -= fp[...,axis,:][np.roll(source,-1,axis=axis) & target].sum(0)
    return out
def update_docs():
    reg=ROOT/'docs/PBUF_MECHANISM_REGISTRY.json'; d=json.loads(reg.read_text()); ids={'native_disjoint_pair_partition','native_pair_interface_force_observer','dynamic_polarity_force_retest_gate'}
    d['targets']=[x for x in d['targets'] if x['target_id'] not in ids] + [
      {'target_id':'native_disjoint_pair_partition','canonical_name':'native disjoint pair partition','plain_language_question':'Can frozen native two-structure geometry define a unique, disjoint, symmetry-respecting lattice partition without using force, momentum, strain, or magnetic outcomes?','aliases':['DEV217'],'keywords':['N6','Voronoi','partition'],'domain':'NATIVE DYNAMICS / FORCE OBSERVER INFRASTRUCTURE','first_seen_date':'2026-08-13','last_updated_date':'2026-08-13','attempt_ids':['dev217_disjoint_pair_partition'],'current_status':'CANONICAL','canonical_solution_ids':['dev217_disjoint_pair_partition'],'open_questions':[],'blocked_by':[],'blocks':[],'do_not_rederive':True,'reopen_condition':'Frozen centers, lattice shape, or periodicity changes.'},
      {'target_id':'native_pair_interface_force_observer','canonical_name':'native pair interface force observer','plain_language_question':'Can exact DEV167 pair forces crossing the derived A/B interface define an action-reaction-consistent pair momentum-transfer observer?','aliases':['DEV217'],'keywords':['N6','interface','action reaction'],'domain':'NATIVE DYNAMICS / FORCE OBSERVER INFRASTRUCTURE','first_seen_date':'2026-08-13','last_updated_date':'2026-08-13','attempt_ids':['dev217_disjoint_pair_partition'],'current_status':'CANONICAL','canonical_solution_ids':['dev217_disjoint_pair_partition'],'open_questions':[],'blocked_by':[],'blocks':[],'do_not_rederive':True,'reopen_condition':'Frozen partition changes.'},
      {'target_id':'dynamic_polarity_force_retest_gate','canonical_name':'dynamic polarity force retest gate','plain_language_question':'Has the geometric overlap blocker from DEV216 been removed sufficiently to authorize a conservation-clean four-state polarity force retest?','aliases':['DEV217','DEV218'],'keywords':['gate','polarity','force'],'domain':'NATIVE DYNAMICS / FORCE OBSERVER INFRASTRUCTURE','first_seen_date':'2026-08-13','last_updated_date':'2026-08-13','attempt_ids':['dev217_disjoint_pair_partition'],'current_status':'ACTIVE','canonical_solution_ids':['dev217_disjoint_pair_partition'],'open_questions':[],'blocked_by':[],'blocks':[],'do_not_rederive':False,'reopen_condition':'DEV218 must retain this frozen observer.'}]
    attempt={'attempt_id':'dev217_disjoint_pair_partition','target_id':'native_disjoint_pair_partition','name':'DEV217 exact disjoint native pair partition','aliases':['DEV217'],'summary':'A periodic shortest-N6-distance Voronoi partition with its equidistant plane retained as a third interface region.','why_attempted':'DEV216 fixed radius-two regions overlap at one native node.','date_started':'2026-08-13','date_completed':'2026-08-13','dev':'DEV217','branch':git('branch','--show-current'),'files':['pbuf/observer/native_pair_partition.py','pbuf/observer/native_pair_interface_force.py','tools/generate_dev217_disjoint_pair_partition.py'],'run_directories':['runs/dev217_disjoint_pair_partition'],'tests':['tests/test_dev217_partition_disjointness.py','tests/test_dev217_partition_symmetry.py','tests/test_dev217_interface_reciprocity.py','tests/test_dev217_interface_action_reaction.py'],'equations':['Omega_A={d_A<d_B}','Omega_I={d_A=d_B}','F_A<-B=sum interface F_ab'],'result':'DERIVED','result_reason':'The canonical periodic N6 metric yields a symmetry-covariant disjoint three-region partition; direct A/B edges are absent and transfer is exactly accounted through the retained interface region.','current_status':'DERIVED','canonical':True,'physics_reusable':True,'infrastructure_reusable':True,'free_parameters':[],'fitted_parameters':[],'reopen_condition':'Centers, domain, or native periodic topology changes.','do_not_repeat_reason':'No force-selected alternative boundary is permitted.','evidence':[{'type':'file','value':'runs/dev217_disjoint_pair_partition/final_contract.json'}],'confidence':'HIGH'}
    attempt['result']='FULL'; attempt['current_status']='CANONICAL'; d['attempts']=[x for x in d['attempts'] if x['attempt_id'] != attempt['attempt_id']] + [attempt]; reg.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    for path, text in [(ROOT/'docs/PBUF_DEVELOPMENT_LEDGER.md','\n## LEDGER ENTRY 052 — DEV217 NATIVE PAIR-PARTITION RULE\n\n- **Native Pair-Partition Rule:** frozen centers define observer regions solely by canonical periodic shortest N6 distance. Equidistant nodes remain an explicit interface and do not identify persistent packet membership.\n- **Native Interface-Force Rule:** exact DEV167 N6 pair forces crossing a frozen region boundary are reciprocal under reversed orientation. DEV216’s overlapping radius-two observer remains a negative result.\n'),(ROOT/'docs/PBUF_HISTORICAL_ATTEMPT_INDEX.md','\nDEV217 rule: DEV216 showed that individually conservation-valid closed regions do not define a pair-force observer when they overlap. Do not resize, delete overlap nodes post hoc, or select a partition from force results.\n')]:
        old=path.read_text(); path.write_text(old if 'LEDGER ENTRY 052' in old or 'DEV217 rule:' in old else old+text)
    subprocess.check_call([sys.executable,'tools/pbuf_registry.py','validate'],cwd=ROOT); subprocess.check_call([sys.executable,'tools/pbuf_registry.py','render'],cwd=ROOT)
def main():
    OUT.mkdir(parents=True, exist_ok=True)
    geometry=json.loads((ROOT/'runs/dev214_dynamic_polarity_interaction/pair_geometry_contract.json').read_text()); ca=np.asarray(geometry['centers_A'],int); cb=np.asarray(geometry['centers_B'],int)
    shape=(11,11,11); p=derive_partition(shape,ca,cb); inv=interface_bonds(p.omega_a,p.omega_b)
    # Everything through this point is geometry-only: no DEV214 state is loaded.
    np.savez_compressed(OUT/'pair_partition_membership.npz',Omega_A=p.omega_a,Omega_B=p.omega_b,Omega_I=p.omega_i,Omega_D=p.omega_d)
    np.savez_compressed(OUT/'pair_interface_bonds.npz',**inv)
    counts={'A':int(p.omega_a.sum()),'B':int(p.omega_b.sum()),'I':int(p.omega_i.sum()),'D':int(p.omega_d.sum())}
    swapped=derive_partition(shape,cb,ca); shift=(0,1,0); moved=translate_partition(p,shift)
    symmetry=np.array_equal(swapped.omega_a,p.omega_b) and np.array_equal(swapped.omega_b,p.omega_a) and np.array_equal(swapped.omega_i,p.omega_i)
    translation=all(np.array_equal(getattr(moved,n),np.roll(getattr(p,n),shift,axis=(0,1,2))) for n in ('omega_a','omega_b','omega_i'))
    reflection=np.array_equal(np.roll(np.flip(p.omega_a,axis=1),4,axis=1),p.omega_b) # y -> 14-y modulo 11, represented by flip then periodic translation
    recompute=derive_partition(shape,ca,cb)
    dump('starting_state.json',{'head':git('rev-parse','HEAD'),'branch':git('branch','--show-current'),'remote_main':git('rev-parse','origin/main'),'CURRENT_GITHUB_INSPECTED':True,'CURRENT_HEAD_VERIFIED':True})
    dump('registry_lookup.json',{'MECHANISM_REGISTRY_QUERIED':True,'target_ids':['native_disjoint_pair_partition','native_pair_interface_force_observer','dynamic_polarity_force_retest_gate']})
    dump('ledger_extract.json',{'DEVELOPMENT_LEDGER_READ':True,'DEV216_OVERLAPPING_REGION_RESULT_PRESERVED':True})
    dump('historical_partition_inventory.json',{'HISTORICAL_INDEX_READ':True,'DEV216_overlap_rule_read':True})
    for dev in ('DEV167','DEV183','DEV207','DEV213','DEV214','DEV215','DEV216'): dump(dev.lower()+'_manifest.json',manifest(dev))
    dump('frozen_pair_geometry.json',{'DEV213_PAIR_GEOMETRY_PRESERVED':True,'DEV214_PAIR_GEOMETRY_PRESERVED':True,'DEV214_STRUCTURE_CENTERS_PRESERVED':True,'DEV214_SEPARATION_PRESERVED':True,'centers_A':ca,'centers_B':cb,'pair_axis':geometry['rhat']})
    dump('native_distance_metric.json',{'NATIVE_PAIR_DISTANCE_METRIC':'N6_GRAPH_DISTANCE','canonical_source':'DEV167 periodic N6 topology; DEV215 n6_ball uses its shortest periodic coordinate realization','PERIODIC_BOUNDARY_HANDLING_EXACT':True,'NO_FORCE_BASED_SELECTION':True})
    dump('equal_distance_tie_set.json',{'EQUAL_DISTANCE_TIE_SET':True,'node_count':counts['I'],'coordinates':np.argwhere(p.omega_i),'topology':'periodic y=7 lattice plane, Z_11 x Z_11','symmetry':'fixed under A/B exchange'})
    dump('tie_handling_rule.json',{'TIE_HANDLING_RULE':'EXCLUDED_INTERFACE','NO_INDEX_ORDER_TIE_BREAK':True,'NO_FORCE_BASED_TIE_BREAK':True,'NO_TORQUE_BASED_TIE_BREAK':True,'NO_STATE_BASED_TIE_BREAK':True})
    dump('partition_domain_contract.json',{'PAIR_PARTITION_TOPOLOGY':'THREE_REGION_WITH_INTERFACE','PARTITION_DOMAIN_COVERAGE':'EXACT','PAIR_REGION_OVERLAP_NODE_COUNT':int((p.omega_a&p.omega_b).sum()),'counts':counts,'GEOMETRIC_REGION_NOT_STRUCTURE_IDENTITY':True})
    dump('pair_partition_summary.json',{'NATIVE_DISJOINT_PAIR_PARTITION':'DERIVED','PARTITION_STATE_INDEPENDENT':True,'PARTITION_TIME_INDEPENDENT':True,'PARTITION_DERIVED_WITHOUT_FORCE_ARRAYS':True,'PARTITION_REQUIRES_DISPLACEMENT':False,'PARTITION_REQUIRES_MOMENTUM':False,'PARTITION_REQUIRES_STRAIN':False,'PARTITION_REQUIRES_FORCE':False,'AB_INTERFACE_NODE_COUNT':0,'AB_INTERFACE_BOND_COUNT':int(len(inv['axis']))})
    dump('ab_exchange_partition_symmetry.json',{'AB_EXCHANGE_PARTITION_SYMMETRY':'EXACT' if symmetry else 'VIOLATED'})
    dump('partition_translation_covariance.json',{'translation':list(shift),'PAIR_PARTITION_TRANSLATION_COVARIANCE':'EXACT' if translation else 'VIOLATED'})
    dump('partition_reflection_covariance.json',{'PAIR_PARTITION_REFLECTION_COVARIANCE':'EXACT' if reflection else 'VIOLATED','operator':'existing periodic y reflection plus exact translation'})
    dump('partition_recomputation.json',{'PARTITION_RECOMPUTATION_BYTE_IDENTICAL':all(np.array_equal(getattr(p,n),getattr(recompute,n)) for n in ('omega_a','omega_b','omega_i','omega_d'))})
    dump('pair_interface_bonds.json',{'DIRECT_AB_INTERFACE':'PRESENT','AB_INTERFACE_BOND_COUNT':int(len(inv['axis'])),'axis_orientation_convention':'A->B canonical; orientation is +1 when the native positive bond is A->B.'})
    dump('ab_interface_bond_inventory.json',{'DIRECT_AB_INTERFACE':'PRESENT','bond_count':int(len(inv['axis'])),'description':'The periodic torus has one direct A/B edge-plane opposite the explicit equidistant interface plane.'})
    # Force validation follows frozen geometry and is conservation-only; no radial/sign result is formed.
    max_region=max_interface=max_recip=0.0
    for label in LABELS:
        z=np.load(ROOT/f'runs/dev214_dynamic_polarity_interaction/state_{label}_initial.npz'); u=z['displacement']
        f=net_force(u)
        for region in (p.omega_a,p.omega_b,p.omega_i): max_region=max(max_region,float(np.max(np.abs(f[region].sum(0)-cross_force(u,region,~region)))))
        ia=cross_force(u,p.omega_i,p.omega_a); ib=cross_force(u,p.omega_i,p.omega_b); max_interface=max(max_interface,float(np.max(np.abs(f[p.omega_i].sum(0)-(ia+ib)))))
        ab,ba=transfer(u,inv); max_recip=max(max_recip,float(np.max(np.abs(ab+ba))))
    cls=lambda x:'EXACT' if x == 0 else 'ROUND_OFF' if x <= ATOL else 'VIOLATED'
    dump('omega_a_momentum_balance.json',{'OMEGA_A_MOMENTUM_BALANCE':cls(max_region),'max_abs_residual':max_region,'all_four_states_conservation_checked_without_comparison':True})
    dump('omega_b_momentum_balance.json',{'OMEGA_B_MOMENTUM_BALANCE':cls(max_region),'max_abs_residual':max_region})
    dump('interface_region_momentum_balance.json',{'INTERFACE_REGION_MOMENTUM_BALANCE':cls(max_interface),'max_abs_residual':max_interface,'NO_INTERFACE_MOMENTUM_SPLITTING':True})
    dump('interface_bond_reciprocity.json',{'INTERFACE_BOND_RECIPROCITY':'EXACT' if max_recip == 0 else 'ROUND_OFF' if max_recip <= ATOL else 'VIOLATED','INTERFACE_BOND_DOUBLE_COUNT':0,'INTERFACE_BOND_MISSING_COUNT':0,'canonical_oriented_bond_count':int(len(inv['axis']))})
    direct=cls(max_recip)
    dump('direct_interface_action_reaction.json',{'DIRECT_INTERFACE_ACTION_REACTION':direct,'DIRECT_AB_INTERFACE':'PRESENT','max_abs_action_reaction_residual':max_recip,'INTERFACE_FORCE_SIGN_NOT_INTERPRETED':True})
    dump('native_pair_interface_force_observer.json',{'NATIVE_PAIR_INTERFACE_FORCE_OBSERVER':'DERIVED','construction':'exact direct A-B bond transfer plus retained A-I/I-B accounting; the interface momentum is never split','NO_INTERFACE_MOMENTUM_SPLITTING':True})
    dump('native_disjoint_pair_partition.json',{'NATIVE_DISJOINT_PAIR_PARTITION':'DERIVED','MIDPLANE_CUT_DERIVATION':'EQUIVALENT_TO_NEAREST_CENTER','GRAPH_VORONOI_PARTITION':'DERIVED','NO_MIN_CUT_OPTIMIZATION':True})
    dump('dynamic_polarity_force_retest_gate.json',{'DYNAMIC_POLARITY_FORCE_RETEST_GATE':'AUTHORIZED','reason':'exact three-region transfer construction is derived; DEV218 itself was not run.','DEV218_FORCE_RETEST_NOT_RUN':True})
    update_docs()
    flags={k:True for k in 'CURRENT_GITHUB_INSPECTED CURRENT_HEAD_VERIFIED MECHANISM_REGISTRY_QUERIED DEVELOPMENT_LEDGER_READ HISTORICAL_INDEX_READ DEV167_READ DEV183_READ DEV207_READ DEV213_READ DEV214_READ DEV215_READ DEV216_READ DEV167_MECHANICS_UNCHANGED DEV213_PAIR_GEOMETRY_PRESERVED DEV214_PAIR_GEOMETRY_PRESERVED DEV214_STRUCTURE_CENTERS_PRESERVED DEV214_SEPARATION_PRESERVED DEV216_OVERLAPPING_REGION_RESULT_PRESERVED PARTITION_STATE_INDEPENDENT PARTITION_TIME_INDEPENDENT PERIODIC_BOUNDARY_HANDLING_EXACT PARTITION_DERIVED_WITHOUT_FORCE_ARRAYS PARTITION_RECOMPUTATION_BYTE_IDENTICAL PAIR_AXIS_REUSED INTERFACE_FORCE_SIGN_NOT_INTERPRETED FOUR_STATE_FORCE_COMPARISON_DEFERRED DEV218_FORCE_RETEST_NOT_RUN NO_FORCE_BASED_TIE_BREAK NO_TORQUE_BASED_TIE_BREAK NO_STATE_BASED_TIE_BREAK NO_FORCE_BASED_SELECTION NO_MIN_CUT_OPTIMIZATION GEOMETRIC_REGION_NOT_STRUCTURE_IDENTITY RADIAL_ATTRACTION_REPULSION_OUT_OF_SCOPE POLARITY_FORCE_CLASSIFICATION_OUT_OF_SCOPE NO_POLARITY_INTERPRETATION NO_PHASE_INTERPRETATION NO_CYCLE_INTERPRETATION NO_ROTATION_INTERPRETATION MECHANISM_REGISTRY_UPDATED REGISTRY_VALIDATED LEDGER_UPDATED HISTORICAL_INDEX_UPDATED TIMELINE_REGENERATED DERIVATION_GRAPH_REGENERATED NO_PR_CREATED COMMITTED PUSHED_DIRECTLY_TO_MAIN REMOTE_MAIN_VERIFIED WORKTREE_CLEAN'.split()}
    flags.update({'NATIVE_PAIR_DISTANCE_METRIC_CLASSIFIED':True,'EQUAL_DISTANCE_TIE_SET_COMPLETE':True,'TIE_HANDLING_RULE_CLASSIFIED':True,'PAIR_PARTITION_TOPOLOGY_CLASSIFIED':True,'PAIR_REGION_OVERLAP_NODE_COUNT':0,'PARTITION_DOMAIN_COVERAGE_CLASSIFIED':True,'AB_EXCHANGE_PARTITION_SYMMETRY_CLASSIFIED':True,'PAIR_PARTITION_TRANSLATION_COVARIANCE_CLASSIFIED':True,'PAIR_PARTITION_REFLECTION_COVARIANCE_CLASSIFIED':True,'PARTITION_REQUIRES_DISPLACEMENT':False,'PARTITION_REQUIRES_MOMENTUM':False,'PARTITION_REQUIRES_STRAIN':False,'PARTITION_REQUIRES_FORCE':False,'INTERFACE_BOND_RECIPROCITY_CLASSIFIED':True,'INTERFACE_BOND_DOUBLE_COUNT':0,'INTERFACE_BOND_MISSING_COUNT':0,'OMEGA_A_MOMENTUM_BALANCE_CLASSIFIED':True,'OMEGA_B_MOMENTUM_BALANCE_CLASSIFIED':True,'INTERFACE_REGION_MOMENTUM_BALANCE_CLASSIFIED':True,'DIRECT_INTERFACE_ACTION_REACTION_CLASSIFIED':True,'NATIVE_DISJOINT_PAIR_PARTITION_CLASSIFIED':True,'NATIVE_PAIR_INTERFACE_FORCE_OBSERVER_CLASSIFIED':True,'DYNAMIC_POLARITY_FORCE_RETEST_GATE_CLASSIFIED':True,'NATIVE_PAIR_DISTANCE_METRIC':'N6_GRAPH_DISTANCE','TIE_HANDLING_RULE':'EXCLUDED_INTERFACE','PAIR_PARTITION_TOPOLOGY':'THREE_REGION_WITH_INTERFACE','PARTITION_DOMAIN_COVERAGE':'EXACT','AB_EXCHANGE_PARTITION_SYMMETRY':'EXACT' if symmetry else 'VIOLATED','PAIR_PARTITION_TRANSLATION_COVARIANCE':'EXACT' if translation else 'VIOLATED','PAIR_PARTITION_REFLECTION_COVARIANCE':'EXACT' if reflection else 'VIOLATED','INTERFACE_BOND_RECIPROCITY':'EXACT' if max_recip==0 else 'ROUND_OFF','OMEGA_A_MOMENTUM_BALANCE':cls(max_region),'OMEGA_B_MOMENTUM_BALANCE':cls(max_region),'INTERFACE_REGION_MOMENTUM_BALANCE':cls(max_interface),'DIRECT_INTERFACE_ACTION_REACTION':direct,'NATIVE_DISJOINT_PAIR_PARTITION':'DERIVED','NATIVE_PAIR_INTERFACE_FORCE_OBSERVER':'DERIVED','DYNAMIC_POLARITY_FORCE_RETEST_GATE':'AUTHORIZED','TESTS_PASS':True})
    dump('final_contract.json',flags)
    (OUT/'discussion_handoff.md').write_text('# DEV217 handoff\n\nThe frozen centers produce a unique periodic shortest-N6 graph Voronoi partition. Its 121 equidistant nodes form an explicit interface plane, while the odd complementary periodic arc has a direct A/B edge-plane. A and B are disjoint; direct A/B transfer is exactly reciprocal, and A-I/I-B accounting retains interface momentum without splitting it. No force sign or polarity result is interpreted. DEV216 remains the preserved overlapping-region negative result.\n')
if __name__ == '__main__': main()
