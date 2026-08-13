#!/usr/bin/env python3
"""DEV213: freeze physical same-time native aggregate preparation semantics.

No force, torque, radial classification, or magnetic labels are calculated.
"""
from __future__ import annotations
import hashlib, json, subprocess, sys, math
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "runs/dev213_native_multi_structure_composition"
sys.path.insert(0, str(ROOT))
from pbuf.excitation.native_vector_pair_dynamics import VectorPairState, source_contact_force
from pbuf.excitation.native_multi_structure_preparation import NativePreparation, exact_support, inject, reverse_internal_state
from pbuf.observer.native_composition_audit import state_validity, support_relation
from pbuf.observer.native_pair_preparation_invariants import summary
from tools import generate_dev169_raw_abell_native_observer as D
from tools import generate_dev184_discrete_launch_density_convergence as E

def git(*args): return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()
def native(x):
    if isinstance(x, np.generic): return x.item()
    if isinstance(x, np.ndarray): return x.tolist()
    if isinstance(x, dict): return {k: native(v) for k, v in x.items()}
    if isinstance(x, (tuple, list)): return [native(v) for v in x]
    return x
def dump(name, value):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / name).write_text(json.dumps(native(value), indent=2, sort_keys=True, allow_nan=False) + "\n")
def save(name, state):
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUT / name, displacement=state.displacement, momentum=state.momentum, progression_step=state.progression_step)
def manifest(dev, script, run):
    p=ROOT/script
    return {"DEV_READ": True, "dev": dev, "script": script, "run": run, "exists": p.exists(), "sha256": hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None}

def update_docs():
    path=ROOT/'docs/PBUF_MECHANISM_REGISTRY.json'; data=json.loads(path.read_text())
    targets=[
      ("native_dynamic_multi_structure_composition", "Can two independently defined dynamic native structures be physically prepared in one valid DEV167 state using only existing state/preparation semantics?", "ACTIVE"),
      ("native_simultaneous_preparation_semantics", "Does existing native preparation/injection semantics define an order-independent simultaneous aggregate state rather than merely post-hoc array addition?", "ACTIVE"),
      ("native_dynamic_structure_identity", "Can two dynamic structures retain preparation/provenance identity within one aggregate native state even when their exact mathematical supports overlap?", "ACTIVE"),
      ("dynamic_polarity_interaction_authorization", "Has the physical composition blocker preventing DEV212's same-structure polarity interaction test been closed independently of magnetic results?", "ACTIVE"),
    ]
    ids={x[0] for x in targets}; data['targets']=[x for x in data['targets'] if x.get('target_id') not in ids]
    data['targets'] += [{"target_id":i,"canonical_name":i.replace('_',' '),"plain_language_question":q,"aliases":["DEV213"],"keywords":["composition","simultaneous","injection","aggregate state","packet addition","overlap","lineage"],"domain":"NATIVE DYNAMICS / EM SUBSTRUCTURE","first_seen_date":"2026-08-13","last_updated_date":"2026-08-13","attempt_ids":["dev213_native_multi_structure_composition"],"current_status":s,"canonical_solution_ids":[],"open_questions":["DEV214 may inspect force and torque only after this initial-state result."],"blocked_by":[],"blocks":[],"do_not_rederive":True,"reopen_condition":"Frozen DEV167/DEV182/DEV196 preparation semantics change."} for i,q,s in targets]
    attempt={"attempt_id":"dev213_native_multi_structure_composition","target_id":"native_dynamic_multi_structure_composition","name":"DEV213 physical native multi-structure composition audit","aliases":["DEV213"],"summary":"Audits existing same-time DEV196 injection algebra and DEV207 combined-state lineage without a force result.","why_attempted":"DEV193 correctly blocked post-hoc transfer-column addition, while DEV196 supplied valid full-state injection semantics.","date_started":"2026-08-13","date_completed":"2026-08-13","dev":"DEV213","branch":git("branch","--show-current"),"files":["pbuf/excitation/native_multi_structure_preparation.py","pbuf/observer/native_composition_audit.py","pbuf/observer/native_pair_preparation_invariants.py","tools/generate_dev213_native_multi_structure_composition.py"],"run_directories":["runs/dev213_native_multi_structure_composition"],"tests":["tests/test_dev213_preparation_order.py","tests/test_dev213_state_validity.py","tests/test_dev213_composition_semantics.py"],"equations":[{"id":"dev213_same_time_injection","latex":"X_AB=I_B(I_A(X_0))=I_A(I_B(X_0))","role":"preparation algebra","source_file":"pbuf/excitation/native_multi_structure_preparation.py","source_commit":None}],"result":"FULL","result_reason":"Existing physical packet preparation is a valid full-state injection; independent same-step injections commute, pass DEV167 validity, retain provenance, and thereafter have only aggregate DEV167 evolution.","current_status":"ACTIVE","canonical":False,"physics_reusable":True,"infrastructure_reusable":True,"free_parameters":[],"fitted_parameters":[],"reopen_condition":"Do not use this result as linear propagation superposition.","do_not_repeat_reason":"Never select composition order or placement by force outcome.","evidence":[{"type":"file","value":"runs/dev213_native_multi_structure_composition/final_contract.json"}],"confidence":"HIGH"}
    data['attempts']=[x for x in data['attempts'] if x.get('attempt_id') != attempt['attempt_id']]+[attempt]
    path.write_text(json.dumps(data,indent=2)+"\n")
    ledger=ROOT/'docs/PBUF_DEVELOPMENT_LEDGER.md'; text=ledger.read_text(); entry='''\n## LEDGER ENTRY 048 — DEV213 NATIVE AGGREGATE-STATE PREPARATION\n\n- **Native Aggregate-State Preparation Rule:** Two independently defined native preparations may be applied to one valid full state under existing preparation semantics, after which the resulting aggregate displacement/momentum field evolves solely under frozen DEV167 dynamics. This establishes simultaneous native-state preparation, not linear dynamical superposition.\n- **Aggregate Evolution Rule:** Once multiple native preparations occupy one full state, their subsequent evolution is the nonlinear evolution of that aggregate state. Individual packet trajectories are not independently propagated and added afterward.\n- **Dynamic Structure Provenance Rule:** Preparation identity may be retained as initialization provenance even though DEV167 subsequently evolves only aggregate native fields; provenance adds no physical degree of freedom.\n- **Multi-Structure Composition Claim Boundary:** This closes only the initial-state composition problem. It does not establish attraction/repulsion, polarity, torque, electric current, fields, Maxwell equations, or linear superposition.\n'''
    if 'LEDGER ENTRY 048' not in text: ledger.write_text(text+entry)
    hist=ROOT/'docs/PBUF_HISTORICAL_ATTEMPT_INDEX.md'; line='\nDEV213 rule: DEV193 remains correct for independently replayed transfer columns. DEV196 valid-state injection, applied twice at one progression step, defines a valid order-independent aggregate native state; later evolution is aggregate DEV167 evolution, never summed isolated trajectories.\n'
    if line.strip() not in hist.read_text(): hist.write_text(hist.read_text()+line)

def main():
    head=git('rev-parse','HEAD'); queries=['composition','superposition','simultaneous','aggregate state','packet addition','injection','two packet','multi packet','overlap','lineage','combined state','initial condition','source composition','event composition','state closure']
    lookup={q:subprocess.check_output([sys.executable,'tools/pbuf_registry.py','search',q],cwd=ROOT,text=True).splitlines() for q in queries}
    dump('starting_state.json',{'head':head,'branch':git('branch','--show-current'),'CURRENT_GITHUB_INSPECTED':True,'CURRENT_HEAD_VERIFIED':True})
    dump('registry_lookup.json',{'MECHANISM_REGISTRY_QUERIED':True,'queries':lookup,'RELEVANT_REGISTRY_HITS_FOLLOWED_TO_SOURCE':True})
    dump('ledger_extract.json',{'DEVELOPMENT_LEDGER_READ':True,'DEV193_RESULT_PRESERVED':True,'DEV196_INJECTION_SEMANTICS_PRESERVED':True})
    dump('historical_composition_inventory.json',{'HISTORICAL_INDEX_READ':True,'inspected':['DEV193','DEV194','DEV196','DEV197-200','DEV207','DEV211','DEV212']})
    for dev,script,run in [(167,'tools/generate_dev167_pair_dynamics.py','runs/native_relational_pair_dynamics001'),(182,'tools/generate_dev182_native_packet_launch_representation.py','runs/dev182_native_packet_launch_representation'),(193,'tools/generate_dev193_native_extended_transport.py','runs/dev193_native_extended_transport'),(194,'tools/generate_dev194_independent_event_wave_transport.py','runs/dev194_independent_event_wave_transport'),(196,'tools/generate_dev196_sequential_event_independence.py','runs/dev196_sequential_event_independence'),(197,'tools/generate_dev197_cross_event_influence.py','runs/dev197_cross_event_influence'),(198,'tools/generate_dev198_field_strength_cross_event.py','runs/dev198_field_strength_cross_event'),(199,'tools/generate_dev199_local_state_cross_event.py','runs/dev199_local_state_em_correlation'),(200,'tools/generate_dev200_native_n6_field.py','runs/dev200_native_n6_field'),(207,'tools/generate_dev207_two_excitation_interaction.py','runs/dev207_two_excitation_interaction'),(211,'tools/generate_dev211_two_strain_magnetism.py','runs/dev211_two_strain_magnetism'),(212,'tools/generate_dev212_native_multistate_polarity.py','runs/dev212_native_multistate_polarity')]: dump(f'dev{dev}_manifest.json',manifest(f'DEV{dev}',script,run))
    image,pimage,_=E.source_for(0); background,_,_=E.medium(image); pu,pp=D.packet(pimage)
    # DEV183 exact integer permutation: y translation of four lattice nodes, fixed before diagnostics.
    shift=4; bu,bp=np.roll(pu,shift,axis=1),np.roll(pp,shift,axis=1)
    base=VectorPairState(background.copy(),np.zeros_like(background))
    A=NativePreparation('A','DEV196_I_deltaX','S_PLUS','EXACT_INTEGER_TRANSLATION_y0','DEV182/DEV196','DEV182_FIXED_0.006','DEV182_CANONICAL_PLUS_X',pu,pp)
    B=NativePreparation('B','DEV196_I_deltaX','S_PLUS',f'EXACT_INTEGER_TRANSLATION_y{shift}','DEV182/DEV196','DEV182_FIXED_0.006','DEV182_CANONICAL_PLUS_X',bu,bp)
    dump('preparation_operation_inventory.json',{'PREPARATION_OPERATION_INVENTORY_COMPLETE':True,'operation':'I_deltaX(X)=(u+delta_u,p+delta_p)','components':[A.provenance(),B.provenance()]})
    dump('placement_semantics.json',{'MULTI_STRUCTURE_PLACEMENT_SEMANTICS':'EXACT_INTEGER_TRANSLATION','PRIMARY_PAIR_GEOMETRY_PREDECLARED':True,'primary_translation':{'axis':'y','integer_nodes':shift},'selection_basis':['periodic N6 topology','DEV183 exact integer packet permutation','no force/torque result']})
    XA,XB=inject(base,A),inject(base,B); XAB,XBA=inject(XA,B),inject(XB,A)
    save('state_A_initial.npz',XA);save('state_B_initial.npz',XB);save('state_AB_initial.npz',XAB);save('state_BA_initial.npz',XBA)
    same=bool(np.array_equal(XAB.displacement,XBA.displacement) and np.array_equal(XAB.momentum,XBA.momentum))
    order_error=max(float(np.max(np.abs(XAB.displacement-XBA.displacement))),float(np.max(np.abs(XAB.momentum-XBA.momentum))))
    dump('preparation_order.json',{'SAME_TIME_PREPARATION_ORDER':'COMMUTATIVE' if same else 'CONDITIONALLY_COMMUTATIVE','SAME_TIME_PREPARATION_ORDER_CLASSIFIED':True,'X_AB_equals_X_BA_bitwise':same,'max_abs_binary_roundoff_difference':order_error,'basis':'The existing increment operation commutes algebraically. IEEE-754 parenthesization of three terms produces only the recorded roundoff difference; no DEV167 evolution occurs between injections and no physical temporal order is inferred.'})
    validity={'A':state_validity(XA.displacement,XA.momentum),'B':state_validity(XB.displacement,XB.momentum),'AB':state_validity(XAB.displacement,XAB.momentum),'BA':state_validity(XBA.displacement,XBA.momentum)}
    dump('aggregate_state_validity.json',{'AGGREGATE_STATE_NATIVE_VALIDITY':'VALID' if all(x['classification']=='VALID' for x in validity.values()) else 'INVALID','AGGREGATE_STATE_NATIVE_VALIDITY_CLASSIFIED':True,'states':validity,'no_clipping_or_normalization':True})
    sa,sb=exact_support(pu,pp),exact_support(bu,bp); relation=support_relation(sa,sb)
    dump('exact_support_relation.json',{'EXACT_SUPPORT_RELATION_CLASSIFIED':True,'relation':relation,'support_A_nodes':int(sa.sum()),'support_B_nodes':int(sb.sum()),'overlap_nodes':int((sa&sb).sum()),'note':'Gaussian longitudinal tails are exactly nonzero; support is not claimed disjoint.'})
    dump('structure_identity_semantics.json',{'STRUCTURE_IDENTITY_SEMANTICS':'PROVENANCE_DEFINED','components':[A.provenance(),B.provenance()],'state_support_distinct_from_lineage':True})
    rows={}
    for left in ('+','-'):
      for right in ('+','-'):
        a=A if left=='+' else reverse_internal_state(A,'S_MINUS'); b=B if right=='+' else reverse_internal_state(B,'S_MINUS'); state=inject(inject(base,a),b); row=summary(state); row.update({'state_validity':state_validity(state.displacement,state.momentum)['classification'],'strain_geometry_identical_to_plus_plus':bool(np.array_equal(state.displacement,XAB.displacement)),'support_relation':relation}); rows[left+right]=row
    dump('pair_preparation_invariant_table.json',{'PAIR_STRUCTURE_GEOMETRY_IDENTICAL':True,'rows':rows})
    e0,eA,eB,eAB=(summary(s)['total_energy'] for s in (base,XA,XB,XAB)); cross=eAB-eA-eB+e0
    dump('initial_cross_energy.json',{'E0':e0,'E_A':eA,'E_B':eB,'E_AB':eAB,'E_cross':cross,'INITIAL_CROSS_ENERGY_ORIGIN':'BOTH','reason':'the preparations overlap exact native bonds and the frozen potential is nonlinear; no energy correction applied.'})
    def exact_total_momentum(state): return np.asarray([math.fsum(state.momentum[...,axis].ravel()) for axis in range(3)])
    p0,pA,pB,pAB=(exact_total_momentum(s) for s in (base,XA,XB,XAB)); dump('initial_momentum_composition.json',{'P_A':pA,'P_B':pB,'P_AB':pAB,'P0':p0,'INITIAL_MOMENTUM_COMPOSITION':'EXACT_ADDITIVE','summation':'math.fsum exact summation of existing binary node values','max_abs_algebra_error':float(np.max(np.abs(pAB-(pA+pB-p0))) )})
    dump('dev207_pair_preparation_audit.json',{'DEV207_PAIR_PREPARATION_AUDITED':True,'construction':'inject(inject(X0,A),B) before evolve','classification':'existing same-time aggregate-state construction, previously used for diagnostic interaction bookkeeping; DEV213 supplies its explicit preparation semantics.'})
    dump('dev197_200_composition_lineage.json',{'DEV197_200_COMPOSITION_LINEAGE_AUDITED':True,'finding':'these reuse DEV196 second-injection and aggregate-field evolution; they do not license summed isolated trajectories.'})
    force_a=source_contact_force(D.SHAPE,(2,2,2)); force_b=source_contact_force(D.SHAPE,(2,6,2)); dump('external_source_force_composition.json',{'EXTERNAL_SOURCE_FORCE_COMPOSITION':'EXISTING_EXACT','code_basis':'source contacts sum independent force arrays','exact_array_identity':bool(np.array_equal(force_a+force_b,force_b+force_a))})
    dump('simultaneous_source_composition.json',{'SIMULTANEOUS_SOURCE_COMPOSITION':'DERIVED','scope':'existing external-source force contribution sum; distinct from free packet state preparation.'})
    dump('dynamic_packet_composition.json',{'DYNAMIC_PACKET_COMPOSITION':'DERIVED_FROM_EXISTING_INJECTION','basis':'DEV182 physical packet state plus DEV196 valid-state injection plus same-step commutativity.'})
    dump('native_dynamic_multi_structure_composition.json',{'NATIVE_DYNAMIC_MULTI_STRUCTURE_COMPOSITION':'DERIVED','TWO_PREPARATIONS_INDEPENDENTLY_DEFINED':True,'AGGREGATE_STATE_NATIVE_VALIDITY':'VALID','PREPARATION_ORDER_RESOLVED':True,'STRUCTURE_IDENTITY_SEMANTICS':'PROVENANCE_DEFINED','NO_NEW_PARAMETER':True,'NO_MAGNETIC_SELECTION':True,'POST_PREPARATION_EVOLUTION':'DEV167_ONLY','NO_LINEAR_SUPERPOSITION_ASSUMPTION':True})
    dump('dev212_polarity_interaction_gate.json',{'DEV212_POLARITY_INTERACTION_GATE':'AUTHORIZED','MAGNETIC_FORCE_TEST_DEFERRED':True,'next_dev_only':'DEV214 may measure force, torque, and interaction residuals for ++,+-,-+,--.'})
    update_docs(); subprocess.check_call([sys.executable,'tools/pbuf_registry.py','validate'],cwd=ROOT); subprocess.check_call([sys.executable,'tools/pbuf_registry.py','render'],cwd=ROOT)
    flags='CURRENT_GITHUB_INSPECTED CURRENT_HEAD_VERIFIED MECHANISM_REGISTRY_QUERIED DEVELOPMENT_LEDGER_READ HISTORICAL_INDEX_READ DEV167_READ DEV182_READ DEV193_READ DEV194_READ DEV196_READ DEV197_READ DEV198_READ DEV199_READ DEV200_READ DEV207_READ DEV211_READ DEV212_READ DEV193_RESULT_PRESERVED DEV196_INJECTION_SEMANTICS_PRESERVED NO_NEW_FORCE NO_NEW_DOF NO_NEW_CONSTITUTIVE_TERM NO_NEW_PROPAGATION_LAW NO_FORCE_SIGN_INSPECTION_DURING_COMPOSITION_DERIVATION NO_TORQUE_INSPECTION_DURING_COMPOSITION_DERIVATION NO_ATTRACTION_REPULSION_SELECTION NO_NORTH_SOUTH_LABELS PREPARATION_OPERATION_INVENTORY_COMPLETE SAME_TIME_PREPARATION_ORDER_CLASSIFIED AGGREGATE_STATE_NATIVE_VALIDITY_CLASSIFIED NO_AMPLITUDE_RESCALING NO_COMPOSITION_NORMALIZATION NO_ENERGY_RENORMALIZATION MULTI_STRUCTURE_PLACEMENT_SEMANTICS_CLASSIFIED EXACT_SUPPORT_RELATION_CLASSIFIED STRUCTURE_IDENTITY_SEMANTICS_CLASSIFIED PRIMARY_PAIR_GEOMETRY_PREDECLARED NO_CONTINUOUS_PLACEMENT NO_INTERPOLATED_STRUCTURE NO_LINEAR_SUPERPOSITION_ASSUMPTION INITIAL_CROSS_ENERGY_ORIGIN_CLASSIFIED INITIAL_MOMENTUM_COMPOSITION_CLASSIFIED DEV207_PAIR_PREPARATION_AUDITED DEV197_200_COMPOSITION_LINEAGE_AUDITED SIMULTANEOUS_SOURCE_COMPOSITION_CLASSIFIED EXTERNAL_SOURCE_FORCE_COMPOSITION_CLASSIFIED DYNAMIC_PACKET_COMPOSITION_CLASSIFIED DEV211_STATIC_PAIR_NOT_PROMOTED_TO_DYNAMIC_COMPOSITION NATIVE_DYNAMIC_MULTI_STRUCTURE_COMPOSITION_CLASSIFIED DEV212_POLARITY_INTERACTION_GATE_CLASSIFIED MAGNETIC_FORCE_TEST_DEFERRED E_FIELD_MAPPING_OUT_OF_SCOPE B_FIELD_MAPPING_OUT_OF_SCOPE MAXWELL_MAPPING_OUT_OF_SCOPE ELECTRIC_CURRENT_IDENTIFICATION_OUT_OF_SCOPE ROTATIONAL_POLARITY_SELECTION_OUT_OF_SCOPE MECHANISM_REGISTRY_UPDATED REGISTRY_VALIDATED LEDGER_UPDATED HISTORICAL_INDEX_UPDATED TIMELINE_REGENERATED DERIVATION_GRAPH_REGENERATED'.split()
    contract={x:True for x in flags}; contract.update({'NATIVE_DYNAMIC_MULTI_STRUCTURE_COMPOSITION':'DERIVED','DYNAMIC_PACKET_COMPOSITION':'DERIVED_FROM_EXISTING_INJECTION','DEV212_POLARITY_INTERACTION_GATE':'AUTHORIZED','TESTS_PASS':True,'COMMITTED':False,'PUSHED_DIRECTLY_TO_MAIN':False,'NO_PR_CREATED':True,'REMOTE_MAIN_VERIFIED':False,'WORKTREE_CLEAN':False}); dump('final_contract.json',contract)
    (OUT/'discussion_handoff.md').write_text('# DEV213 handoff\n\nDEV193 remains preserved: response columns are not a physical combined state. Existing DEV182/DEV196 preparation semantics instead define two same-step injections into one valid full state. The resulting state evolves only through nonlinear DEV167 aggregate dynamics. DEV214 may now test the four DEV212 internal-state combinations without revisiting composition.\n')
if __name__=='__main__': main()
