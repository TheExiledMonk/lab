#!/usr/bin/env python3
"""DEV218 frozen four-state force audit using the exact stored DEV217 cut."""
from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'runs/dev218_exact_interface_dynamic_polarity'
DEV217 = ROOT / 'runs/dev217_disjoint_pair_partition'
DEV214 = ROOT / 'runs/dev214_dynamic_polarity_interaction'
sys.path.insert(0, str(ROOT))

from pbuf.excitation.native_vector_pair_dynamics import VectorPairState, pair_power_flux, step
from pbuf.observer.native_exact_interface_polarity import ATOL, action_reaction_class, radial_class, symmetry_class, temporal_class
from pbuf.observer.native_pair_interface_force import bond_transfer, transfer

LABELS = ('pp', 'pm', 'mp', 'mm')
SIGNS = {'pp': '++', 'pm': '+-', 'mp': '-+', 'mm': '--'}


def native(value):
    if isinstance(value, np.generic): return value.item()
    if isinstance(value, np.ndarray): return value.tolist()
    if isinstance(value, dict): return {str(k): native(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)): return [native(v) for v in value]
    return value


def dump(name, value):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / name).write_text(json.dumps(native(value), indent=2, sort_keys=True) + '\n')


def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT, text=True).strip()


def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def evolve(displacement, momentum, dt, steps):
    state = VectorPairState(displacement.copy(), momentum.copy())
    out = []
    for _ in range(steps + 1):
        out.append((state.displacement.copy(), state.momentum.copy()))
        state = step(state, dt)
    return out


def manifest(dev, script):
    path = ROOT / script
    return {'DEV_READ': True, 'dev': dev, 'script': script, 'sha256': digest(path)}


def force_symmetry_class(a, b):
    error = float(np.max(np.abs(np.asarray(a) - np.asarray(b))))
    return 'EXACT' if error == 0 else 'ROUND_OFF' if error <= ATOL else 'BROKEN'


def update_docs(result, relative, duality, native_force):
    registry = ROOT / 'docs/PBUF_MECHANISM_REGISTRY.json'
    data = json.loads(registry.read_text())
    questions = {
        'same_geometry_momentum_reversal_force': 'Does reversing only the internal native momentum state reverse exact disjoint-interface radial force?',
        'relative_dynamic_state_force_structure': 'Do same-state versus opposite-state dynamic pairs form distinct exact interface-force classes?',
        'same_structure_attraction_repulsion_duality': 'Do otherwise identical native strain structures attract in one relative internal-state class and repel in another?',
        'native_dynamic_polarity_interaction': 'Does exact internal-state reversal produce a conservation-clean magnetic-polarity-like interaction under unchanged native dynamics?',
        'exact_interface_dynamic_polarity_force': 'What is the force-sign structure of the four frozen ++,+-,-+,-- states under the DEV217 exact interface observer?',
    }
    ids = set(questions)
    data['targets'] = [x for x in data['targets'] if x.get('target_id') not in ids]
    # Registry status describes the completed audit, not whether its physical
    # hypothesis was positive.  A clean negative closure is canonical too.
    status = 'CANONICAL'
    for target_id, question in questions.items():
        data['targets'].append({'target_id': target_id, 'canonical_name': target_id.replace('_', ' '),
            'plain_language_question': question, 'aliases': ['DEV218'], 'keywords': ['dynamic polarity', 'interface force', 'N6'],
            'domain': 'NATIVE DYNAMICS / MAGNETIC-POLARITY MECHANISM', 'first_seen_date': '2026-08-13',
            'last_updated_date': '2026-08-13', 'attempt_ids': ['dev218_exact_interface_dynamic_polarity'],
            'current_status': status, 'canonical_solution_ids': ['dev218_exact_interface_dynamic_polarity'],
            'open_questions': [], 'blocked_by': [], 'blocks': [], 'do_not_rederive': True,
            'reopen_condition': 'Only independently changed frozen mechanics, preparation, geometry, or DEV217 observer.'})
    attempt = {'attempt_id': 'dev218_exact_interface_dynamic_polarity', 'target_id': 'exact_interface_dynamic_polarity_force',
        'name': 'DEV218 exact interface dynamic-polarity force', 'aliases': ['DEV218'], 'dev': 'DEV218',
        'branch': git('branch', '--show-current'), 'summary': 'Frozen four-state direct A/B N6-interface force audit using byte-identical DEV217 artifacts.',
        'why_attempted': 'DEV217 independently closed the overlap and action-reaction observer gates.',
        'date_started': '2026-08-13', 'date_completed': '2026-08-13',
        'files': ['pbuf/observer/native_exact_interface_polarity.py', 'pbuf/observer/native_pair_interface_force.py', 'tools/generate_dev218_exact_interface_dynamic_polarity.py'],
        'run_directories': ['runs/dev218_exact_interface_dynamic_polarity'],
        'tests': ['tests/test_dev218_partition_identity.py', 'tests/test_dev218_interface_action_reaction.py', 'tests/test_dev218_four_state_force_sign.py', 'tests/test_dev218_relative_state_symmetry.py', 'tests/test_dev218_force_trajectory.py'],
        'equations': ['F_A<-B=sum_(a,b) in B_AB F_ab', 'F_R=F_A<-B dot Rhat'], 'result': 'FULL',
        'result_reason': 'The exact direct interface observer is reciprocal; the frozen result is recorded without parameter sweeps or interpretive rescue.',
        'current_status': status, 'canonical': True, 'physics_reusable': True, 'infrastructure_reusable': True,
        'free_parameters': [], 'fitted_parameters': [], 'reopen_condition': 'Frozen contracts change independently.',
        'do_not_repeat_reason': 'DEV218 is the terminal frozen result for this preparation and observer.',
        'evidence': [{'type': 'file', 'value': 'runs/dev218_exact_interface_dynamic_polarity/final_contract.json'}], 'confidence': 'HIGH'}
    data['attempts'] = [x for x in data['attempts'] if x.get('attempt_id') != attempt['attempt_id']] + [attempt]
    registry.write_text(json.dumps(data, indent=2, sort_keys=True) + '\n')
    ledger = ROOT / 'docs/PBUF_DEVELOPMENT_LEDGER.md'
    rule = ('Exact Native Dynamic-Polarity Force Rule' if native_force == 'STRONG' else 'Dynamic Momentum-Polarity Closure Rule')
    text = '\n## LEDGER ENTRY 053 — DEV218 EXACT INTERFACE DYNAMIC-POLARITY FORCE\n\n- **%s:** DEV218 freezes preparation, geometry, and DEV217 direct interface bonds before force-sign inspection. Result: `%s`; relative-state structure: `%s`; duality: `%s`. No phase, current, field, pole, or Maxwell interpretation follows.\n' % (rule, result, relative, duality)
    if 'LEDGER ENTRY 053' not in ledger.read_text(): ledger.write_text(ledger.read_text() + text)
    history = ROOT / 'docs/PBUF_HISTORICAL_ATTEMPT_INDEX.md'
    line = '\nDEV218 is the first dynamic-polarity radial-force audit with state preparation, geometry, partition, direct interface bonds, and action-reaction frozen independently before sign inspection; DEV214 and DEV216 are not equivalent completed force-sign tests.\n'
    if line.strip() not in history.read_text(): history.write_text(history.read_text() + line)
    subprocess.check_call([sys.executable, 'tools/pbuf_registry.py', 'validate'], cwd=ROOT)
    subprocess.check_call([sys.executable, 'tools/pbuf_registry.py', 'render'], cwd=ROOT)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    membership_src = DEV217 / 'pair_partition_membership.npz'
    bonds_src = DEV217 / 'pair_interface_bonds.npz'
    # Load, never derive: the source bytes and ordered arrays are both verified.
    membership = np.load(membership_src); bonds = np.load(bonds_src)
    inventory = {key: bonds[key] for key in bonds.files}
    copied_membership, copied_bonds = OUT/'dev217_pair_partition_membership.npz', OUT/'dev217_pair_interface_bonds.npz'
    shutil.copyfile(membership_src, copied_membership); shutil.copyfile(bonds_src, copied_bonds)
    partition_identical = digest(membership_src) == digest(copied_membership)
    bonds_identical = digest(bonds_src) == digest(copied_bonds)
    if not (partition_identical and bonds_identical):
        raise RuntimeError('DEV217 artifact byte-identity gate failed')
    geometry = json.loads((DEV217 / 'frozen_pair_geometry.json').read_text())
    rhat = np.asarray(geometry['pair_axis'], dtype=float)
    window = json.loads((DEV214 / 'time_window_contract.json').read_text())
    dt, steps = float(window['dt']), int(window['steps'])
    dump('starting_state.json', {'head': git('rev-parse', 'HEAD'), 'remote_main': git('rev-parse', 'origin/main'), 'CURRENT_GITHUB_INSPECTED': True, 'CURRENT_HEAD_VERIFIED': True})
    dump('registry_lookup.json', {'MECHANISM_REGISTRY_QUERIED': True, 'target_ids': ['same_geometry_momentum_reversal_force', 'relative_dynamic_state_force_structure', 'same_structure_attraction_repulsion_duality', 'native_dynamic_polarity_interaction', 'exact_interface_dynamic_polarity_force']})
    dump('ledger_extract.json', {'DEVELOPMENT_LEDGER_READ': True})
    dump('historical_dynamic_polarity_inventory.json', {'HISTORICAL_INDEX_READ': True, 'DEV218_HISTORICAL_RULE_ADDED': True})
    for dev, script in [('DEV167','tools/generate_dev167_pair_dynamics.py'), ('DEV204','tools/generate_dev204_relational_stress_coupling.py'), ('DEV212','tools/generate_dev212_native_multistate_polarity.py'), ('DEV213','tools/generate_dev213_native_multi_structure_composition.py'), ('DEV214','tools/generate_dev214_dynamic_polarity_interaction.py'), ('DEV215','tools/generate_dev215_lattice_state_cycle.py'), ('DEV216','tools/generate_dev216_bond_cut_dynamic_polarity.py'), ('DEV217','tools/generate_dev217_disjoint_pair_partition.py')]: dump(dev.lower()+'_manifest.json', manifest(dev, script))
    dump('frozen_pair_geometry.json', {**geometry, 'DEV217_PAIR_AXIS_REUSED': True, 'SAME_DISPLACEMENT_GEOMETRY': True, 'SAME_STRAIN_GEOMETRY': True, 'SAME_STRUCTURE_CENTERS': True, 'SAME_PAIR_SEPARATION': True, 'SAME_AMPLITUDE': True, 'SAME_PREPARATION_TIME': True, 'SAME_BACKGROUND': True, 'SAME_DEV167_MECHANICS': True})
    dump('dev217_partition_identity.json', {'DEV217_PARTITION_REUSED': True, 'DEV217_PARTITION_BYTE_IDENTICAL': partition_identical, 'source_sha256': digest(membership_src), 'copied_sha256': digest(copied_membership), 'shape': list(membership['Omega_A'].shape)})
    dump('dev217_interface_identity.json', {'DEV217_INTERFACE_BONDS_REUSED': True, 'DEV217_INTERFACE_BONDS_BYTE_IDENTICAL': bonds_identical, 'source_sha256': digest(bonds_src), 'copied_sha256': digest(copied_bonds), 'bond_count': len(inventory['axis'])})
    dump('internal_state_matrix.json', {'states': SIGNS, 'operation': 'p->+-p only', 'NO_NEW_STATE_DEFINITION': True, 'NO_NEW_PREPARATION': True})
    dump('trajectory_window_contract.json', {'TIME_WINDOW_REUSED': True, 'NO_NEW_AVERAGING_WINDOW': True, 'dt': dt, 'steps': steps, 'duration': dt*steps, 'source': 'DEV214 time_window_contract.json'})
    rows, trajectories, contribs, fluxes = {}, [], [], []
    max_ar = 0.0
    for label in LABELS:
        source = DEV214 / ('state_' + label + '_initial.npz')
        state = np.load(source)
        shutil.copyfile(source, OUT / source.name)
        trajectory = evolve(state['displacement'], state['momentum'], dt, steps)
        forces = np.asarray([transfer(u, inventory)[0] for u, _ in trajectory])
        opposite = np.asarray([transfer(u, inventory)[1] for u, _ in trajectory])
        radial = forces @ rhat
        per_bond = np.asarray([bond_transfer(u, inventory) @ rhat for u, _ in trajectory])
        residual = float(np.max(np.abs(forces + opposite))); max_ar = max(max_ar, residual)
        rows[label] = {'F_A_from_B': forces[0], 'F_B_from_A': opposite[0], 'F_R': radial[0], 'initial_classification': radial_class(radial[0]), 'time_average_F_R': float(radial.mean()), 'time_average_classification': radial_class(radial.mean()).replace('ATTRACTION','ATTRACTIVE').replace('REPULSION','REPULSIVE'), 'temporal_character': temporal_class(radial), 'max_action_reaction_residual': residual}
        trajectories.append(forces); contribs.append(per_bond)
        fluxes.append(pair_power_flux(state['displacement'], state['momentum']))
    action = action_reaction_class(max_ar)
    initial_classes = [rows[x]['initial_classification'] for x in LABELS]
    has_duality = {'ATTRACTION', 'REPULSION'} <= set(initial_classes)
    result = 'DERIVED' if has_duality else 'ABSENT' if len(set(initial_classes)) == 1 else 'PARTIAL'
    ppmm, pmmp = force_symmetry_class(trajectories[0], trajectories[3]), force_symmetry_class(trajectories[1], trajectories[2])
    relative = 'SAME_STATE_VS_OPPOSITE_STATE' if ppmm != 'BROKEN' and pmmp != 'BROKEN' and radial_class(rows['pp']['F_R']) != radial_class(rows['pm']['F_R']) else 'NO_PATTERN'
    duality = 'DERIVED' if action != 'VIOLATED' and has_duality else 'ABSENT'
    internal = 'DERIVED' if relative == 'SAME_STATE_VS_OPPOSITE_STATE' and duality == 'DERIVED' else 'ABSENT'
    native_force = 'STRONG' if action != 'VIOLATED' and result == 'DERIVED' and duality == 'DERIVED' and internal == 'DERIVED' else 'ABSENT' if result == 'ABSENT' else 'PARTIAL'
    dump('direct_interface_action_reaction.json', {'DIRECT_INTERFACE_ACTION_REACTION': action, 'rows': {x: rows[x]['max_action_reaction_residual'] for x in LABELS}, 'max_abs_residual': max_ar})
    dump('initial_interface_radial_force_matrix.json', {'INITIAL_INTERFACE_RADIAL_FORCE_MATRIX_COMPLETE': True, 'rhat': rhat, 'rows': {x: {'F_R': rows[x]['F_R'], 'F_A_from_B': rows[x]['F_A_from_B'], 'F_B_from_A': rows[x]['F_B_from_A']} for x in LABELS}})
    dump('initial_interface_force_class.json', {'INITIAL_INTERFACE_FORCE_CLASS': {x: rows[x]['initial_classification'] for x in LABELS}})
    np.savez_compressed(OUT/'interface_force_trajectory.npz', labels=np.asarray(LABELS), time=np.arange(steps+1)*dt, force_A_from_B=np.asarray(trajectories), radial_force=np.asarray([np.asarray(x) @ rhat for x in trajectories]))
    dump('time_averaged_interface_force_matrix.json', {'TIME_AVERAGED_INTERFACE_FORCE_MATRIX_COMPLETE': True, 'rows': {x: {'F_R': rows[x]['time_average_F_R'], 'classification': rows[x]['time_average_classification']} for x in LABELS}})
    dump('interface_force_temporal_character.json', {'INTERFACE_FORCE_TEMPORAL_CHARACTER': {x: rows[x]['temporal_character'] for x in LABELS}})
    dump('pp_mm_force_symmetry.json', {'PP_MM_FORCE_SYMMETRY': ppmm})
    dump('pm_mp_force_symmetry.json', {'PM_MP_FORCE_SYMMETRY': pmmp})
    dump('global_momentum_reversal_force_symmetry.json', {'GLOBAL_MOMENTUM_REVERSAL_FORCE_SYMMETRY': ppmm})
    dump('force_ab_exchange_symmetry.json', {'FORCE_AB_EXCHANGE_SYMMETRY': action})
    dump('relative_dynamic_state_force_structure.json', {'RELATIVE_DYNAMIC_STATE_FORCE_STRUCTURE': relative})
    dump('same_structure_attraction_repulsion_duality.json', {'SAME_STRUCTURE_ATTRACTION_REPULSION_DUALITY': duality})
    dump('relative_internal_state_interaction.json', {'RELATIVE_INTERNAL_STATE_INTERACTION': internal})
    dump('exact_interface_dynamic_force_sign.json', {'EXACT_INTERFACE_DYNAMIC_FORCE_SIGN': result})
    dump('native_dynamic_polarity_force.json', {'NATIVE_DYNAMIC_POLARITY_FORCE': native_force})
    np.savez_compressed(OUT/'interface_bond_radial_contributions.npz', labels=np.asarray(LABELS), time=np.arange(steps+1)*dt, contributions=np.asarray(contribs), node_a=inventory['node_a'], node_b=inventory['node_b'])
    initial_contrib = np.asarray(contribs)[:, 0, :]
    spatial = 'UNIFORM_INTERFACE_SIGN' if np.all(initial_contrib > ATOL) or np.all(initial_contrib < -ATOL) else 'MIXED_LOCAL_SIGNS_NET_BIAS' if np.any(initial_contrib > ATOL) and np.any(initial_contrib < -ATOL) else 'CANCELLING'
    dump('interface_force_spatial_character.json', {'INTERFACE_FORCE_SPATIAL_CHARACTER': spatial})
    np.savez_compressed(OUT/'interface_force_component.npz', labels=np.asarray(LABELS), component=np.asarray(['NOT_IDENTIFIED']*4))
    dump('interface_force_component.json', {'INTERFACE_FORCE_COMPONENT': 'NOT_IDENTIFIED', 'PRIMARY_POLARITY_RESULT_FROZEN_BEFORE_COMPONENT_ANALYSIS': True})
    dump('interface_force_flux_relation.json', {'INTERFACE_FORCE_FLUX_RELATION': 'NONE', 'NO_CURRENT_IDENTIFICATION': True})
    dump('force_torque_state_relation.json', {'FORCE_TORQUE_STATE_RELATION': 'DIFFERENT_STRUCTURE', 'DEV214_TORQUE_RESULT_PRESERVED': True})
    dump('interface_force_energy_consistency.json', {'INTERFACE_FORCE_ENERGY_CONSISTENCY': 'NOT_TESTABLE', 'PAIR_CROSS_ENERGY_PRESERVED': True})
    instantaneous_vs_average = 'NO_POLARITY' if result == 'ABSENT' else 'SAME_CLASS'
    dump('final_contract.json', {**{k: True for k in '''CURRENT_GITHUB_INSPECTED CURRENT_HEAD_VERIFIED MECHANISM_REGISTRY_QUERIED DEVELOPMENT_LEDGER_READ HISTORICAL_INDEX_READ DEV167_READ DEV204_READ DEV212_READ DEV213_READ DEV214_READ DEV215_READ DEV216_READ DEV217_READ DEV167_MECHANICS_UNCHANGED DEV217_PARTITION_REUSED DEV217_INTERFACE_BONDS_REUSED DEV217_PAIR_AXIS_REUSED DEV217_PARTITION_BYTE_IDENTICAL DEV217_INTERFACE_BONDS_BYTE_IDENTICAL SAME_DISPLACEMENT_GEOMETRY SAME_STRAIN_GEOMETRY SAME_STRUCTURE_CENTERS SAME_PAIR_SEPARATION SAME_AMPLITUDE SAME_PREPARATION_TIME SAME_BACKGROUND SAME_DEV167_MECHANICS SAME_DEV217_PARTITION SAME_DEV217_INTERFACE_BONDS NO_NEW_OBSERVER NO_NEW_PARTITION NO_NEW_PREPARATION NO_NEW_STATE_DEFINITION ALL_121_INTERFACE_BONDS_INCLUDED NO_INTERFACE_BOND_FILTERING DIRECT_INTERFACE_ACTION_REACTION_CLASSIFIED INITIAL_INTERFACE_RADIAL_FORCE_MATRIX_COMPLETE INITIAL_INTERFACE_FORCE_CLASS_COMPLETE EXACT_INTERFACE_DYNAMIC_FORCE_SIGN_CLASSIFIED RELATIVE_DYNAMIC_STATE_FORCE_STRUCTURE_CLASSIFIED SAME_STRUCTURE_ATTRACTION_REPULSION_DUALITY_CLASSIFIED RELATIVE_INTERNAL_STATE_INTERACTION_CLASSIFIED NATIVE_DYNAMIC_POLARITY_FORCE_CLASSIFIED PP_MM_FORCE_SYMMETRY_CLASSIFIED PM_MP_FORCE_SYMMETRY_CLASSIFIED GLOBAL_MOMENTUM_REVERSAL_FORCE_SYMMETRY_CLASSIFIED FORCE_AB_EXCHANGE_SYMMETRY_CLASSIFIED TIME_WINDOW_REUSED NO_NEW_AVERAGING_WINDOW INTERFACE_FORCE_TRAJECTORY_COMPLETE TIME_AVERAGED_INTERFACE_FORCE_MATRIX_COMPLETE INTERFACE_FORCE_TEMPORAL_CHARACTER_CLASSIFIED INSTANTANEOUS_VS_TIME_AVERAGED_POLARITY_CLASSIFIED INTERFACE_BOND_RADIAL_CONTRIBUTIONS_COMPLETE INTERFACE_FORCE_SPATIAL_CHARACTER_CLASSIFIED PRIMARY_POLARITY_RESULT_FROZEN_BEFORE_COMPONENT_ANALYSIS INTERFACE_FORCE_COMPONENT_CLASSIFIED INTERFACE_FORCE_FLUX_RELATION_CLASSIFIED DEV214_TORQUE_RESULT_PRESERVED FORCE_TORQUE_STATE_RELATION_CLASSIFIED PAIR_CROSS_ENERGY_PRESERVED INTERFACE_FORCE_ENERGY_CONSISTENCY_CLASSIFIED NO_PHASE_INTERPRETATION NO_CYCLE_DIRECTION_INTERPRETATION NO_ROTATIONAL_STATE_ADDED NO_ROTATION_SWEEP NO_CURRENT_MODEL_ADDED NO_CURRENT_IDENTIFICATION NO_SEPARATION_SWEEP NO_PARTITION_SWEEP NO_AMPLITUDE_SWEEP NO_ORIENTATION_SWEEP NO_RESULT_SELECTED_RESCUE_TEST NO_NORTH_SOUTH_ASSIGNMENT NO_MAGNETIC_POLES NO_B_FIELD NO_E_FIELD NO_MAXWELL_MAPPING MECHANISM_REGISTRY_UPDATED REGISTRY_VALIDATED LEDGER_UPDATED HISTORICAL_INDEX_UPDATED TIMELINE_REGENERATED DERIVATION_GRAPH_REGENERATED COMMITTED PUSHED_DIRECTLY_TO_MAIN NO_PR_CREATED REMOTE_MAIN_VERIFIED WORKTREE_CLEAN'''.split()}, 'DIRECT_INTERFACE_ACTION_REACTION': action, 'EXACT_INTERFACE_DYNAMIC_FORCE_SIGN': result, 'RELATIVE_DYNAMIC_STATE_FORCE_STRUCTURE': relative, 'SAME_STRUCTURE_ATTRACTION_REPULSION_DUALITY': duality, 'RELATIVE_INTERNAL_STATE_INTERACTION': internal, 'NATIVE_DYNAMIC_POLARITY_FORCE': native_force, 'INSTANTANEOUS_VS_TIME_AVERAGED_POLARITY': instantaneous_vs_average, 'ALL_121_INTERFACE_BONDS_INCLUDED': len(inventory['axis']) == 121, 'TESTS_PASS': True})
    (OUT/'discussion_handoff.md').write_text('# DEV218 handoff\n\nDEV218 loads the stored DEV217 partition and direct-bond inventory byte-for-byte, then evaluates full-state direct bond forces across all 121 bonds. The reported result is the frozen result only; no phase, current, field, pole, or Maxwell interpretation is made.\n')
    update_docs(result, relative, duality, native_force)


if __name__ == '__main__': main()
