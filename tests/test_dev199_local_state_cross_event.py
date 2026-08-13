import numpy as np
import json
from pathlib import Path
from pbuf.observer.local_state_cross_event import four_state_cross_term, sigma_prime, sigma_second
from pbuf.excitation.native_vector_pair_dynamics import pair_forces


def test_exact_bond_inclusion_exclusion_reconstructs_pair_force():
    rng=np.random.default_rng(199); x0=rng.normal(scale=1e-3,size=(3,3,3,3)); a=rng.normal(scale=1e-4,size=x0.shape); b=rng.normal(scale=1e-4,size=x0.shape)
    q=four_state_cross_term(x0,x0+a,x0+b,x0+a+b)
    direct=pair_forces(x0+a+b)-pair_forces(x0+a)-pair_forces(x0+b)+pair_forces(x0)
    assert np.array_equal(q['force_cross'],direct)
    assert np.array_equal(q['force_cross'],q['constitutive_cross']+q['geometric_cross'])


def test_constitutive_derivatives_at_equilibrium():
    assert sigma_prime(0.) == 1.
    assert sigma_second(0.) == 0.


def test_phase_a_freeze_is_native_only_and_force_reconstruction_is_exact():
    root=Path(__file__).resolve().parents[1]; out=root/'runs/dev199_local_state_em_correlation/phase_a'
    freeze=json.loads((out/'phase_A_freeze.json').read_text()); recon=json.loads((out/'measured_vs_reconstructed_cross_force.json').read_text())
    assert freeze['PHASE_A_EM_BLIND'] and recon['DeltaF_measured_equals_DeltaF_constitutive']
    assert 'em_structural' not in (root/'tools/generate_dev199_local_state_cross_event.py').read_text().lower()
