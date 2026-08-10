import numpy as np
from pbuf.matter.native_excitation_invariants import conservation_audit, norm_audit, scalar_identity_invariant

def test_norms_are_audited_without_asserting_joint_invariant():
    assert len(norm_audit())==9 and len(conservation_audit())==7
    result=scalar_identity_invariant(np.ones(8))
    assert result['conserved'] and not result['physical_energy_norm']

