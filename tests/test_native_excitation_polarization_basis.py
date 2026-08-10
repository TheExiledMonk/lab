from pbuf.quantum.native_excitation_modes import rotating_mode
from pbuf.quantum.native_excitation_interference import basis_invariance_audit,handedness
def test_basis_norm_and_handedness():
    assert basis_invariance_audit(rotating_mode(128,16))['invariant']; assert handedness(rotating_mode(128,16,handedness=1))!=handedness(rotating_mode(128,16,handedness=-1))

