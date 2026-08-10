from pbuf.excitation.native_excitation_invariants import invariant_audit
from pbuf.excitation.native_excitation_state import NativeExcitationState, localized_packet
from pbuf.excitation.native_excitation_transfer import progress_source_free

def test_quadratic_norm_is_not_renormalized_and_is_conserved():
    s=NativeExcitationState(localized_packet(amplitude=3)); progress_source_free(s,40)
    a=invariant_audit(s.history); assert a['conserved'] and a['relative_drift']<1e-14

