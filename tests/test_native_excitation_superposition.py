from pbuf.excitation.native_excitation_invariants import superposition_audit
from pbuf.excitation.native_excitation_state import localized_packet

def test_exact_weak_superposition():
    assert superposition_audit(localized_packet(polarization=(1,0)),localized_packet(center=50,polarization=(0,1)))['passes']

