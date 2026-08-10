from pbuf.excitation.native_excitation_invariants import reversibility_audit
from pbuf.excitation.native_excitation_state import localized_packet

def test_exact_reverse_progression(): assert reversibility_audit(localized_packet())['reversible']

