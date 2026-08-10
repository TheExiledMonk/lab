from pbuf.excitation.native_emergent_em_observer import maxwell_structure_comparison, observe_effective_pair
from pbuf.excitation.native_excitation_state import localized_packet

def test_observer_is_one_way_and_em_gate_remains_closed():
    o=observe_effective_pair(localized_packet()); assert o['observer_only'] and not o['feeds_back'] and not o['mapping_established']
    assert not maxwell_structure_comparison()['compatible']

