from pbuf.excitation.native_excitation_state import NativeExcitationState, localized_packet
from pbuf.excitation.native_excitation_transfer import dependency_contract, progress_source_free
from pbuf.excitation.native_excitation_invariants import centroid

def test_packet_continues_after_source_removal():
    x=localized_packet(center=24); c0=centroid(x); s=NativeExcitationState(x); progress_source_free(s,8)
    assert abs(centroid(s.values)-c0-8)<1e-10
    assert dependency_contract()['source_present_after_launch'] is False

