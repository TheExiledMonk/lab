from pbuf.matter.native_excitation_state import NativeExcitationState, energy_like_classification, excitation_registry

def test_x01_x20_attempted_and_q_remains_neutral():
    rows=excitation_registry(); assert [r['id'] for r in rows]==[f'X{i:02d}' for i in range(1,21)]
    s=NativeExcitationState(2); s.identity_progress(.5)
    assert s.q_state==2 and s.history[0]['spatial_position']==.5
    assert not energy_like_classification()['EXCITATION_ENERGY_LIKE_STATE_ESTABLISHED']

