import numpy as np
from pbuf.excitation.native_excitation_state import NativeExcitationState, localized_packet
from pbuf.excitation.native_excitation_transfer import progress_source_free

def test_next_state_depends_on_previous_state():
    x=localized_packet(); s=NativeExcitationState(x.copy()); progress_source_free(s,2)
    np.testing.assert_array_equal(s.values,np.roll(x,2,axis=0)); assert len(s.history)==3

