import numpy as np
from pbuf.foundation.native_neighbor_state import NativeNeighborState
from pbuf.foundation.native_neighbor_static_projection import static_parity
def test_frozen_bounded_law_parity():
    assert static_parity(NativeNeighborState(np.linspace(-.8,.8,21),np.zeros((21,2))))['status']=='PARITY_ESTABLISHED'
