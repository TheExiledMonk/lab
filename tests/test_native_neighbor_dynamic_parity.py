import numpy as np
from pbuf.foundation.native_neighbor_state import NativeNeighborState
from pbuf.foundation.native_neighbor_dynamic_projection import dynamic_parity
def test_exact_shift_parity():
    rng=np.random.default_rng(3); assert dynamic_parity(NativeNeighborState(np.zeros(32),rng.normal(size=(32,2))),7)['status']=='PARITY_ESTABLISHED'
