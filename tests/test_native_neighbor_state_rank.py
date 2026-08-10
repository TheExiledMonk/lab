import numpy as np
from pbuf.foundation.native_neighbor_state import NativeNeighborState,state_registry,local_link_frame,decompose,compose
def test_minimum_rank_and_round_trip():
    f=local_link_frame([1,2,3]); d=np.array([.2,-.3,.4]); q,x=decompose(d,f)
    assert np.allclose(compose(q,x,f),d)
    assert NativeNeighborState(np.zeros(4),np.zeros((4,2))).physical_rank==3
    assert len(state_registry())==20
