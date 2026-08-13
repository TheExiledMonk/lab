import numpy as np
from pbuf.observer.native_pattern_boundary import neighbor_mismatch
def test_mismatch_is_directed_neighbor_difference():
    s=np.arange(3*3*3*6,dtype=float).reshape(3,3,3,6); d=neighbor_mismatch(s)
    assert np.array_equal(d[1,1,1,0,:],s[2,1,1]-s[1,1,1])
    assert np.array_equal(d[1,1,1,1,:],s[0,1,1]-s[1,1,1])
