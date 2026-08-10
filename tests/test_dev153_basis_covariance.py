import numpy as np
from pbuf.foundation.native_transverse_transfer_capacity import orthogonal_transport
def test_frame_transport_covariance():
    B=orthogonal_transport(.7); R=orthogonal_transport(.2); x=np.array([1.,2.])
    assert np.allclose(B@(R@x),(B@R@B.T)@(B@x))
