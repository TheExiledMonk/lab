import numpy as np
from pbuf.wl.hst_acs_geometry import AffineDetectorTransform

def test_direction_independent_affine_reverse_is_float64():
    t=AffineDetectorTransform([[2.,.1],[-.2,3.]],[11.,-7.])
    p=np.linspace(-2,2,100).reshape(50,2);q=t.reverse(t.forward(p))
    assert q.dtype==np.float64 and np.max(np.abs(q-p))<2e-15

