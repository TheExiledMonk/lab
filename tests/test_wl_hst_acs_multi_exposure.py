import numpy as np
from pbuf.wl.hst_acs_geometry import AffineDetectorTransform

def test_same_event_different_exposure_coordinates():
    p=np.array([[1.,2.]])
    a=AffineDetectorTransform(np.eye(2),[0,0],"A");b=AffineDetectorTransform([[0,-1],[1,0]],[0,0],"B")
    assert not np.array_equal(a.forward(p),b.forward(p))
    assert np.allclose(a.reverse(a.forward(p)),b.reverse(b.forward(p)))

