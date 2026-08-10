import numpy as np
import pytest
from pbuf.wl.hst_acs_geometry import AffineDetectorTransform, RectangularChip, classify_detector_points

def test_rotation_roundtrip_and_reflection_preserved():
    p=np.array([[1.,2.],[-3.,4.]])
    for degrees in (0,90,180,270):
        a=np.deg2rad(degrees);m=np.array([[np.cos(a),-np.sin(a)],[np.sin(a),np.cos(a)]])
        t=AffineDetectorTransform(m,np.array([4.,5.]));assert np.allclose(t.reverse(t.forward(p)),p,atol=1e-14)
    assert np.linalg.det(AffineDetectorTransform(np.diag([-1.,1.]),np.zeros(2)).matrix)<0

def test_singular_inverse_rejected():
    with pytest.raises(ValueError,match="NOT_REVERSIBLE"):
        AffineDetectorTransform([[1,0],[0,0]],[0,0]).reverse([[1,1]])

def test_two_chip_gap_and_boundaries():
    chips=(RectangularChip("CHIP_1",0,10,0,4),RectangularChip("CHIP_2",0,10,6,10))
    status,ids=classify_detector_points([[5,2],[5,8],[5,5],[-1,2],[0,2]],chips)
    assert status.tolist()==["ACTIVE_CHIP","ACTIVE_CHIP","INTER_CHIP_GAP","OUTSIDE_DETECTOR","BOUNDARY"]
    assert ids[:2].tolist()==["CHIP_1","CHIP_2"]

