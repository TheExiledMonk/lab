import numpy as np
from pbuf.observer.native_em_candidates import axial_dual, transverse_relative_motion


def test_axial_dual_is_fixed_antisymmetric_tensor_identity():
    a=np.zeros((3,3)); a[1,2]=4; a[2,1]=-4
    assert np.allclose(axial_dual(a),[-4,0,0])


def test_transverse_candidate_has_no_fixed_axis_component():
    x=np.array([1.,2.,3.])
    assert np.allclose(transverse_relative_motion(x,[1,0,0]),[0,2,3])
