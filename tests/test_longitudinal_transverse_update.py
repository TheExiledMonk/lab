import numpy as np
from pbuf.foundation.trajectory_normalization_audit import decompose_update

def test_decomposition_reconstructs_delta():
    out=decompose_update(np.array([0.,0.,1.]),np.array([2.,3.,4.]),.25)
    np.testing.assert_allclose(out['parallel_vector']+out['transverse_vector'],out['delta_vector'])
    np.testing.assert_allclose(np.dot(out['transverse_vector'],np.array([0.,0.,1.])),0,atol=1e-15)

