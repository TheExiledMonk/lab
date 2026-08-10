import numpy as np
from pbuf.foundation.trajectory_normalization_audit import decompose_update

def test_raw_normalized_identity():
    out=decompose_update(np.array([[0.,0.,1.]]),np.array([[2.,0.,0.]]),.1)
    np.testing.assert_allclose(out['normalized_vector'],out['raw_vector']/out['raw_magnitude'][:,None])
    np.testing.assert_allclose(np.linalg.norm(out['normalized_vector'],axis=1),1)

