import numpy as np
from pbuf.foundation.trajectory_normalization_audit import decompose_update

def test_raw_magnitude_is_path_step_dependent():
    n=np.array([0.,0.,1.]); r=np.array([1.,0.,0.])
    a=decompose_update(n,r,.5)['raw_magnitude']; b=decompose_update(n,r,2.)['raw_magnitude']
    assert b>a

