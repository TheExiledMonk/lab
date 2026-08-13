import numpy as np
from pbuf.observer.native_mode_classification import longitudinal_transverse, transverse_rank


def test_exact_l_t_decomposition():
    k=np.array([1.,0.,0.])
    assert longitudinal_transverse(np.array([1.,0.,0.]),k)['sector']=='EXACT_LONGITUDINAL'
    assert longitudinal_transverse(np.array([0.,1.,0.]),k)['sector']=='EXACT_TRANSVERSE'


def test_transverse_rank_is_algebraic():
    assert transverse_rank(np.eye(3),np.array([1.,0.,0.]))==2
