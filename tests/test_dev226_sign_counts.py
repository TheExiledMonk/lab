import numpy as np
from pbuf.analysis.native_staggered_order import sign_counts

def test_exact_sign_counts_have_no_tolerance():
    c=sign_counts(np.array([[-1.,0.,1.,1e-300]]),np.array([0,0,1,2]))
    assert (c['opposed'][0],c['aligned'][0],c['zero'][0]) == (1,2,1)
