import numpy as np
from pbuf.analysis.native_staggered_order import sign_counts

def test_axis_counts_close_to_axis_bond_count():
    c=sign_counts(np.array([[-1,1,0,-1,1,0]]),np.array([0,0,1,1,2,2]))
    assert all(v['opposed'][0]+v['aligned'][0]+v['zero'][0] == v['total'][0] for v in c['by_axis'].values())
