import numpy as np
from pbuf.analysis.native_staggered_order import contract_bonds

def test_tensor_contraction_is_raw_frobenius_product():
    a=np.zeros((1,2,1,1,3,3)); a[0,0,0,0,0,1]=2; a[0,1,0,0,0,1]=-3
    p=np.array([[[0,0,0],[1,0,0]]])
    assert contract_bonds(a,p)[0,0] == -6
