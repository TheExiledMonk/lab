import numpy as np
from pbuf.analysis.native_staggered_order import axial_equivalence_error, contract_bonds

def test_axial_identity_matches_tensor_contraction():
    a=np.zeros((1,2,1,1,3,3)); a[0,0,0,0,0,1]=2; a[0,0,0,0,1,0]=-2; a[0,1,0,0,0,1]=3; a[0,1,0,0,1,0]=-3
    p=np.array([[[0,0,0],[1,0,0]]]); assert axial_equivalence_error(a,p,contract_bonds(a,p)) == 0
