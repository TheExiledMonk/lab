import numpy as np
from pbuf.analysis.native_staggered_order import unique_n6_bonds, contract_bonds

def test_simultaneous_tensor_reflection_preserves_contraction():
    a=np.random.default_rng(3).normal(size=(1,2,2,2,3,3)); p,_,_=unique_n6_bonds((2,2,2)); q=np.diag([-1.,1.,1.])
    assert np.allclose(contract_bonds(a,p),contract_bonds(np.einsum('ij,t...jk,kl->t...il',q,a,q.T),p))
