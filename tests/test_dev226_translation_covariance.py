import numpy as np
from pbuf.analysis.native_staggered_order import unique_n6_bonds, contract_bonds

def test_translation_preserves_periodic_bond_values_as_a_multiset():
    a=np.random.default_rng(1).normal(size=(1,3,3,3,3,3)); p,_,_=unique_n6_bonds((3,3,3))
    assert np.array_equal(np.sort(contract_bonds(a,p)),np.sort(contract_bonds(np.roll(a,1,axis=1),p)))
