import numpy as np
from pbuf.observer.native_pair_partition import derive_partition
from pbuf.observer.native_pair_interface_force import interface_bonds, transfer

def test_direct_interface_has_exact_reciprocal_transfer():
    p=derive_partition((11,11,11),(1,5,5),(1,9,5)); b=interface_bonds(p.omega_a,p.omega_b)
    u=np.random.default_rng(217).normal(0,1e-4,(11,11,11,3))
    fa,fb=transfer(u,b)
    assert np.array_equal(fa,-fb)
