import numpy as np
from pbuf.observer.native_pair_partition import derive_partition

def test_frozen_n6_voronoi_partition_is_complete_and_disjoint():
    p=derive_partition((11,11,11),(1,5,5),(1,9,5))
    assert not np.any(p.omega_a & p.omega_b)
    assert np.array_equal(p.omega_a | p.omega_b | p.omega_i,p.omega_d)
    assert (p.omega_a.sum(),p.omega_b.sum(),p.omega_i.sum()) == (605,605,121)
    assert np.all(p.distance_a[p.omega_i] == p.distance_b[p.omega_i])
