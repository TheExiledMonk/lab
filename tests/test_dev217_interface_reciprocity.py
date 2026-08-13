import numpy as np
from pbuf.observer.native_pair_partition import derive_partition
from pbuf.observer.native_pair_interface_force import interface_bonds

def test_direct_periodic_ab_edge_plane_is_canonical_and_not_double_counted():
    p=derive_partition((11,11,11),(1,5,5),(1,9,5)); b=interface_bonds(p.omega_a,p.omega_b)
    assert len(b['axis']) == 121
    assert np.all(b['axis'] == 1) and np.all(b['orientation'] == -1)
    assert len({tuple(x) for x in np.c_[b['node_a'],b['node_b']]}) == 121
