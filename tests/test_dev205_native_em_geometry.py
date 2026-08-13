import numpy as np
from pbuf.observer.native_em_geometry import pair_geometry, native_source_diagnostic


def test_pair_geometry_broadcasts_node_axial_value_to_n6_bonds():
    p=np.zeros((1,1,1,6,3)); p[...,1]=1
    q=np.zeros((1,1,1,3)); q[...,2]=2
    r=pair_geometry(p,q)
    assert r['dot'].shape[-1]==6 and np.allclose(r['cross'][...,0],2)


def test_native_source_diagnostic_is_zero_for_constant_vector():
    assert np.array_equal(native_source_diagnostic(np.ones((3,3,3,3))),np.zeros((3,3,3)))
