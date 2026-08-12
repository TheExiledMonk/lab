import numpy as np
from pbuf.observer.native_pair_interaction import reflected_x, orientation_packets, four_state_trajectory

def test_reflection_is_an_exact_native_lattice_operation():
    a=np.arange(4*3*3*3.,dtype=float).reshape(4,3,3,3)
    q=reflected_x(a)
    assert np.array_equal(reflected_x(q),a)
    assert np.array_equal(q[...,1:],np.flip(a,axis=0)[...,1:])

def test_four_state_state_residual_is_inclusion_exclusion():
    z=lambda x:{'displacement':x,'momentum':x,'force':x}
    x=np.zeros((1,3,3,3,3)); a=.001*np.ones_like(x); b=2*a
    r=four_state_trajectory(z(x),z(a),z(b),z(a+b))
    assert np.array_equal(r['momentum'],x)
