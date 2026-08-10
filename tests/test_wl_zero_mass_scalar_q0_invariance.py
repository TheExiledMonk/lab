import numpy as np
from pbuf.wl.native_scalar_transport_audit import q0_fractional_invariance
from pbuf.wl.native_zero_mass_scalar_transport import transport_on_frozen_trajectory

def test_fractional_multiplicative_law_is_q0_independent():
    assert q0_fractional_invariance([1.1,.8,1.2])["max_cv"] < 1e-14

def test_scalar_does_not_change_frozen_trajectory():
    p=np.c_[np.arange(5),np.zeros(5)]; d=np.tile([1.,0.],(5,1))
    p1,d1,_=transport_on_frozen_trajectory(1,p,d,np.ones(4))
    p4,d4,_=transport_on_frozen_trajectory(4,p,d,np.ones(4))
    assert np.array_equal(p1,p4) and np.array_equal(d1,d4)
