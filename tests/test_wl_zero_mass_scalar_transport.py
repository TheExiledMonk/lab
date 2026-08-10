import numpy as np
from pbuf.wl.native_zero_mass_scalar_transport import CANDIDATES, apply_factors, apply_state_ratio

def test_all_u_candidates_registered():
    assert set(CANDIDATES) == {f"U{i:02d}" for i in range(1,21)}

def test_identity_exact():
    for q in (.25,.5,1,2,4,8): assert apply_factors(q,np.ones(20)).q_receive == q

def test_telescoping_ratio():
    x=np.array([.5,2,1,8,4.]); s=apply_state_ratio(3,x)
    assert np.isclose(s.q_ratio,x[-1]/x[0])
