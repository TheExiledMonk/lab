import numpy as np
from pbuf.foundation.fast_slow_excitation_audit import classify_pair_transfer, persistence_test, terminal_pair_transfer

def test_pair_is_static_and_antisymmetric():
    a=np.array([1.,-2.]); b=np.array([.5,3.])
    np.testing.assert_allclose(terminal_pair_transfer(-a,-b),-terminal_pair_transfer(a,b))
    assert classify_pair_transfer(a,b)['candidate_status']=='STATIC_MEDIUM_STATE'
    assert persistence_test(a)['classification']=='STATIC_MEDIUM_SAMPLE'

