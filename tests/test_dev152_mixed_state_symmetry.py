import numpy as np
from pbuf.foundation.native_neighbor_mixed_state import construct_case, progress_case

def test_excitation_sign_reversal():
    a=construct_case(2,2); b={**a,"X":-a["X"]}
    assert np.allclose(progress_case(a)["history"][...,1:],-progress_case(b)["history"][...,1:])

