import numpy as np
from pbuf.foundation.native_neighbor_mixed_state import construct_case, progress_case

def test_load00_preserves_dynamic_norm_and_static_state():
    c=construct_case(0,0); r=progress_case(c)
    assert r["relative_norm_drift"] < 1e-12
    assert np.array_equal(r["history"][0,:,0],r["history"][-1,:,0])

