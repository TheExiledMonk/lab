import numpy as np
from pbuf.observer.native_spatial_winding import reflect_x

def test_x_reflection_is_an_involution_with_polar_vector_rule():
    u = np.arange(11 * 11 * 11 * 3, dtype=float).reshape(11, 11, 11, 3)
    assert np.array_equal(reflect_x(reflect_x(u, 1), 1), u)
