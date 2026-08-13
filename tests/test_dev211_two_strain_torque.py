import numpy as np

from pbuf.observer.native_two_strain_magnetism import torque_support


def test_dev211_support_torque_is_cross_product_sum():
    m = np.zeros((5, 5, 5), dtype=bool); m[2, 3, 2] = True
    f = np.zeros((5, 5, 5, 3)); f[2, 3, 2] = (1, 0, 0)
    assert np.array_equal(torque_support(f, m, (2, 2, 2)), np.array([0., 0., -1.]))
