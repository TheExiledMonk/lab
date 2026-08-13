import numpy as np

from pbuf.observer.native_two_strain_magnetism import support_mask


def test_dev211_support_is_exact_six_neighbor_contact():
    assert support_mask((13, 13, 13), (3, 6, 6)).sum() == 6
