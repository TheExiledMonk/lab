import numpy as np
from pbuf.observer.native_internal_state_inventory import axial_dual


def test_axial_dual_reverses_with_antisymmetric_tensor():
    a = np.zeros((3, 3)); a[0, 1] = 2; a[1, 0] = -2
    assert np.array_equal(axial_dual(-a), -axial_dual(a))
