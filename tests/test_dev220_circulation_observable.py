import numpy as np
from pbuf.observer.native_spatial_winding import circulation, square_yz_loop

def test_zero_bond_is_undefined_without_threshold():
    u = np.zeros((11, 11, 11, 3))
    result = circulation(u, square_yz_loop((11, 11, 11), (1, 5, 5), 3))
    assert result['defined'] is False
