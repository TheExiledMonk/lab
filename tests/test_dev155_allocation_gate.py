import numpy as np
from pbuf.excitation.native_excitation_n6 import gaussian_packet,neighbor_mean,quadratic_norm
def test_equal_neighbor_mean_is_not_silently_promoted():
    x=gaussian_packet((12,12,12),center=(4,6,6)); assert not np.isclose(quadratic_norm(neighbor_mean(x)),quadratic_norm(x))
