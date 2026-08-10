import numpy as np
from pbuf.excitation.native_excitation_n6 import N6_OFFSETS,neighbor_stack
def test_all_six_neighbors_are_executed():
    x=np.arange(7*8*9*2.).reshape(7,8,9,2); y=neighbor_stack(x)
    assert len(N6_OFFSETS)==6 and y.shape==(6,7,8,9,2)
