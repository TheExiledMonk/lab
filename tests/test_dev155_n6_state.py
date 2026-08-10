import numpy as np
from pbuf.excitation.native_excitation_n6 import NativeExcitationN6State,gaussian_packet
def test_rank2_state_lives_on_3d_grid():
    x=gaussian_packet((9,8,7)); s=NativeExcitationN6State(x)
    assert s.values.shape==(9,8,7,2) and s.rank==2 and s.topology=="N6_3D_PERIODIC"
