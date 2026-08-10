import numpy as np
import pytest
from pbuf.matter.native_loaded_excitation_modes import center_diagnostics, gaussian_packet, progression_ratio

def test_packet_centroid_and_measured_progression_ratio():
    x=np.linspace(-5,5,101); p=gaussian_packet(x,center=1,width=.5)
    assert abs(center_diagnostics(x,p)['C03']-1)<1e-8
    assert progression_ratio(2,1)==.5
    with pytest.raises(ValueError): progression_ratio(1,1.1)

