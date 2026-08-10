import numpy as np
from pbuf.wl.native_zero_mass_energy import scale_cancellation
from pbuf.wl.native_energy_redshift import redshift_from_energy_ratio

def test_L0_and_absolute_energy_cancel():
    assert all(scale_cancellation(.2,.3,a)["exact_cancellation"] for a in (.5,1,2,4))
    r=np.array([1,.9,.5]); beta=7.3
    assert np.allclose(redshift_from_energy_ratio((beta*r)/beta),redshift_from_energy_ratio(r))
