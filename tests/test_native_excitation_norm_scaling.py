import numpy as np
from pbuf.quantum.native_excitation_modes import carrier_mode,quadratic_norm
def test_quadratic_amplitude_scaling():
    q=carrier_mode(128,16); assert np.isclose(quadratic_norm(.25*q),.25**2*quadratic_norm(q))

