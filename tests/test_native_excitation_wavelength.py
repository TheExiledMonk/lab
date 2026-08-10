import numpy as np
from pbuf.quantum.native_excitation_modes import carrier_mode,estimate_wavelengths,native_k
def test_independent_wavelength_estimators_agree():
    e=estimate_wavelengths(carrier_mode(192,16)); assert np.isclose(e['L01'],16,rtol=.1); assert np.isclose(e['L03'],16,rtol=.1); assert np.isclose(e['L04'],16)
def test_k_requires_positive_wavelength():
    assert np.isclose(native_k(8),np.pi/4)

