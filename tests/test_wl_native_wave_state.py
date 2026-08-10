from pbuf.wl.native_wave_state import *

def test_frequency_wavelength_and_circularity():
    assert frequency_closure(2,4)["T0"] == .5
    assert wavelength_closure(2,6)["L0"] == 3
    assert frequency_closure(2,4,derived_from_T0=True)["status"] == "CIRCULAR"
    assert wavelength_closure(2,6,derived_from_L0=True)["status"] == "CIRCULAR"

def test_ratio_routes_agree():
    assert abs(redshift_from_wavelength(2,3)-.5)<1e-15
    assert abs(redshift_from_frequency(3,2)-.5)<1e-15
