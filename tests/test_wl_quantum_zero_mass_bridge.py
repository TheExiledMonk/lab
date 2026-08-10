import numpy as np
from pbuf.wl.quantum_zero_mass_bridge import *

def test_absolute_external_relations_are_consistent():
    p=2.5e-27; k=wave_number_from_momentum(p); e=energy_from_momentum(p)
    assert np.isclose(e,HBAR_SI*C_SI*k)
    assert np.isclose(wave_number_from_energy(e),k)
    assert np.isclose(wavelength_from_energy(e),2*np.pi/k)

def test_ratio_cancellation():
    r=np.array([1,.99,.5,.2]); b=ratio_bridge(r)
    assert np.array_equal(b["momentum_ratio"],r)
    assert np.array_equal(b["k_ratio"],r)
    assert np.allclose(b["wavelength_ratio"],1/r)
    assert np.allclose(b["one_plus_z"],1/r)
