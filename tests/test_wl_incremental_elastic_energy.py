import numpy as np
from pbuf.wl.native_incremental_elastic_energy import *

def test_nonzero_background_exact_subtraction():
    bg=np.array([.2,-.15]); de=np.array([.01,.025])
    assert np.allclose(incremental_elastic_energy(bg,de),bounded_strain_energy(bg+de)-bounded_strain_energy(bg))

def test_small_strain_taylor_and_full_law():
    bg=.2; de=1e-5
    exact=incremental_elastic_energy(bg,de)
    assert np.isclose(exact,taylor_increment(bg,de),rtol=1e-8)

def test_signed_increment_and_positive_excitation_are_distinct():
    audit=positivity_audit(.4,np.array([-.1,.1]))
    assert audit["signed_increment_can_be_negative"]
    assert audit["positive_excitation_nonnegative"]

