import numpy as np
from pbuf.wl.native_zero_mass_strain_mode import *

def test_decomposition_and_identity_are_deterministic():
    bg=np.full(32,.1); de=np.exp(-((np.arange(32)-16)/3)**2)*.05
    got_bg,got_de=decompose_strain(bg+de,bg)
    assert np.allclose(got_bg,bg) and np.allclose(got_de,de)
    a=extract_packet(de,bg,event_uid="e",trajectory_uid="t")
    b=extract_packet(de,bg,event_uid="e",trajectory_uid="t")
    assert a.packet_uid==b.packet_uid and np.array_equal(a.mask,b.mask)

def test_same_direction_different_amplitude_scalar():
    x=np.arange(64); shape=np.exp(-.5*((x-32)/4)**2)
    a=scalar_candidates(extract_packet(.05*shape))["A09"]
    b=scalar_candidates(extract_packet(.10*shape))["A09"]
    assert b>a

