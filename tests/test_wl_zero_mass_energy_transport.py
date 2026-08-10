import numpy as np
from pbuf.wl.native_zero_mass_energy import propagate_energy_ratio,direction_magnitude_audit

def test_uniform_and_forward_reverse_controls():
    s=np.linspace(0,2,101); q=np.sin(np.pi*s)
    assert np.allclose(propagate_energy_ratio(s,np.zeros_like(s)),1)
    f=propagate_energy_ratio(s,q)[-1]
    r=propagate_energy_ratio(s[::-1],q[::-1])[-1]
    assert np.isclose(f*r,1)

def test_curved_direction_is_not_magnitude_change():
    t=np.linspace(0,1,10); out=direction_magnitude_audit(np.c_[np.cos(t),np.sin(t)])
    assert out["MOMENTUM_DIRECTION_CHANGED"] and not out["MOMENTUM_MAGNITUDE_CHANGED"]
