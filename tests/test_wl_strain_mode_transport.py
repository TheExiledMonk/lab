import numpy as np
from pbuf.wl.native_strain_mode_transport import *

def test_direction_only_bending():
    t=np.linspace(0,1,20); d=np.c_[np.cos(t),np.sin(t)]
    r=direction_scalar_control(d,np.ones(20))
    assert r["direction_changed"] and r["direction_independent"]

def test_entry_exit_and_reverse():
    assert entry_exit_control(1,.8,1)["classification"]=="TEMPORARY_EXCHANGE"
    assert forward_reverse_control(2,2)["conservative"]

