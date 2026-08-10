import numpy as np
from pbuf.wl.native_energy_redshift import *

def test_ratio_redshift_roundtrip_and_stops():
    r=np.array([1,.9,.8]); assert np.allclose(energy_ratio_from_redshift(redshift_from_energy_ratio(r)),r)
    s=np.linspace(0,4,401); ratio=1/(1+np.sin(np.pi*s/2)**2)
    out=energy_redshift_stop(s,ratio,.5)
    assert out["classification"]=="MULTIPLE_STOPS" and len(out["stop_candidates"])==4

def test_no_stop():
    s=np.arange(3.); out=energy_redshift_stop(s,np.ones(3),1)
    assert out["ambiguity"]=="NO_ENERGY_REDSHIFT_STOP"
