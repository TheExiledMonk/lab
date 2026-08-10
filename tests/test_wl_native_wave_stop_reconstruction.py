import numpy as np
from pbuf.wl.native_redshift_stopping import footprint,stopping_depths

def test_known_and_synthetic_wave_stop_footprints():
    src=np.array([[-1,0],[1,0],[0,.5],[0,-.5]],float); depth=3
    received=src+np.array([2.,-1.])*depth
    known=received-np.array([2.,-1.])*depth
    stop=stopping_depths(np.linspace(0,5,501),np.linspace(0,1,501),.6)["stop_candidates"][0]
    wave=received-np.array([2.,-1.])*stop
    assert footprint(known)["RMS_radius"] == footprint(src)["RMS_radius"]
    assert abs(footprint(wave)["RMS_radius"]-footprint(src)["RMS_radius"])<1e-12
