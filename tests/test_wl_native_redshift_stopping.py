import numpy as np
from pbuf.wl.native_redshift_stopping import *

def test_monotonic_multiple_and_no_stop():
    s=np.linspace(0,4,401)
    assert abs(stopping_depths(s,s/4,.5)["stop_candidates"][0]-2)<1e-12
    z=(s-1)*(s-2)*(s-3)
    r=stopping_depths(s,z,0)
    assert r["status"] == "MULTIPLE_REDSHIFT_STOP_CANDIDATES" and len(r["stop_candidates"]) == 3
    assert stopping_depths(s,s,9)["status"] == "NO_REDSHIFT_STOP_SOLUTION"
