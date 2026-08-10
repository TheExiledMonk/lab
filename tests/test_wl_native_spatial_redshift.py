import numpy as np
from pbuf.wl.native_spatial_redshift import stopping_candidates,spatial_redshift_stop

def test_all_nonmonotonic_roots_preserved():
    s=np.linspace(0,4,401); z=np.sin(np.pi*s)
    roots=stopping_candidates(s,z,.5); assert len(roots)==4
    assert spatial_redshift_stop(s,z,.5,mechanism="synthetic",scale_free=True)["ambiguity"]=="MULTIPLE_SPATIAL_REDSHIFT_STOPS"
