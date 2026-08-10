import numpy as np
from pbuf.wl.native_source_reconstruction import *
def test_forward_constrained_depth_direction_advantage():
 src=np.array([[-1.,0],[1.,0],[0.,1]])
 def f(p,z):return p+z*np.array([.2,-.1]),np.tile([z*.1,1.],(len(p),1))
 obs,dirs=f(src,2.)
 def builder(z):return obs-z*np.array([.2,-.1])
 r,q=forward_constrained_inverse(np.linspace(1,3,81),obs,builder,f,received_directions=dirs)
 assert abs(r.source_depth_native-2)<1e-12 and r.roundtrip_score<1e-12
