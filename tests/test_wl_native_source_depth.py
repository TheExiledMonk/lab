import numpy as np
from pbuf.wl.native_source_depth import *
def test_registry_and_simple_focus():
 assert len(estimator_registry())==35
 z=np.linspace(0,4,256);base=np.array([[-1.,0],[1,0],[0,-1],[0,1]])
 p=np.array([(q-2)*base for q in z]);curves=primitive_score_curves(z,p)
 assert abs(candidates_from_curve(z,curves["D03"],"D03")[0].depth_native-2)<.02
def test_secondary_minima_preserved():
 z=np.linspace(0,4,401);q=((z-1)*(z-3))**2
 assert len(candidates_from_curve(z,q,"D03",.001))==2
