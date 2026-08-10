import numpy as np
from pbuf.wl.native_source_reconstruction import ambiguity_area
def test_rich_lane_reduces_ambiguity():
 z,r=np.mgrid[-1:1:101j,-1:1:101j]
 apparent=(z+r)**2;rich=apparent+4*(z-r)**2
 assert ambiguity_area(rich)<ambiguity_area(apparent)
