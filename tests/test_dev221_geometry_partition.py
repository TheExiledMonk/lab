import numpy as np
from pbuf.observer.native_extended_geometry import signed_axis_coordinates
def test_partition_has_explicit_center_plane():
 s=signed_axis_coordinates(11,1); assert (s==0).sum()==1 and (s>0).sum()==5 and (s<0).sum()==5
