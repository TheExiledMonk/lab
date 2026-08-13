import numpy as np
from pbuf.observer.native_extended_geometry import geometry_moments
def test_shape_tensor_is_symmetric():
 g=geometry_moments(np.ones((11,11,11)),(1,5,5))['shape_tensor']; assert np.array_equal(g,g.T)
