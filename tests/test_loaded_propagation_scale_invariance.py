import numpy as np
from pbuf.matter.native_mass_loading_state import strain_loading_fraction

def test_dimensionless_strain_loading_is_coordinate_scale_invariant():
    e=np.array([.1,.25,.75])
    baseline=strain_loading_fraction(e)
    for alpha in (.5,1,2,4): np.testing.assert_allclose(strain_loading_fraction(e),baseline)
