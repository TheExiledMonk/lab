import numpy as np
from pbuf.foundation.native_link_constitutive_derivatives import tangent_stiffness, constitutive_curvature, validate_analytic_derivatives
def test_exact_derivatives():
    e=np.linspace(-.9,.9,101); assert np.array_equal(tangent_stiffness(e),constitutive_curvature(e))
    result=validate_analytic_derivatives(e); assert result["tangent_valid"] and result["curvature_valid"]
