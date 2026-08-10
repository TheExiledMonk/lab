import numpy as np
from pbuf.foundation.native_link_constitutive_derivatives import descriptors
def test_descriptor_relations_are_algebraic():
    e=np.linspace(0,.9,50); d=descriptors(e)
    assert np.allclose(d["link_stretch_ratio"],1+e) and np.all(np.diff(d["stress"])>0)
    assert np.array_equal(d["tangent_stiffness"],d["constitutive_curvature"])
