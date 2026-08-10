import numpy as np
from pbuf.wl.optical_bundle_transport import compose_transport, derivative_invariants, finite_difference_derivatives


def test_combined_jacobian_order_and_incompatibility():
    a=np.array([[2.,1.],[0.,3.]]);b=np.array([[1.,4.],[2.,0.]])
    total,status=compose_transport(a,b)
    assert status=="TRANSPORT_COMPOSITION_VALID" and np.array_equal(total,a@b)
    assert compose_transport(a,b,compatible=False)==(None,"TRANSPORT_COMPOSITION_NOT_VALID")


def test_affine_jacobian_and_quadratic_hessian():
    def f(x): return np.column_stack((2*x[:,0]+3*x[:,1]+x[:,0]**2, -x[:,0]+4*x[:,1]+2*x[:,1]**2))
    j,h=finite_difference_derivatives(f,np.array([[0.,0.]]),step=1e-4)
    assert np.allclose(j[0],[[2,3],[-1,4]],atol=1e-9)
    assert np.allclose(h[0,0],[[2,0],[0,0]],atol=1e-7)
    assert np.allclose(h[0,1],[[0,0],[0,4]],atol=1e-7)


def test_singular_classification():
    x=derivative_invariants(np.array([[[1.,0.],[0.,0.]]]))
    assert x["classification"][0]=="LOCALLY_SINGULAR_CANDIDATE"
    assert np.isinf(x["condition_number"][0])
