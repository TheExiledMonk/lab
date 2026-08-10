import numpy as np
from pbuf.wl.multiscale_transport_relations import (canonical_manifest, derivatives,
    local_mean_variance, spatial_quadrupole, spin2_rotate)


def grid(n=65):
    y,x=np.mgrid[-(n//2):n//2+1, -(n//2):n//2+1]
    return x.astype(float),y.astype(float)


def test_manifest_is_frozen_rank_29():
    m=canonical_manifest(); assert len(m)==29; assert [r["index"] for r in m]==list(range(29))


def test_constant_field():
    z=np.ones((65,65))*7; mean,var,valid=local_mean_variance(z,4); d=derivatives(mean); q=spatial_quadrupole(z,4)
    assert np.all(var[valid] == 0); assert np.allclose(d["gradient_magnitude"][8:-8,8:-8],0)
    assert np.allclose(d["eigenvalue_difference"][8:-8,8:-8],0); assert np.allclose(q["q_abs"][valid],0)


def test_linear_gradient_and_zero_curvature():
    x,y=grid(); z=2*x-3*y; d=derivatives(z)
    assert np.allclose(d["gradient_u"],2); assert np.allclose(d["gradient_v"],-3)
    assert np.allclose(d["trace"],0); assert np.allclose(d["eigenvalue_difference"],0)


def test_isotropic_and_anisotropic_quadratics():
    x,y=grid(); iso=derivatives(x*x+y*y); an=derivatives(x*x-y*y); rot=derivatives(2*x*y)
    c=(32,32); assert abs(iso["trace"][c])>0 and abs(iso["hessian_q1"][c])<1e-14 and abs(iso["hessian_q2"][c])<1e-14
    assert abs(an["hessian_q1"][c]-4)<1e-14 and abs(an["hessian_q2"][c])<1e-14
    assert abs(rot["hessian_q1"][c])<1e-14 and abs(rot["hessian_q2"][c]-4)<1e-14


def test_spin2_rotation_exact():
    q1,q2=0.37,-1.2; phi=.317; a,b=spin2_rotate(q1,q2,phi)
    expected=(q1+1j*q2)*np.exp(2j*phi)
    assert abs((a+1j*b)-expected)/(abs(expected)+np.finfo(float).eps)<=1e-12


def test_scale_sanity_constant_and_linear():
    x,y=grid();
    for z in (np.ones_like(x),2*x-y):
        for r in (2,4,8):
            q=spatial_quadrupole(z,r); assert np.allclose(q["q_abs"][r:-r,r:-r],0,atol=1e-12)
