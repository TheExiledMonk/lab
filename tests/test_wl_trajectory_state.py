import numpy as np

from pbuf.wl.trajectory_state import (
    PATH_FRACTIONS, TrajectoryAccumulator, bundle_history,
    sample_path_fractions, summarize_trajectory,
)


def test_straight_path_and_known_polyline():
    x=np.column_stack((np.zeros(11),np.zeros(11),np.linspace(0,3,11)))
    d=np.tile([0.,0.,1.],(11,1));s,_=summarize_trajectory(x,d)
    assert s["path_length"] == 3
    assert s["path_excess"] < 1e-14
    assert s["net_direction_change"] == s["total_direction_change"] == 0
    assert s["path_curvature_integral"] == s["curvature_max"] == 0
    q=np.array([[0.,0.,0.],[3.,0.,0.],[3.,4.,0.]])
    dirs=np.array([[1.,0.,0.],[1.,0.,0.],[0.,1.,0.]])
    assert summarize_trajectory(q,dirs)[0]["path_length"] == 7


def test_circular_arc_has_known_curvature_and_turn():
    radius=4.; angle=np.linspace(0,np.pi/2,2001)
    x=np.column_stack((radius*np.cos(angle),radius*np.sin(angle),np.zeros_like(angle)))
    d=np.column_stack((-np.sin(angle),np.cos(angle),np.zeros_like(angle)))
    s,_=summarize_trajectory(x,d)
    np.testing.assert_allclose(s["total_direction_change"],np.pi/2,rtol=2e-8)
    np.testing.assert_allclose(s["curvature_mean"],1/radius,rtol=2e-7)
    np.testing.assert_allclose(s["total_direction_change"],s["path_length"]/radius,rtol=2e-7)


def test_s_curve_retains_variation_hidden_by_endpoint():
    t=np.linspace(0,1,1001);phase=.45*np.sin(2*np.pi*t)
    d=np.column_stack((np.sin(phase),np.zeros_like(t),np.cos(phase)))
    dt=np.diff(t);x=np.vstack(([0.,0.,0.],np.cumsum(.5*(d[:-1]+d[1:])*dt[:,None],axis=0)))
    s,_=summarize_trajectory(x,d)
    assert s["net_direction_change"] < 1e-12
    assert s["total_direction_change"] > 1


def _bundle(mapper):
    y,x=np.mgrid[-1:1:9j,-1:1:9j];uv=np.stack((x,y),axis=-1);states=[]
    for f in PATH_FRACTIONS:
        u,v=mapper(f,x,y);states.append(np.stack((u,v,np.full_like(u,f)),axis=-1))
    return bundle_history(np.asarray(states),uv)


def test_bundle_expansion_anisotropy_and_rotation():
    z=_bundle(lambda f,x,y:((1+f)*x,(1+f)*y))
    np.testing.assert_allclose(z["area_ratio"][-1,4,4],4,atol=1e-12)
    z=_bundle(lambda f,x,y:((1+f)*x,(1+.5*f)*y))
    np.testing.assert_allclose(z["singular_value_1"][-1,4,4],2,atol=1e-12)
    np.testing.assert_allclose(z["singular_value_2"][-1,4,4],1.5,atol=1e-12)
    z=_bundle(lambda f,x,y:(np.cos(f)*x-np.sin(f)*y,np.sin(f)*x+np.cos(f)*y))
    np.testing.assert_allclose(z["determinant"][:,4,4],1,atol=1e-12)
    np.testing.assert_allclose(z["singular_value_1"][:,4,4],1,atol=1e-12)


def test_quadratic_second_order_and_determinism():
    z=_bundle(lambda f,x,y:(x+f*(2*x*x+3*x*y+4*y*y),y))
    np.testing.assert_allclose(z["hessian"][-1,4,4,0],[4,3,8],atol=2e-12)
    q=_bundle(lambda f,x,y:(x+f*x*x,y))
    for key in q: np.testing.assert_equal(q[key],_bundle(lambda f,x,y:(x+f*x*x,y))[key])


def test_accumulator_native_integrals_and_fraction_sampling():
    a=TrajectoryAccumulator()
    for z in range(3):a.update([0,0,z],[0,0,1],{"response":z})
    r=a.finalize();assert r.native_path_summary["response"]["integral"] == 2
    sampled=sample_path_fractions([[0,0,0],[0,0,2]],[[0,0,1],[0,0,1]])
    np.testing.assert_allclose(sampled["positions"][:,2],2*PATH_FRACTIONS)
