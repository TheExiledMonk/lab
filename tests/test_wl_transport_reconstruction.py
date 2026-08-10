import numpy as np
from pbuf.wl.reconstructed_geometry import geometry_from_derivatives
from pbuf.wl.transport_reconstruction import quadratic_taylor_error

def banks(J,shape=(7,9),curvature=False):
 z=np.zeros(shape);f={"d_u_delta_u":z+(J[0,0]-1),"d_v_delta_u":z+J[0,1],"d_u_delta_v":z+J[1,0],"d_v_delta_v":z+(J[1,1]-1)}
 s={f"d_{a}_{b}":z.copy() for b in ("delta_u","delta_v","wf") for a in ("uu","uv","vv")} if curvature else None
 return f,s

def test_identity_transport():
 f,s=banks(np.eye(2),curvature=True);g=geometry_from_derivatives(f,s)
 assert np.allclose(g["local_area_change"],0) and np.allclose(g["orientation_change"],0)
 assert np.allclose(g["local_curvature"],0) and np.allclose(g["local_anisotropy"],0)

def test_translation_has_identity_derivative_geometry():
 g=geometry_from_derivatives(banks(np.eye(2))[0]);assert np.allclose(g["transport_area_ratio"],1);assert np.allclose(g["spin2_shape_q1"],0)

def test_isotropic_scale():
 s=1.7;g=geometry_from_derivatives(banks(np.eye(2)*s)[0]);assert np.allclose(g["transport_area_ratio"],s*s);assert np.allclose(g["local_anisotropy"],0)

def test_anisotropic_stretch():
 g=geometry_from_derivatives(banks(np.diag([1.4,.8]))[0]);assert np.all(g["spin2_shape_q1"]>0)

def test_rotation():
 a=.37;J=np.array([[np.cos(a),-np.sin(a)],[np.sin(a),np.cos(a)]]);g=geometry_from_derivatives(banks(J)[0]);assert np.allclose(g["transport_area_ratio"],1);assert np.allclose(g["local_anisotropy"],0,atol=1e-15);assert np.allclose(g["orientation_change"],a)

def test_quadratic_second_order_is_exact_and_better():
 first,second=quadratic_taylor_error();assert second<1e-14;assert second<first
