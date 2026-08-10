import numpy as np
import pytest
from pbuf.wl.arrival_formation import (ReceiverPlane, native_receiver_plane,
    intersect_rays, form_arrival_events, arrival_relations, canonical_sha, channel_manifest)


def plane(): return ReceiverPlane(np.zeros(3),np.array([1.,0,0]),np.array([0.,1,0]),np.array([0.,0,1]))


def test_forward_oblique_on_surface_parallel_backward():
    p=plane();x=np.array([[0,0,-1],[1,2,-2],[3,4,0],[0,0,-1],[0,0,1]],float)
    d=np.array([[0,0,1],[1,0,1],[0,0,1],[1,0,0],[0,0,1]],float)
    q,t,_,_,s=intersect_rays(x,d,p)
    assert t[0]==pytest.approx(1);assert np.allclose(q[0,:2],0);assert s[0]=="FORWARD_INTERSECTION"
    assert t[1]==pytest.approx(2*np.sqrt(2));assert np.allclose(q[1,:2],[3,2]);assert s[2]=="ON_SURFACE"
    assert s[3]=="PARALLEL" and np.isnan(t[3]);assert s[4]=="BACKWARD_INTERSECTION" and t[4]<0


def test_invalid_basis_and_handedness():
    with pytest.raises(ValueError):ReceiverPlane(np.zeros(3),np.array([2.,0,0]),np.array([0.,1,0]),np.array([0.,0,1]))
    with pytest.raises(ValueError):ReceiverPlane(np.zeros(3),np.array([1.,0,0]),np.array([0.,-1,0]),np.array([0.,0,1]))


def test_translation_and_basis_rotation_invariance():
    p=plane();x=np.array([[2.,3.,-1.]]);d=np.array([[0.,0.,1.]])
    a=form_arrival_events(x,d,p).event_geometry;shift=np.array([7.,-4.,8.])
    moved=ReceiverPlane(p.origin+shift,p.e_u,p.e_v,p.normal)
    b=form_arrival_events(x+shift,d,moved).event_geometry
    assert a["arrival_u"]==pytest.approx(b["arrival_u"]);assert a["arrival_v"]==pytest.approx(b["arrival_v"])
    c=np.cos(.4);s=np.sin(.4);rot=ReceiverPlane(np.zeros(3),np.array([c,s,0]),np.array([-s,c,0]),p.normal)
    r=form_arrival_events(x,d,rot).event_geometry
    assert np.hypot(a["arrival_u"],a["arrival_v"])==pytest.approx(np.hypot(r["arrival_u"],r["arrival_v"]))


def test_same_position_different_direction_and_latent_reference():
    x=np.array([[0,0,-1],[0,0,-2]],float);d=np.array([[0,0,1],[.2,0,1]],float);d[1,0]=0 # same hit, distinct speed-normalized input follows
    d[1]=[0,.2,1];x[1]=[0,-.4,-2]
    e=form_arrival_events(x,d,plane())
    assert np.allclose(e.event_geometry["arrival_u"],0);assert np.allclose(e.event_geometry["arrival_v"],0)
    assert not np.allclose(e.event_geometry["arrival_dir_v"][0],e.event_geometry["arrival_dir_v"][1])
    assert np.array_equal(e.receiver_reference["receiver_row_index"],[0,1])


def test_affine_covariance_cross_covariance_and_determinism():
    side=9;y,x=np.mgrid[-4:5,-4:5];launch=np.column_stack((x.ravel(),y.ravel())).astype(float)
    arrival=launch@np.array([[2.,0.],[0.,3.]])
    direction=np.column_stack((.1*arrival[:,0],-.2*arrival[:,1]))
    rel,m=arrival_relations(launch,arrival,direction,side,scales=(1,))
    center=4*side+4
    assert rel["s1_arrival_cov_uv"][center]==pytest.approx(0,abs=1e-14)
    assert rel["s1_arrival_area_ratio"][center]==pytest.approx(6)
    assert rel["s1_cov_u_du"][center]>0 and rel["s1_cov_v_dv"][center]<0
    assert rel["s1_cov_u_dv"][center]==pytest.approx(0,abs=1e-14)
    assert canonical_sha(channel_manifest())==canonical_sha(channel_manifest())


def test_native_plane_contract():
    p=native_receiver_plane();assert p.origin[2]==pytest.approx(4.77);assert p.manifest()["orthogonality_errors"]["eu_ev"]==0
