import numpy as np
import pytest

from pbuf.wl.geometric_optics import (OpticalSurface, aperture_test,
    ideal_direction_transform, make_interaction_records, optical_record_uid,
    propagate_to_surface, system_sha256)


def surface(z=10., **kw):
    args=dict(surface_id="S",surface_type="PLANE",origin=np.array([0.,0.,z]),
              basis_u=np.array([1.,0.,0.]),basis_v=np.array([0.,1.,0.]),normal=np.array([0.,0.,1.]))
    args.update(kw);return OpticalSurface(**args)


def test_free_and_oblique_propagation_reuses_plane_semantics():
    p=np.array([[0.,0.,0.],[1.,2.,0.]])
    d=np.array([[0.,0.,1.],[.6,0.,.8]])
    r=propagate_to_surface(p,d,surface())
    assert np.allclose(r["intersection_position"],[[0,0,10],[8.5,2,10]])
    assert np.allclose(r["outgoing_direction"],d)
    assert np.all(r["intersection_status"]=="FORWARD_INTERSECTION")


def test_surface_rejects_silent_basis_repair():
    with pytest.raises(ValueError,match="BASIS_INVALID"):
        surface(basis_v=np.array([1.,1.,0.]))


def test_aperture_pass_block_and_included_boundary():
    assert list(aperture_test([0.,2.,1.],[0.,0.,0.],1.)) == [
        "PASS_APERTURE","BLOCKED_BY_APERTURE","ON_APERTURE_BOUNDARY"]
    s=surface(aperture_type="CIRCULAR",aperture_parameters={"radius":1.},interaction_type="APERTURE_TEST",
              reverse_classification="EXACTLY_REVERSIBLE")
    records,r=make_interaction_records(["a","b"],np.array([[0,0,0],[2,0,0]],float),np.tile([0,0,1.],(2,1)),s,1)
    assert r["transmitted"].tolist()==[True,False]
    assert records[1].reverse_classification=="INFORMATION_LOSSY"
    assert records[1].reverse_metadata["upstream_state_retained"]


def test_ideal_axis_offaxis_and_symmetry():
    s=surface(z=0.,interaction_type="IDEAL_DIRECTION_TRANSFORM",
              interaction_parameters={"focal_distance":2.},reverse_classification="REVERSIBLE_WITH_METADATA")
    pts=np.array([[0,0,0],[1,0,0],[-1,0,0],[0,1,0],[0,-1,0]],float)
    d=ideal_direction_transform(pts,s)
    assert np.allclose(d[0],[0,0,1])
    assert np.allclose(d[1],[-d[2,0],0,d[2,2]])
    assert np.allclose(d[3],[0,-d[4,1],d[4,2]])


def test_uids_and_system_sha_are_deterministic():
    assert optical_record_uid("e","S1",1)==optical_record_uid("e","S1",1)
    assert optical_record_uid("e","S1",1)!=optical_record_uid("e","S1",2)
    assert system_sha256([surface()])==system_sha256([surface()])

