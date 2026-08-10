import numpy as np
import pytest
from pbuf.wl.reverse_transport import (ReverseCandidateSet, affine_inverse, arrival_knn,
    correspondence_index, reconstruct_receiver, roundtrip_errors, transport_diagnostics)

@pytest.mark.parametrize("a,b",[(np.eye(2),[0,0]),(np.eye(2),[2,-3]),
    (np.array([[0,-1],[1,0]]),[0,0]),(np.array([[2,.3],[-.1,3]]),[1,2]),
    (np.diag([1e4,1.]),[0,0])])
def test_analytic_roundtrips(a,b):
    x=np.array([[-2.,1.],[0,0],[3,4]]);y=x@a.T+np.asarray(b);got,status=affine_inverse(y,a,b)
    assert status=="UNIQUE_INVERSE";assert np.allclose(got,x,rtol=1e-12,atol=1e-12)

def test_singular_and_many_to_one_preserve_ambiguity():
    got,status=affine_inverse([[1,0]],[[1,0],[0,0]],[0,0]);assert got is None and status=="NON_UNIQUE_INVERSE"
    c=ReverseCandidateSet.from_candidates("event",[4,9]);assert c.uniqueness_classification=="MULTIPLE_LAUNCH_CANDIDATES"

def test_surface_intersection_reversal_and_direction_deletion_control():
    x=np.array([[1.,2.,-3.],[-2.,4.,-1.]])
    d=np.array([[0.,0.,1.],[.2,-.1,1.]]);d/=np.linalg.norm(d,axis=1)[:,None]
    t=np.array([3.,2.]);p=x+t[:,None]*d;got=reconstruct_receiver(p,d,t)
    assert roundtrip_errors(got,x)["classification"] in ("EXACT_ROUNDTRIP","NUMERICALLY_EXACT_ROUNDTRIP")
    # Position+t alone leaves tangential endpoint coordinates unconstrained.
    alt_d=d.copy();alt_d[:,:2]=0;alt_d[:,2]=1
    assert np.linalg.norm(reconstruct_receiver(p,alt_d,t)-x)>np.linalg.norm(got-x)

def test_conditioning_orientation_correspondence_and_neighbors_deterministic():
    a=np.array([np.eye(2),np.diag([1e6,1]),np.array([[-1,0],[0,1]]),np.array([[1,0],[0,0]])])
    q=transport_diagnostics(a);assert q["detJ"][2]<0 and q["classification"][3]=="LOCALLY_SINGULAR_CANDIDATE"
    c=correspondence_index([2,1,2],[4,5,6]);assert c["launch_event_indices"].tolist()==[1,0,2]
    p=np.array([[0,0],[1,0],[-1,0],[0,1],[0,-1]],float)
    assert np.array_equal(arrival_knn(p,ks=(4,))["k4"],arrival_knn(p,ks=(4,))["k4"])

def test_bundle_deletion_control():
    a=np.diag([100.,1.]);single=transport_diagnostics(a[None])["transport_condition_number"][0]
    # An independent bundle constraint on the weak axis makes the stacked normal matrix isotropic.
    stacked=np.vstack((a,np.array([[0,100.]])));s=np.linalg.svd(stacked,compute_uv=False)
    assert s[-1]/s[0] > single
