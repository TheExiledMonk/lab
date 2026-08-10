import numpy as np
from pbuf.wl.geometric_optics import OpticalSurface, make_interaction_records
from pbuf.wl.reverse_transport import reverse_free_propagation, reverse_optical_history


def s(radius=10.): return OpticalSurface("S1","PLANE",np.array([0,0,2.]),np.array([1,0,0.]),np.array([0,1,0.]),np.array([0,0,1.]),
    "CIRCULAR",{"radius":radius},"APERTURE_TEST",{},"EXACTLY_REVERSIBLE")


def test_free_roundtrip():
    p=np.array([[1.,2.,3.]]);d=np.array([[0.,0.,1.]]);t=np.array([7.]);f=p+t[:,None]*d
    assert np.array_equal(reverse_free_propagation(f,d,t),p)


def test_metadata_reverse_and_blocked_loss():
    records,_=make_interaction_records(["e"],np.array([[0.,0.,0.]]),np.array([[0.,0.,1.]]),s(),1)
    c=reverse_optical_history(records)[0]
    assert c.uniqueness_classification=="UNIQUE_INVERSE"
    assert np.array_equal(c.candidate_states[0]["position"],[0,0,0])
    blocked,_=make_interaction_records(["b"],np.array([[2.,0.,0.]]),np.array([[0.,0.,1.]]),s(1.),1)
    c=reverse_optical_history(blocked)[0]
    assert c.uniqueness_classification=="BLOCKED_INFORMATION_NOT_PRESENT_DOWNSTREAM"

