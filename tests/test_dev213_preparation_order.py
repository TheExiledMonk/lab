import numpy as np
from pbuf.excitation.native_vector_pair_dynamics import VectorPairState
from pbuf.excitation.native_multi_structure_preparation import NativePreparation, inject

def test_same_step_preparations_commute():
    z=np.zeros((3,3,3,3)); base=VectorPairState(z,z)
    a=NativePreparation('A','I','+','id','DEV196','fixed','fixed',np.ones_like(z),np.ones_like(z))
    b=NativePreparation('B','I','+','id','DEV196','fixed','fixed',2*np.ones_like(z),-np.ones_like(z))
    ab=inject(inject(base,a),b); ba=inject(inject(base,b),a)
    assert np.array_equal(ab.displacement,ba.displacement) and np.array_equal(ab.momentum,ba.momentum)
