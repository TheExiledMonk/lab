import numpy as np
from pbuf.excitation.native_vector_pair_dynamics import VectorPairState
from pbuf.excitation.native_multi_structure_preparation import NativePreparation, inject, reverse_internal_state

def test_reversal_changes_only_preparation_momentum():
    z=np.zeros((3,3,3,3)); p=np.ones_like(z)
    prep=NativePreparation('A','I','S_PLUS','id','DEV212','fixed','fixed',z,p)
    rev=reverse_internal_state(prep,'S_MINUS')
    assert np.array_equal(prep.displacement_increment,rev.displacement_increment)
    assert np.array_equal(prep.momentum_increment,-rev.momentum_increment)
