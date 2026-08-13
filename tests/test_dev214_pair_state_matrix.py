import numpy as np
from pbuf.excitation.native_multi_structure_preparation import NativePreparation, reverse_internal_state
def test_dev212_reversal_keeps_geometry_and_changes_only_momentum():
 z=np.zeros((3,3,3,3)); p=np.ones_like(z); a=NativePreparation('A','I','+','fixed','DEV212','fixed','fixed',z,p); b=reverse_internal_state(a,'-')
 assert np.array_equal(a.displacement_increment,b.displacement_increment) and np.array_equal(a.momentum_increment,-b.momentum_increment)
