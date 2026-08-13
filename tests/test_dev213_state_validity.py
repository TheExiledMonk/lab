import numpy as np
from pbuf.observer.native_composition_audit import state_validity

def test_unstrained_state_is_valid():
    result=state_validity(np.zeros((3,3,3,3)),np.zeros((3,3,3,3)))
    assert result['classification']=='VALID' and result['max_abs_strain']==0.0
