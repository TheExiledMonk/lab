import pytest
from pbuf.matter.native_loaded_propagation import LoadedPropagationState

def test_state_rejects_super_c_and_negative_beta():
    for beta in (-.01,1.01):
        with pytest.raises(ValueError): LoadedPropagationState(.5,1,(1,0,0),beta)

