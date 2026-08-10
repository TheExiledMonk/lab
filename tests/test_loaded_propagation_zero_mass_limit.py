import pytest
from pbuf.matter.native_loaded_propagation import LoadedPropagationState, diagnostic_surface

def test_zero_loading_endpoint_is_beta_one_for_every_q():
    _,_,b=diagnostic_surface(); assert (b[0]==1).all()
    LoadedPropagationState(0,4,(1,0,0),1,zero_mass_limit=True)
    with pytest.raises(ValueError): LoadedPropagationState(0,4,(1,0,0),.9,zero_mass_limit=True)

