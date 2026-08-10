from pbuf.matter.native_loaded_propagation import LoadedPropagationState, loading_only_audit

def test_same_loading_can_represent_different_beta_but_loading_only_law_fails():
    a=LoadedPropagationState(.5,1,(1,0,0),.2); b=LoadedPropagationState(.5,4,(1,0,0),.8)
    assert a.rest_loading_state==b.rest_loading_state and a.propagation_fraction_beta!=b.propagation_fraction_beta
    assert loading_only_audit()['status']=='FAILS_VARIABLE_SPEED_SAME_MASS'

