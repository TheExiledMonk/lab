import numpy as np
from pbuf.quantum.native_localized_excitation_states import construct_composite,state_registry,state_observables
def test_composite_is_diagnostic_not_bound():
    s=construct_composite(np.ones(16),np.ones((16,2)))
    assert len(state_registry())==20 and state_observables(s)['classification']=='SYNTHETIC_OVERLAP_NOT_BOUND'
