import numpy as np
from pbuf.quantum.native_localized_excitation_states import stability_analysis
def test_free_history_does_not_claim_localized_stability():
    r=stability_analysis(np.zeros((2,3,8,2)))
    assert not r['stable_localized_state'] and r['classification']=='FREE'
