import numpy as np
from pbuf.quantum.native_excitation_modes import carrier_mode
from pbuf.quantum.native_excitation_interference import interference_audit
def test_constructive_and_destructive_interference():
    q=carrier_mode(128,16); assert np.isclose(interference_audit(q,q)['norm_superposition'],4*interference_audit(q,q)['norm_a']); assert np.isclose(interference_audit(q,-q)['norm_superposition'],0)

