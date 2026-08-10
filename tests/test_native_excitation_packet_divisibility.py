import numpy as np
from pbuf.quantum.native_excitation_modes import carrier_mode
from pbuf.quantum.native_excitation_quantization import divisibility_audit
def test_arbitrary_requested_norm_fractions_are_preserved():
    r=divisibility_audit(carrier_mode(128,16)); assert all(np.isclose(x['requested_norm_fraction'],x['measured_norm_fraction']) for x in r['fractions'])

