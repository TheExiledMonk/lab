import numpy as np
from pbuf.quantum.native_excitation_modes import carrier_mode,propagate
def test_all_modes_advance_by_exact_permutation():
    for w in (4,8,16,32):
        q=carrier_mode(128,w); assert np.allclose(propagate(q,7)[-1],np.roll(q,7,axis=0))

