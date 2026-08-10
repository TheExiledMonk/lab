import numpy as np
from pbuf.excitation.native_excitation_n6 import N6_OFFSETS,NativeExcitationN6State,gaussian_packet,propagate_directional,quadratic_norm
def test_each_n6_permutation_is_reversible_and_norm_preserving():
    x=gaussian_packet((12,11,10),center=(4,5,5))
    for d in N6_OFFSETS:
      s=NativeExcitationN6State(x.copy()); propagate_directional(s,3,d); assert np.isclose(quadratic_norm(s.values),quadratic_norm(x),rtol=0,atol=1e-12)
      propagate_directional(s,3,tuple(-v for v in d)); assert np.array_equal(s.values,x)
