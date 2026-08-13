import numpy as np
from pbuf.observer.native_pattern_boundary import N6_ORDER, n6_signature
def test_ordered_signature_is_signed_and_canonical():
    assert N6_ORDER == ('+x','-x','+y','-y','+z','-z')
    assert n6_signature(np.zeros((3,3,3,3))).shape == (3,3,3,6)
