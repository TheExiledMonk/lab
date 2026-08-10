import numpy as np
from pbuf.excitation.excitation_static_medium_coupling import coupling_audit, expose_static_medium

def test_static_medium_is_read_only_and_coupling_unresolved():
    x=expose_static_medium(np.ones(8)); assert not x.flags.writeable
    assert not coupling_audit(x)['coefficient_free_coupling_available']

