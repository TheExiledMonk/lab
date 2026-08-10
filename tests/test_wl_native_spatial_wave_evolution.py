import numpy as np
import pytest
from pbuf.wl.native_spatial_wave_evolution import candidate_registry,accumulate_log_wavelength,evolve_wavelength

def test_all_candidates_audited_without_manufacturing_state():
    x=candidate_registry(); assert len(x)==12 and all(r["status"]=="MISSING_NATIVE_STATE" for r in x)
def test_forward_reverse_and_nonreversible_control():
    s=np.linspace(0,2,101); q=np.sin(s); f=evolve_wavelength(s,q); r=evolve_wavelength(s,q,f[-1],orientation="reverse")
    assert abs(r[-1]-1)<1e-12
    with pytest.raises(ValueError,match="NON_REVERSIBLE"): evolve_wavelength(s,q,orientation="reverse",reversible=False)
