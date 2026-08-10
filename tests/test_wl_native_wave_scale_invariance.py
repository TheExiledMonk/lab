from pbuf.wl.native_propagation_units import *
from pbuf.wl.native_wave_state import redshift_from_wavelength

def test_native_rescaling_invariance():
    z=redshift_from_wavelength(2,3)
    for alpha in (.5,1,2,4):
        assert redshift_from_wavelength(2*alpha,3*alpha) == z
        a=PropagationParameterAudit("tau","NATIVE_TIME_EXPLICIT","clock","tau+=dt","dx/dtau","d/dtau",True,True)
        assert l0_over_t0(2*alpha,a)["value"] * alpha == C_SI/2
