from pbuf.quantum.native_excitation_modes import carrier_mode,propagate,stability_audit
def test_exact_shift_is_stable_native_mode(): assert stability_audit(propagate(carrier_mode(128,16),20),16)['classification']=='STABLE_NATIVE_MODE'

