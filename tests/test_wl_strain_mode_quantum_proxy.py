import numpy as np, pytest
from pbuf.wl.quantum_zero_mass_bridge import strain_mode_ratio_bridge

def test_proxy_gate_and_ratios():
    with pytest.raises(ValueError): strain_mode_ratio_bridge([1,.5],proxy_established=False)
    r=strain_mode_ratio_bridge(np.array([1,.5]),proxy_established=True)
    assert np.allclose(r["momentum_ratio"],r["k_ratio"])
    assert np.allclose(r["wavelength_ratio"],[1,2])

