import pytest
from pbuf.wl.native_energy_redshift import redshift_from_mode_energy_proxy

def test_no_redshift_without_proxy():
    with pytest.raises(ValueError): redshift_from_mode_energy_proxy([1,.9],proxy_established=False)
    assert redshift_from_mode_energy_proxy([1,.5],proxy_established=True)[1]==1
