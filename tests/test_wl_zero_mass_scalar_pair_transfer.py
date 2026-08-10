from pbuf.wl.native_zero_mass_scalar_transport import apply_factors, inverse_pair_factor, ratio_transfer

def test_synthetic_pair_roundtrip_exact():
    assert apply_factors(1,[2,inverse_pair_factor(2)]).q_receive == 1

def test_orientation_inverse():
    assert ratio_transfer(2,5)*ratio_transfer(5,2) == 1
