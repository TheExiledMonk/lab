from pbuf.foundation.native_loaded_link_observer import compare_path
def test_missing_frozen_ray_is_not_claimed_as_parity(): assert compare_path([[0,0]])["status"]=="UNDERDETERMINED"
