from pbuf.foundation.native_neighbor_loaded_excitation import run_case
def test_current_law_does_not_fake_localization(): assert run_case(2,6)['interaction']=='GEOMETRIC_ONLY'
