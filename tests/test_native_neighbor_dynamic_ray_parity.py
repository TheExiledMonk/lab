from pbuf.foundation.native_neighbor_loaded_excitation import run_case
def test_path_is_native_and_available_without_ray_solver(): assert isinstance(run_case(8,5)['centroid'],float)
