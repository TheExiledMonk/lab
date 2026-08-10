from pbuf.excitation.native_excitation_transfer import dependency_contract

def test_no_ray_or_static_response_dependency():
    c=dependency_contract(); assert not c['trajectory_solver_used_to_move_excitation']
    assert not c['static_response_used_as_dynamic_state'] and c['free_coefficients']==0
