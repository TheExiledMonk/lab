from pbuf.foundation.native_neighbor_loaded_excitation import run_matrix
def test_full_matrix_and_no_cross_coefficient():
    rows=run_matrix(32); assert len(rows)==72 and all(r['new_interaction_coefficients']==0 for r in rows)
