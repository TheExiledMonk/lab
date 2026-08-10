import pytest
from pbuf.wl.medium_dimensional_closure import emergent_time_mapping,reject_fundamental_time_dimension

def test_fundamental_time_and_iteration_rejected():
    with pytest.raises(ValueError,match="REJECT_FUNDAMENTAL_TIME_DIMENSION"): reject_fundamental_time_dimension(T0=1)
    with pytest.raises(ValueError,match="REJECT_FUNDAMENTAL_TIME_DIMENSION"): reject_fundamental_time_dimension(solver_iterations=10)

def test_emergent_time_only_from_length():
    assert emergent_time_mapping(2,3,c=6)["elapsed_seconds"]==1
    assert emergent_time_mapping(None,3)["available"] is False
