import pytest
from pbuf.wl.native_propagation_units import *

def test_false_time_and_solver_iteration_rejected():
    a=current_propagation_parameter_audit()
    assert l0_over_t0(2,a)["outcome"] == "NATIVE_TIME_NOT_ESTABLISHED"
    with pytest.raises(ValueError): reject_solver_iteration_as_time("solver_iteration")

def test_known_speed_and_messenger_independence():
    a=PropagationParameterAudit("tau","NATIVE_TIME_EXPLICIT","native clock","tau+=dt","dx/dtau","direction derivative",True,True)
    assert l0_over_t0(2,a)["value"] == C_SI/2
    assert compare_messenger_speeds(1,1,shared_operator=True)["classification"] == "ALGEBRAICALLY_IDENTICAL_SPEED_CONSTRAINT"
    assert compare_messenger_speeds(1,1)["classification"] == "INDEPENDENT_COMMON_SPEED_CONSTRAINT"
