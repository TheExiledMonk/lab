import numpy as np
import pytest
from pbuf.wl.native_zero_mass_scalar import ZeroMassScalarSource, ZeroMassScalarState, source_scalar_ontology_contract

def test_source_scalar_is_neutral_initial_condition():
    source=ZeroMassScalarSource(4.0); state=ZeroMassScalarState(source.q_scalar)
    assert state.q_ratio == 1.0
    c=source_scalar_ontology_contract()
    assert c["source_scalar_generated_by_medium"] is False
    assert c["identified_physical_quantity"] is None

def test_relative_state_rejects_zero_emit():
    with pytest.raises(ZeroDivisionError): _=ZeroMassScalarState(0.0).q_ratio
