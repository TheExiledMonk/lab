import numpy as np
import pytest
from pbuf.matter.native_mass_loading_state import loading_inventory, normalization_audit, strain_loading_fraction

def test_all_loading_candidates_and_native_strain_fraction():
    rows=loading_inventory()
    assert [r['id'] for r in rows] == [f'L{i:02d}' for i in range(1,21)]
    np.testing.assert_allclose(strain_loading_fraction([0,.5,.9]),[0,.5,.9])
    assert normalization_audit()['energy']['outcome']=='NO_NATIVE_FINITE_W_NORMALIZATION'

def test_strain_fraction_respects_open_constitutive_bound():
    with pytest.raises(ValueError): strain_loading_fraction([1.0])

