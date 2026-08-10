import numpy as np
import pytest
from pbuf.wl.optical_interaction_state import availability_masks, state_from_event, validate_derivation_graph


def base():
    return {"event_uid":"x","ray_index":1,"receiver_row_index":2,"launch_grid_index":3,
            "arrival_u":.1,"arrival_v":.2,"arrival_dir_u":0.,"arrival_dir_v":0.,"arrival_dir_n":1.}


def test_optional_physical_state_exact_preservation():
    e=base()|{"physical_weight":2.5,"wavelength":814e-9,"frequency":3.68e14,
              "path_length":4.,"arrival_time":8.,"relative_delay":.2,"phase":.3,
              "stokes":[1.,.1,.2,.3],"final_J11":1.}
    s=state_from_event(e)
    assert s.carried_signal["physical_weight"] == 2.5
    assert s.spectral["wavelength"] == 814e-9
    assert s.temporal == {"path_length":4.,"arrival_time":8.,"relative_delay":.2}
    assert s.phase_coherence["phase"] == .3 and s.polarization["stokes"] == [1.,.1,.2,.3]
    assert all(s.availability_metadata.values())


def test_missing_state_has_masks_and_no_defaults():
    s=state_from_event(base())
    assert not s.availability_metadata["has_physical_weight"]
    assert not s.availability_metadata["has_spectral_state"]
    assert s.carried_signal == {} and s.spectral == {}
    m=availability_masks([base(),base()|{"physical_weight":1.}])
    assert np.array_equal(m["has_physical_weight"],[False,True])


def test_cycle_detection():
    validate_derivation_graph({"edges":[{"source":"a","derived":"b"}]})
    with pytest.raises(ValueError,match="OPTICAL_STATE_DERIVATION_CYCLE"):
        validate_derivation_graph({"edges":[{"source":"a","derived":"b"},{"source":"b","derived":"a"}]})
