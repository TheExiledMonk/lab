import numpy as np
from pbuf.wl.instrument_interaction import (CONTRACT_VERSION, default_contract, event_uid,
    event_uids, input_availability)

def test_event_uid_determinism_and_identity_sensitivity():
    assert event_uid(1,2,3)==event_uid(1,2,3)
    assert len({event_uid(1,2,3),event_uid(1,2,4),event_uid(2,2,3)})==3
    u=event_uids([1,2],[3,4]);assert len(np.unique(u))==2

def test_contract_canonical_determinism_and_no_optical_invention():
    availability=input_availability({"arrival_u","arrival_v","arrival_dir_u","arrival_dir_v","arrival_dir_n",
        "receiver_incidence_cosine","receiver_incidence_angle","path_length","path_excess",
        "final_J11","final_J12","final_J21","final_J22","ray_index","receiver_row_index","launch_grid_index"})
    c=default_contract(availability);assert c.version==CONTRACT_VERSION and c.sha256==default_contract(availability).sha256
    a={x["category"]:x["availability"] for x in availability}
    assert a["I4"]==a["I6"]==a["I7"]=="NOT_DEFINED_IN_CURRENT_MODEL"
    assert not c.target_access and not c.hst_pixel_access
