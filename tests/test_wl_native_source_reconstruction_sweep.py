from pbuf.wl.native_source_reconstruction_sweep import blind_reconstruct, synthetic_observation, trial_matrix

def test_matrix_and_information_lane_isolation():
    rows=trial_matrix(validation=True)
    assert len(rows)==24
    obs=synthetic_observation(rows[0])
    assert blind_reconstruct(obs,"C1")["information_used"] == ["position"]
    assert blind_reconstruct(obs,"C2")["information_used"] == ["position","direction"]
    assert "trajectory" not in blind_reconstruct(obs,"C3")["information_used"]

def test_trial_order_does_not_change_identity():
    rows=trial_matrix(validation=True)
    assert sorted(x["trial_id"] for x in rows)==sorted(x["trial_id"] for x in reversed(rows))
