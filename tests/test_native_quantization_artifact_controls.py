from pbuf.quantum.native_transition_quantization import cluster_final_states,artifact_controls
def test_continuous_family_not_quantized():
    r=cluster_final_states([1,2,3],[1,2,3]); assert r['classification']=='CONTINUOUS_FINAL_STATE_FAMILY' and not r['discrete_attractor']
def test_nonconverged_grid_is_artifact():
    assert artifact_controls([1,2,3],[1,3,6])['classification']=='NUMERICAL_QUANTIZATION_ONLY'
