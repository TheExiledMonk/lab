from pbuf.wl.native_scale_candidates import NativeScaleEstimate,dependency_graph

def test_shared_trajectory_candidates_not_double_counted():
    a=NativeScaleEstimate("S01","curvature",independence_class="TRAJECTORY_LOCAL")
    b=NativeScaleEstimate("S01b","curvature variant",independence_class="TRAJECTORY_LOCAL")
    g=dependency_graph([a,b])
    assert len(g["edges"])==1
    assert g["effective_independent_support"]==0
