from pbuf.wl.spatial_lineage_audit import SpatialTransform, coordinate_lineage

def edge(a,b,physical=True):
    return SpatialTransform("synthetic",a,"m","identity",b,"m",True,True,physical)

def test_complete_synthetic_lineage():
    names=["source_grid","network_grid","propagation_grid","receiver_coordinates"]
    graph=coordinate_lineage(names,[edge(names[i],names[i+1]) for i in range(3)])
    assert graph["complete"] and graph["broken_spatial_lineage_edges"]==0

def test_broken_lineage_is_not_guessed():
    names=["source_grid","network_grid","propagation_grid","receiver_coordinates"]
    graph=coordinate_lineage(names,[edge("source_grid","network_grid"),edge("propagation_grid","receiver_coordinates")])
    assert not graph["complete"] and graph["broken_spatial_lineage_edges"]==1
