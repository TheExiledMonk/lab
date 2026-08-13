from pbuf.observer.native_spatial_winding import square_yz_loop

def test_source_centered_multi_cell_loop_is_closed_and_contractible():
    edges = square_yz_loop((11, 11, 11), (1, 5, 5), 3)
    assert len(edges) == 24
    assert {axis for _, axis, _ in edges} == {1, 2}
    assert all(0 <= node[axis] < 11 for node, axis, _ in edges)
