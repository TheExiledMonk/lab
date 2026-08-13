from pbuf.analysis.native_staggered_order import unique_n6_bonds

def test_unique_periodic_n6_bonds_have_three_positive_bonds_per_node():
    pairs, axes, _ = unique_n6_bonds((3, 4, 5))
    assert len(pairs) == 3 * 3 * 4 * 5 and set(axes) == {0, 1, 2}
