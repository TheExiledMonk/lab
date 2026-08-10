from pbuf.audit.native_dependency_graph import graph
def test_no_inferred_cross_edge():
    g=graph(); assert all(e["kind"]=="CODE_DEPENDENCY" for e in g["edges"])
    assert not any(e["from"]=="u accumulated" and e["to"]=="excitation progression" for e in g["edges"])
