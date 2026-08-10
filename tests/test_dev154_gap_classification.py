from pbuf.audit.native_dependency_graph import graph
def test_missing_relations_are_conceptual_not_code_edges():
    g=graph(); pairs={(e["from"],e["to"]) for e in g["conceptual_only"]}; assert ("u accumulated","excitation progression") in pairs and ("rbar history","direction") in pairs
