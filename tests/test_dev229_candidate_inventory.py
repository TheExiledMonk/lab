from test_dev229_repo_first import load
def test_candidate_routes_are_inventory_only():
 x=load('persistent_source_candidate_inventory.json'); assert len(x['candidates'])==5; assert x['NO_SOURCE_COUNT_SWEEP'] and x['NO_RESULT_SELECTED_COMPOSITION']
