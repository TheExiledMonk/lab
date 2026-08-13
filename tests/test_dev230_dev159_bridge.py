from test_dev230_repo_first import load
def test_dev159_bridge_is_narrow():
    assert load('dynamic_source_to_propagating_disturbance.json')['DYNAMIC_SOURCE_TO_PROPAGATING_DISTURBANCE']=='PARTIAL'
