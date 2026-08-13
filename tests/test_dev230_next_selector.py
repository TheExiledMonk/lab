from test_dev230_repo_first import load
def test_selector_follows_representation_block():
    d=load('dev231_test_selection.json'); assert d['DEV231_TEST_SELECTION']=='SOURCE_TO_WAVE_REPRESENTATION_BRIDGE' and d['DEV231_TEST_SELECTION_FROZEN']
