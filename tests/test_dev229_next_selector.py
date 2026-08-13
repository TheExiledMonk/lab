from test_dev229_repo_first import load
def test_release_block_selects_dev230_release_gate():
 x=load('dev230_test_selection.json'); assert x['DEV230_TEST_SELECTION']=='NATIVE_SOURCE_RELEASE_REPRESENTATION_GATE' and x['DEV230_TEST_SELECTION_FROZEN']
