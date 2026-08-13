from test_dev229_repo_first import load
def test_no_result_selected_source_parameters():
 c=load('final_contract.json'); assert c['NO_RESULT_SELECTED_DURATION'] and c['NO_RESULT_SELECTED_SOURCE_SHAPE'] and c['NO_RESULT_SELECTED_SOURCE_COUNT']
