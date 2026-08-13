from test_dev229_repo_first import load
def test_provenance_is_not_identity(): assert load('source_identity_representation.json')['PREPARATION_PROVENANCE_ALONE_NOT_SOURCE_IDENTITY']
