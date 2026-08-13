from test_dev229_repo_first import load
def test_release_block_is_preserved(): assert load('native_source_release_semantics.json')['NATIVE_SOURCE_RELEASE_SEMANTICS']=='NOT_DERIVED'
