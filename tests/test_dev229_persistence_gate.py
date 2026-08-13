from test_dev229_repo_first import load
def test_persistence_is_blocked_by_release_not_duration():
 x=load('persistent_native_source_derivation.json'); assert x['PERSISTENT_NATIVE_SOURCE_DERIVATION']=='BLOCKED_SOURCE_RELEASE'
