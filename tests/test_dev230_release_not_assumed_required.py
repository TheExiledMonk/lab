from test_dev230_repo_first import load
def test_release_is_not_an_em_prerequisite():
    d=load('source_release_dependency_for_em_generation.json'); assert d['SOURCE_RELEASE_DEPENDENCY_FOR_EM_GENERATION']=='NOT_REQUIRED'; assert d['NO_RELEASE_EVENT_REQUIRED_BY_ASSUMPTION']
