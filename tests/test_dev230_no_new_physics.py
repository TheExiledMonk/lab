from test_dev230_repo_first import load
def test_no_new_physics():
    c=load('final_contract.json'); assert all(c[k] for k in ('NO_NEW_PHYSICS','NO_NEW_FORCE','NO_NEW_DOF','NO_NEW_SOURCE_LAW','NO_NEW_EM_FIELD','NO_NEW_MAGNETIC_PRIMITIVE'))
