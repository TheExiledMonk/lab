from test_dev229_repo_first import load
def test_no_native_physics_was_added():
 c=load('final_contract.json'); assert all(c[x] for x in ['NO_NEW_FORCE','NO_NEW_DOF','NO_NEW_SOURCE_LAW','NO_NEW_MAGNETIC_PRIMITIVE'])
