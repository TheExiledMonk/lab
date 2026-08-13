from test_dev230_repo_first import load
def test_off_and_dynamic_controls_are_distinct():
    assert load('native_emission_off_state.json')['NATIVE_EMISSION_OFF_STATE']=='DERIVED_ZERO'; assert load('dynamic_excited_emission_state.json')['DYNAMIC_EXCITED_EMISSION_STATE']=='ON_PROPAGATING'
