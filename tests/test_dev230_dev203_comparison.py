from test_dev230_repo_first import load
def test_no_forced_wave_identity():
    assert load('source_generated_residual_to_dev203_wave_relation.json')['SOURCE_GENERATED_RESIDUAL_TO_DEV203_WAVE_RELATION']=='NOT_COMPARABLE'
