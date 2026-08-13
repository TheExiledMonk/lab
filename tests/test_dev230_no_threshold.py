from test_dev230_repo_first import load
def test_no_threshold_or_fit():
    c=load('final_contract.json'); assert c['NO_EMISSION_THRESHOLD'] and c['NO_AMPLITUDE_FIT'] and c['NO_FREQUENCY_FIT']
