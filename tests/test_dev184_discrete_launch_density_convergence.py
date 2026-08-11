import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'runs/dev184_discrete_launch_density_convergence'
def test_dev184_complete_finite_launch_ladder():
 f=json.loads((OUT/'final_contract.json').read_text())
 assert f['DEV184_COMPLETE'] and f['ALL_EIGHT_REALIZATIONS_INCLUDED']
 assert f['EXACT_RESET_REPLAY'] and f['PACKET_AWARE_LINEAGE_USED']
 assert json.loads((OUT/'dev183_subset_hash_verification.json').read_text())['counts']=={'BASELINE':1,'C25_DISCRETE':30,'C50_DISCRETE':60,'C100_DISCRETE':121}
