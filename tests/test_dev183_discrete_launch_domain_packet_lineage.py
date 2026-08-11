import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'runs/dev183_discrete_launch_domain_packet_lineage'
def test_dev183_domain_and_lineage_closure():
 final=json.loads((OUT/'final_contract.json').read_text()); domain=json.loads((OUT/'domain_cardinality.json').read_text())
 assert final['OUTCOME']=='OUTCOME_A' and final['DENSITY_CONVERGENCE_AUTHORIZED']
 assert domain['N_admissible']==121 and domain['launch_domain_dimension']=='2D'
 assert json.loads((OUT/'serialization_roundtrip.json').read_text())['SERIALIZATION_ROUNDTRIP_PASS']
