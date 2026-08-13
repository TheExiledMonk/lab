import json
from pathlib import Path
RUN=Path(__file__).parents[1]/'runs/dev216_bond_cut_dynamic_polarity'
def test_pair_gate_blocks_classification_when_not_equal_and_opposite():
 assert json.loads((RUN/'conservation_clean_dynamic_force_sign.json').read_text())['CONSERVATION_CLEAN_DYNAMIC_FORCE_SIGN']=='UNRESOLVED'
