import json
from pathlib import Path
RUN=Path(__file__).parents[1]/'runs/dev216_bond_cut_dynamic_polarity'
def test_all_four_frozen_states_are_recorded():
 assert set(json.loads((RUN/'bond_cut_radial_force_matrix.json').read_text())['rows'])=={'++','+-','-+','--'}
