import numpy as np
from pathlib import Path
RUN=Path(__file__).parents[1]/'runs/dev216_bond_cut_dynamic_polarity'
def test_trajectory_contains_four_states_and_frozen_window():
 z=np.load(RUN/'bond_cut_force_trajectory.npz'); assert z['radial_A'].shape==(4,13)
