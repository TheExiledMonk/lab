from pathlib import Path
import numpy as np

RUN = Path(__file__).parents[1] / 'runs/dev218_exact_interface_dynamic_polarity'


def test_force_trajectory_covers_frozen_dev214_window():
    data = np.load(RUN / 'interface_force_trajectory.npz')
    assert set(data['labels']) == {'pp', 'pm', 'mp', 'mm'}
    assert data['radial_force'].shape == (4, 13)
