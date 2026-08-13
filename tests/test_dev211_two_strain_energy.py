import numpy as np

from pbuf.observer.native_interaction_energy import interaction_residual


def test_dev211_four_state_inclusion_exclusion():
    assert np.array_equal(interaction_residual(np.array([10]), np.array([3]), np.array([5]), np.array([1])), np.array([3]))
