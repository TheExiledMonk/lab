from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pbuf.cosmology.native_background_stress import homogeneous_potential, homogeneous_stress, stress_derivative


def test_frozen_stress_is_strictly_monotone_and_potential_increases_in_extension():
    assert stress_derivative(0.0) == 1.0
    assert stress_derivative(0.8) > 0
    assert homogeneous_stress(1.5) > homogeneous_stress(1.1) > 0
    assert homogeneous_potential(1.5) > homogeneous_potential(1.1) > 0
