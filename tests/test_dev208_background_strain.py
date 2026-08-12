import json
from pathlib import Path
import numpy as np
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pbuf.cosmology.native_background_strain import extension, homogeneous_relations


def test_uniform_n6_scaling_has_six_isotropic_relations():
    relations = homogeneous_relations(1.25)
    assert relations.shape == (6, 3)
    assert np.allclose(np.linalg.norm(relations, axis=1), 1.25)
    assert extension(1.25) == 0.25


def test_dev208_does_not_assume_cosmological_mapping():
    run = Path(__file__).parents[1] / "runs/dev208_native_cosmic_turnaround"
    result = json.loads((run / "scale_factor_native_spacing_mapping.json").read_text())
    assert result["COSMOLOGICAL_SCALE_FACTOR_TO_NATIVE_SPACING"] == "NOT_DERIVED"
    assert result["NO_ASSUMED_SCALE_FACTOR_MAPPING"]
