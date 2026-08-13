import json
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pbuf.observer.native_transfer_distortion import affine_centroids, periodic_centroid_jacobian, polar_fields

OUT = ROOT / "runs/dev190_source_independent_transfer_distortion"


def test_affine_recovery_away_from_periodic_seam():
    a = np.array([[1.7, .3], [.2, .8]])
    j = periodic_centroid_jacobian(affine_centroids(a)).reshape(11, 11, 2, 2)
    assert np.allclose(j[1:-1, 1:-1], a)
    _, _, eigen, aniso, det, _ = polar_fields(j.reshape(-1, 2, 2))
    assert np.allclose(det.reshape(11, 11)[1:-1, 1:-1], np.linalg.det(a))
    assert np.isfinite(eigen).all() and np.isfinite(aniso).all()


def test_dev190_contract_and_rerun_are_deterministic():
    subprocess.check_call([sys.executable, "tools/generate_dev190_source_independent_transfer_distortion.py"], cwd=ROOT)
    c = json.loads((OUT / "final_contract.json").read_text())
    for key in ("SCALAR_OPERATOR_RANK_33_VERIFIED", "MULTICHANNEL_RANK_121_VERIFIED",
                "SYNTHETIC_AFFINE_CONTROLS_PASS", "PIPELINE_DETERMINISTIC",
                "NO_ARBITRARY_CHANNEL_NORMALIZATION"):
        assert c[key]
    assert c["SOURCE_INDEPENDENT_TRANSFER_DISTORTION"] == "DERIVED_LOCAL_TENSOR"
