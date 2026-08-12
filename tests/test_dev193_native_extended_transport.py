import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pbuf.observer.native_extended_transport import kernel_moment_decomposition

OUT = ROOT / "runs/dev193_native_extended_transport"

def load(name): return json.loads((OUT / name).read_text())

def test_contract_records_the_physical_boundary_and_closed_gates():
    c = load("final_contract.json")
    assert c["NATIVE_EXTENDED_TRANSPORT"] == "BLOCKED_INITIAL_STATE_COMPOSITION"
    assert c["PHYSICAL_EXTENDED_SOURCE_SUPERPOSITION"] == "BLOCKED_NONLINEAR_INTERACTION"
    assert c["SPIN2_OBSERVABLE_GATE"] == "CLOSED"
    assert c["NO_ASTRONOMICAL_SOURCE_IMAGE"]

def test_exact_finite_moment_decomposition():
    z = np.load(ROOT / "runs/dev188_native_source_distribution_pushforward/transfer_kernel_weight.npz", allow_pickle=False)
    s = np.ones(121)
    d = kernel_moment_decomposition(z["R00_weight"], z["R00_coordinates"], s)
    assert np.allclose(d["output_covariance"], d["centroid_covariance"] + d["response_covariance"], atol=1e-12)

def test_required_controls_pass():
    for name in ("delta_control.json", "two_source_linearity_control.json", "affine_kernel_control.json",
                 "response_broadening_control.json", "throughput_gradient_control.json", "nonaffine_control.json"):
        assert load(name)["pass"]
