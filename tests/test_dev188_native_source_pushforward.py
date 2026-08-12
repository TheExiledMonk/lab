import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pbuf.observer.native_transfer import NativeIncidentDistribution, NativeTransferOperator

OUT = ROOT / "runs/dev188_native_source_distribution_pushforward"


def test_dev188_transfer_contract_and_kernel_shape():
    contract = json.loads((OUT / "final_contract.json").read_text())
    assert contract["NATIVE_SOURCE_PUSHFORWARD"] == "DERIVED"
    assert contract["OBSERVATIONAL_COMPARISON_GATE"] == "CLOSED"
    assert contract["ALL_EIGHT_TRANSFER_KERNELS_BUILT"]
    kernel = np.load(OUT / "transfer_kernel_weight.npz", allow_pickle=False)
    for realization in range(8):
        weight = kernel[f"R{realization:02d}_weight"]
        assert weight.shape == (33, 121)
        assert np.all(weight >= 0)


def test_native_transfer_is_explicit_id_keyed_linear_pushforward():
    kernel = np.load(OUT / "transfer_kernel_weight.npz", allow_pickle=False)
    source_manifest = json.loads((OUT / "source_launch_domain.json").read_text())
    ids = tuple(row["launch_id"] for row in source_manifest["states"])
    operator = NativeTransferOperator("R00", ids, kernel["R00_cell_ids"],
                                      kernel["R00_coordinates"], kernel["R00_weight"], {})
    source = NativeIncidentDistribution(ids, np.eye(1, 121, 60)[0])
    assert np.array_equal(operator.pushforward(source), kernel["R00_weight"][:, 60])
    with np.testing.assert_raises(ValueError):
        operator.pushforward(NativeIncidentDistribution(tuple(reversed(ids)), source.values))
