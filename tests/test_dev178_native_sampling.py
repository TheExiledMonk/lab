import numpy as np
from tools.generate_dev178_high_density_native_vulkan import quarter_packet
from tools.pbuf_native_receipt_viewer import load_receipt

def test_atomic_quarter_lane_is_deterministic_and_never_adaptive():
    packet = np.ones((7, 7))
    first, mask = quarter_packet(packet); second, again = quarter_packet(packet)
    assert np.array_equal(first, second) and np.array_equal(mask, again)
    assert mask.sum() == 12 and np.count_nonzero(first) == 12

def test_viewer_loader_preserves_arrays(tmp_path):
    path = tmp_path / "receipt.npz"; source = np.array([[1., np.nan, 3.]])
    np.savez_compressed(path, source_positions=source)
    got = load_receipt(path)
    assert np.isnan(got["source_positions"][0, 1])
    got["source_positions"][0, 0] = 99
    assert np.load(path, allow_pickle=False)["source_positions"][0, 0] == 1
