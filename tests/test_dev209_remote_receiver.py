import numpy as np

from pbuf.observer.native_remote_receiver import receiver_planes


def test_receivers_are_fixed_reflection_paired_planes():
    yz = np.ones((11, 11), dtype=bool)
    receivers = receiver_planes((11, 11, 11), yz)
    assert [receivers[x]["same_plane_x"] for x in ("R1", "R2", "R3")] == [4, 5, 6]
    assert [receivers[x]["reversed_plane_x"] for x in ("R1", "R2", "R3")] == [6, 5, 4]
    assert all(not np.any(v["SAME"] & v["REVERSED"]) or k == "R2" for k, v in receivers.items())
