import numpy as np

from pbuf.observer.native_orientation_transport import orientation_packets, reflected_x


def test_dev207_reflection_is_reused_exactly():
    packet = np.arange(4 * 3 * 3 * 3.0).reshape(4, 3, 3, 3)
    assert np.array_equal(reflected_x(reflected_x(packet)), packet)
    states = orientation_packets(packet, packet)
    assert np.array_equal(states["SAME"][0], packet)
    assert np.array_equal(states["REVERSED"][0], reflected_x(packet))
