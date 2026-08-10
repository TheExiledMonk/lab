import numpy as np
import pytest

from pbuf.wl.receiver_sky_bridge import (BRIDGE_VERSION, PBUFReceiverSkyBridge,
                                         canonical_sha256, roundtrip_audit)


def bridge(rotation=((1.0, 0.0), (0.0, 1.0)), reflection=False):
    return PBUFReceiverSkyBridge(
        BRIDGE_VERSION, "synthetic origin", (1, 0, 0), (0, 1, 0), (0, 0, 1),
        "EXPLICIT_ANGULAR_SCALE", 2.0, "arcsec/A0_unit", {"synthetic": True},
        (3.0, -4.0), "synthetic tangent plane", rotation, reflection,
        {"translation": "UPSTREAM_FIXED", "rotation": "UPSTREAM_FIXED",
         "scale": "UPSTREAM_FIXED", "reflection": "UPSTREAM_FIXED"},
        ("test",), canonical_sha256({"test": True}))


@pytest.mark.parametrize("degrees", [0, 90, 180, 270])
def test_scale_rotation_translation_roundtrip(degrees):
    angle = np.deg2rad(degrees)
    r = ((np.cos(angle), -np.sin(angle)), (np.sin(angle), np.cos(angle)))
    b = bridge(r)
    points = np.array([[0., 0.], [1., 0.], [-2., 5.]])
    assert roundtrip_audit(b, points)["max_error"] < 1e-12


def test_known_scale():
    b = bridge()
    np.testing.assert_allclose(b.forward([[1., 0.]]), [[5., -4.]])


def test_reflection_is_explicit():
    b = bridge(((-1., 0.), (0., 1.)), True)
    assert b.reflection_status
    np.testing.assert_allclose(b.reverse(b.forward([[2., 3.]])), [[2., 3.]])


def test_missing_scale_never_defaults_to_one():
    with pytest.raises(ValueError, match="PBUF_RECEIVER_PHYSICAL_SCALE_NOT_ESTABLISHED"):
        PBUFReceiverSkyBridge(BRIDGE_VERSION, "x", (1,0,0), (0,1,0), (0,0,1),
          "NUMERICAL_ONLY_SCALE", 1., "native", {}, (0,0), "x", ((1,0),(0,1)), False, {}, (), "x")
