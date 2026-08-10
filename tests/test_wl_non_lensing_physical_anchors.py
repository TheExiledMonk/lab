import pytest
from pbuf.wl.non_lensing_physical_anchors import FrozenGlobalScale,physical_anchor_manifest

def test_manifest_has_only_calibration_or_validation_roles():
    assert all(x["role"] in {"CALIBRATION","VALIDATION"} for x in physical_anchor_manifest()["records"])

def test_global_scale_immutable_after_freeze():
    s=FrozenGlobalScale(); s.calibrate(1,"LAB")
    with pytest.raises(RuntimeError,match="GLOBAL_SCALE_IMMUTABLE_AFTER_FREEZE"):s.calibrate(2,"EARTH")
