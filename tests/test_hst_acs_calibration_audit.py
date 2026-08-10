from __future__ import annotations

import unittest
import numpy as np

from pbuf.data import hst_acs_calibration_audit as audit


class HstAcsCalibrationAuditTests(unittest.TestCase):
    def setUp(self):
        self.geometry = {"shape": [16, 16], "CCDCHIP": 1, "LTV1": 0, "LTV2": 0,
                         "LTM1_1": 1, "LTM1_2": 0, "LTM2_1": 0, "LTM2_2": 1,
                         "BINAXIS1": 1, "BINAXIS2": 1}

    def test_identity(self):
        a = np.arange(256.0).reshape(16, 16)
        result = audit.difference_stats(a, a)
        self.assertEqual(result["fraction_exactly_unchanged"], 1)
        self.assertEqual(result["rms_difference"], 0)
        self.assertEqual(audit.geometry_classification(self.geometry, self.geometry), "PIXEL_PRESERVING")

    def test_uniform_additive_bias(self):
        raw = np.arange(256.0).reshape(16, 16); delta = (raw - 7) - raw
        self.assertEqual(audit.classify_difference(delta)["pattern"], "UNIFORM")
        self.assertEqual(audit.difference_stats(raw, raw-7)["mean_difference"], -7)

    def test_flat_field_is_pixel_preserving_multiplication(self):
        raw = np.arange(256.0).reshape(16, 16)+1; flat = np.linspace(.8,1.2,256).reshape(16,16)
        calibrated = raw/flat
        self.assertEqual(audit.geometry_classification(self.geometry, self.geometry), "PIXEL_PRESERVING")
        np.testing.assert_allclose(calibrated*flat, raw)

    def test_resampling(self):
        source = np.arange(256.0).reshape(16,16)
        translated = .75*source + .25*np.roll(source, 1, axis=1)
        self.assertEqual(audit.detect_resampling(source, translated), "RESAMPLED")

    def test_pixel_replacement(self):
        a=np.zeros((16,16)); b=a.copy(); b[3:6,4:8]=9
        result=audit.classify_difference(b-a)
        self.assertEqual(result["pattern"], "LOCALIZED")
        self.assertEqual(result["information"], "INFORMATION_MODIFYING")

    def test_dq_exact(self):
        a=np.zeros((4,4),dtype=np.uint16); b=a.copy(); b[1,2]=4; b[3,3]=8
        result=audit.dq_comparison(a,b)
        self.assertEqual(result["changed_pixel_count"],2)
        self.assertEqual(result["new_value_histogram"],{"0":14,"4":1,"8":1})

    def test_cte_like_vertical_trail(self):
        a=np.zeros((32,32)); a[8,10]=100
        b=a.copy(); b[8,10]-=15; b[9,10]+=10; b[10,10]+=5
        self.assertAlmostEqual(a.sum(),b.sum())
        self.assertEqual(audit.classify_difference(b-a)["information"],"INFORMATION_MODIFYING")
        result=audit.classify_difference(b-a)
        self.assertGreater(result["row_profile_activity"],0)

    def test_unit_guard(self):
        self.assertEqual(audit.direct_difference_status(self.geometry,self.geometry,"COUNTS","ELECTRONS"),"DIRECT_DIFFERENCE_NOT_VALID")

    def test_deterministic_hash(self):
        self.assertEqual(audit.canonical_sha256({"b":2,"a":1}),audit.canonical_sha256({"a":1,"b":2}))

if __name__ == "__main__": unittest.main()
