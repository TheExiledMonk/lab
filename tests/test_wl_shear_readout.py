import unittest

import numpy as np

from pbuf.wl.shear_readout import (
    build_shear_candidates, candidate_bank_sha256, rotate_pair_components,
    spin2_transform, synthetic_gate_report, traceless_pair,
)


class ShearReadoutTests(unittest.TestCase):
    def test_bank_is_deterministic_bounded_and_paired(self):
        a, b = build_shear_candidates(), build_shear_candidates()
        self.assertEqual(a, b)
        self.assertEqual(len(a), 45)
        self.assertLessEqual(len(a), 64)
        self.assertEqual(candidate_bank_sha256(a), candidate_bank_sha256(b))
        self.assertTrue(any(x.requires_3d for x in a))
        self.assertTrue(any(x.requires_kde for x in a))
        self.assertTrue(all(x.primitive1.endswith("q1") for x in a))

    def test_spin_two_rotation_and_reflection(self):
        x = np.array([-.8, -.1, .4, 1.2])
        y = np.array([.3, -1., .7, .2])
        q = traceless_pair(x, y)
        for angle in (0, 45, 90, 135):
            theta = np.deg2rad(angle)
            actual = traceless_pair(*rotate_pair_components(x, y, theta))
            expected = spin2_transform(*q, theta)
            np.testing.assert_allclose(actual, expected, atol=1e-12)
        reflected = traceless_pair(-x, y)
        np.testing.assert_allclose(reflected[0], q[0])
        np.testing.assert_allclose(reflected[1], -q[1])

    def test_all_structural_gates_pass(self):
        report = synthetic_gate_report()
        for name, passed in report.items():
            if name != "rotation_max_abs_error":
                self.assertTrue(passed, name)


if __name__ == "__main__":
    unittest.main()
