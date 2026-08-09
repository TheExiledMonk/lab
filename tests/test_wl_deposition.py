import unittest

import numpy as np

from pbuf.wl.deposition import METHODS, METHOD_BY_NAME


class DepositionTests(unittest.TestCase):
    def test_inventory_exact(self):
        self.assertEqual([m.name for m in METHODS], [
            "hard_bin_current", "hard_bin_half_open", "nearest_center",
            "bilinear_cic", "tsc_3x3", "gaussian_sigma_half_cell",
        ])

    def test_one_and_many_ray_conservation(self):
        u = np.array([-1.0, -0.75, 0.0, 0.25, 1.0], dtype=np.float64)
        v = np.array([-1.0, 0.75, 0.0, -0.25, 1.0], dtype=np.float64)
        for method in METHODS:
            with self.subTest(method=method.name):
                self.assertAlmostEqual(method.deposit(u[:1], v[:1], None,
                                                       bins=4, extent=1).sum(), 1.0)
                self.assertTrue(np.isclose(method.deposit(u, v, None, bins=4, extent=1).sum(),
                                           u.size, rtol=1e-12, atol=1e-12))

    def test_boundary_cases_conserve(self):
        # center, vertical edge, horizontal edge, four-cell corner, min, max
        u = np.array([-0.75, 0.0, -0.25, 0.0, -1.0, 1.0])
        v = np.array([-0.75, -0.25, 0.0, 0.0, -1.0, 1.0])
        for method in METHODS:
            with self.subTest(method=method.name):
                image = method.deposit(u, v, None, bins=4, extent=1.0)
                self.assertEqual(image.dtype, np.float64)
                self.assertTrue(np.isclose(image.sum(), 6.0, rtol=1e-12, atol=1e-12))

    def test_weighted_symmetry_and_equal_split(self):
        for name in ("bilinear_cic", "tsc_3x3", "gaussian_sigma_half_cell"):
            method = METHOD_BY_NAME[name]
            centered = method.deposit(np.array([0.0]), np.array([0.0]), None,
                                      bins=5, extent=1.0)
            self.assertTrue(np.allclose(centered, centered[::-1, ::-1]))
            split = method.deposit(np.array([0.2]), np.array([0.0]), None,
                                   bins=5, extent=1.0)
            self.assertTrue(np.allclose(split[:, 2], split[:, 3]))

    def test_nearest_center_tie_goes_lower(self):
        image = METHOD_BY_NAME["nearest_center"].deposit(
            np.array([0.0]), np.array([0.0]), None, bins=4, extent=1.0)
        self.assertEqual(image[1, 1], 1.0)

    def test_signed_values_remain_signed(self):
        for method in METHODS:
            with self.subTest(method=method.name):
                image = method.deposit(np.array([0.0]), np.array([0.0]),
                                       np.array([-2.5]), bins=4, extent=1.0)
                self.assertAlmostEqual(image.sum(), -2.5)
                self.assertLessEqual(float(np.max(image)), 0.0)

    def test_out_of_range_and_non_square_rejection(self):
        for method in METHODS:
            self.assertEqual(method.deposit(np.array([2.0]), np.array([0.0]), None,
                                            bins=4, extent=1.0).sum(), 0.0)
            with self.assertRaises(ValueError):
                method.deposit(np.array([0.0]), np.array([0.0]), None,
                               bins=(3, 4), extent=1.0)


if __name__ == "__main__":
    unittest.main()
