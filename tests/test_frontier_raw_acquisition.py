from __future__ import annotations

import random
import tempfile
import unittest
from pathlib import Path

from astropy.io import fits

from pbuf.data import frontier_raw_acquisition as acq


class FrontierRawAcquisitionTests(unittest.TestCase):
    def test_product_classification_exact(self):
        expected = {
            "jXXXX_raw.fits": acq.RAW_DETECTOR,
            "jXXXX_flt.fits": acq.FLT_CONTROL,
            "jXXXX_flc.fits": acq.FLC_CONTROL,
            "foo_drz.fits": "HST_DRIZZLED",
            "foo_wht.fits": "UNKNOWN_HST_PRODUCT",
        }
        self.assertEqual({name: acq.classify_hst_product(name) for name in expected}, expected)
        self.assertEqual(acq.classify_frontier_file(acq.CONTROL_FILENAME), "HLSP_DRIZZLED_CONTROL")

    def test_drz_never_substitutes_for_missing_raw(self):
        products = acq.build_product_manifest(
            [{"obsid": "1", "productFilename": "foo_drz.fits", "dataURI": "mast:HST/foo"}], {"1"}
        )
        families = acq.exposure_families(products)
        self.assertFalse(products[0]["included"])
        self.assertEqual(products[0]["classification"], "HST_DRIZZLED")
        self.assertEqual(families["foo"]["status"], "RAW_PRODUCT_MISSING")

    def test_field_uses_geometry_not_name(self):
        common = {"obs_collection": "HST", "instrument_name": "ACS/WFC", "filters": "F814W",
                  "intentType": "science", "dataRights": "PUBLIC"}
        main = dict(common, target_name="abell2744-hffpar", s_ra=acq.TARGET_RA, s_dec=acq.TARGET_DEC)
        parallel = dict(common, target_name="abell2744", s_ra=acq.TARGET_RA + 0.13, s_dec=acq.TARGET_DEC)
        self.assertEqual(acq.classify_field(main), "MAIN_CLUSTER")
        self.assertEqual(acq.classify_field(parallel), "PARALLEL_FIELD")
        self.assertTrue(acq.observation_inclusion(main)[0])
        self.assertFalse(acq.observation_inclusion(parallel)[0])

    def test_selection_hash_is_order_independent(self):
        rows = []
        for i, suffix in enumerate(("raw", "flt", "flc")):
            rows.append({"obsid": "10", "rootname": "j123", "filename": f"j123_{suffix}.fits",
                         "product_uri": f"mast:HST/j123_{suffix}.fits", "size": i + 1,
                         "classification": acq.classify_hst_product(f"j123_{suffix}.fits"), "included": True})
        expected = acq.selection_sha256(acq.canonical_selection(rows))
        rng = random.Random(125)
        for _ in range(100):
            rng.shuffle(rows)
            self.assertEqual(acq.selection_sha256(acq.canonical_selection(rows)), expected)

    def test_resume_contract(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "j123_raw.fits"
            path.write_bytes(b"valid")
            digest = acq.sha256_file(path)
            self.assertEqual(acq.resume_action(path, 5, digest), "SKIP_VALID")
            self.assertEqual(acq.resume_action(path, 4, digest), "REDOWNLOAD_CORRUPT")
            self.assertEqual(acq.resume_action(path, 5, "0" * 64), "REDOWNLOAD_CORRUPT")
            path.unlink()
            path.with_name(path.name + ".part").write_bytes(b"partial")
            self.assertEqual(acq.resume_action(path, 20, None), "REDOWNLOAD_PARTIAL")

    def test_synthetic_fits_integrity_and_headers(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "j123_raw.fits"
            primary = fits.PrimaryHDU()
            primary.header.update({
                "TELESCOP": "HST", "INSTRUME": "ACS", "DETECTOR": "WFC",
                "FILTER1": "CLEAR1L", "FILTER2": "F814W", "ROOTNAME": "j123",
                "PROPOSID": 12345, "DATE-OBS": "2014-01-01", "TIME-OBS": "00:00:00",
                "EXPTIME": 100.0, "RA_TARG": acq.TARGET_RA, "DEC_TARG": acq.TARGET_DEC,
                "PA_V3": 42.0,
            })
            fits.HDUList([primary, fits.ImageHDU(data=[[1, 2], [3, 4]])]).writeto(path)
            result = acq.validate_fits(path, acq.RAW_DETECTOR)
            self.assertTrue(result["valid"], result)
            self.assertTrue(result["metadata_valid"])
            self.assertEqual(result["field_classification"], "MAIN_CLUSTER")

    def test_unknown_products_are_retained(self):
        products = acq.build_product_manifest(
            [{"obsid": "1", "productFilename": "j123_x1d.fits", "dataURI": "mast:HST/x"}], {"1"}
        )
        self.assertEqual(len(products), 1)
        self.assertEqual(products[0]["classification"], "UNKNOWN_HST_PRODUCT")
        self.assertFalse(products[0]["included"])


if __name__ == "__main__":
    unittest.main()
