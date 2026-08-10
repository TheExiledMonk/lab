import unittest

from pbuf.wl.channel_compatibility import component_availability, validate_final_audit
from pbuf.wl.final_audit_reader import FinalAuditError, load_final_audit


class ChannelCompatibilityTests(unittest.TestCase):
    def test_load_final_audit_uses_last_structured_result(self):
        with self.subTest("last result wins"):
            import tempfile
            from pathlib import Path
            with tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / "final.log"
                path.write_text('RESULT_JSON\n{"old": true}\nnoise\nRESULT_JSON\n{"final": true}\n', encoding="utf-8")
                self.assertEqual(load_final_audit(path), {"final": True})


    def test_load_final_audit_rejects_missing_structured_result(self):
        import tempfile
        from pathlib import Path
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "final.log"
            path.write_text("human readable only\n", encoding="utf-8")
            with self.assertRaises(FinalAuditError):
                load_final_audit(path)


    def test_component_availability_is_not_inferred(self):
        diagnostics = {"candidate": {"kappa": {"pearson": 1.0}}}
        audit = {"clusters": {cluster: {"methods": {method: {
            "observational_diagnostics": diagnostics} for method in (
            "hard_bin_half_open", "nearest_center", "bilinear_cic", "tsc_3x3",
            "gaussian_sigma_half_cell")}} for cluster in (
            "Abell2744", "MACS0416", "MACS1149", "AbellS1063", "Abell370")}}
        self.assertEqual(component_availability(audit),
                         {"kappa": True, "gamma1": False, "gamma2": False})


if __name__ == "__main__":
    unittest.main()
