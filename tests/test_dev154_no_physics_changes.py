from pathlib import Path
def test_dev154_modules_are_audit_only():
    root=Path(__file__).parents[1]; lab=(root/"pbuf/labs/audit/native_microphysics_reconstruction001.py").read_text()
    assert "source_hashes_unchanged" in lab
    assert all(p.parent.name=="audit" for p in (root/"pbuf/audit").glob("*.py"))
