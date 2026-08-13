from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_dev227_is_a_documentary_audit_only():
    source = (ROOT / "tools/generate_dev227_magnetic_candidate_exhaustion.py").read_text()
    assert "pbuf/" not in source
    assert "NO_NEW_FORCE" in source
    assert "NO_NEW_DOF" in source
    assert "NO_NEW_MAGNETIC_FIELD" in source
