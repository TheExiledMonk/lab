from pathlib import Path

from pbuf.wl.hst_acs_reference_closure import acquire_references


def test_offline_missing_reference(tmp_path: Path):
    result = acquire_references(["jref$does_not_exist.fits"], tmp_path, offline=True)
    assert not result["reference_closure_complete"]
    assert result["missing_references"] == ["jref$does_not_exist.fits"]


def test_offline_corrupt_fits(tmp_path: Path):
    (tmp_path / "bad_idc.fits").write_bytes(b"not fits")
    result = acquire_references(["jref$bad_idc.fits"], tmp_path, offline=True)
    assert not result["reference_closure_complete"]
    assert result["records"][0]["status"] == "CORRUPT_FITS"
