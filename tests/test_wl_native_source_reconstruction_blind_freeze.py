import json,pytest
from pbuf.wl.native_source_reconstruction_sweep import TruthVault

def test_truth_cannot_load_before_prediction_freeze(tmp_path):
    sealed=tmp_path/"truth.json"; sealed.write_text(json.dumps({"truth":[{"z":2}]}))
    vault=TruthVault(sealed,tmp_path/"blind_prediction_manifest.json")
    with pytest.raises(PermissionError,match="DEV139_BLINDNESS_FAILURE"): vault.load()
    (tmp_path/"blind_prediction_manifest.json").write_text("{}")
    assert vault.load()==[{"z":2}]
