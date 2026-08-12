import json
from pathlib import Path


def test_phase_b_is_frozen_and_categorical():
    root=Path(__file__).resolve().parents[1]; a=root/'runs/dev199_local_state_em_correlation/phase_a/phase_A_freeze.json'; b=root/'runs/dev199_local_state_em_correlation/phase_b/em_structural_correlation_result.json'
    assert a.exists() and b.exists()
    result=json.loads(b.read_text())
    assert result['PHASE_A_HASHES_VERIFIED'] and result['NO_PBUF_TO_SI_MAPPING'] and result['NO_QED_COEFFICIENT_FIT']
