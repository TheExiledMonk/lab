from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[1]

def test_dev202_contract_and_required_artifacts():
    out=ROOT/'runs/dev202_self_loaded_transverse'
    c=json.loads((out/'final_contract.json').read_text())
    assert c['BOND_TRANSVERSE_STIFFNESS_DERIVED']
    assert c['TWO_TRANSVERSE_PROBE_RESULTS_RECORDED']
    assert (out/'transverse_probe_1_trace.npz').exists()
