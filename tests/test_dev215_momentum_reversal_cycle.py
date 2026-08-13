import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; RUN=ROOT/'runs/dev215_lattice_state_cycle'
def test_prepared_momentum_reversal_is_not_time_reversal():
    a=json.loads((RUN/'momentum_reversal_cycle_direction.json').read_text())
    b=json.loads((RUN/'momentum_reversal_vs_time_reversal.json').read_text())
    assert a['MOMENTUM_REVERSAL_CYCLE_DIRECTION'] == 'NOT_APPLICABLE'
    assert b['MOMENTUM_REVERSAL_VS_TIME_REVERSAL'] == 'DISTINCT'
    assert b['inverse_step_endpoint_max_abs'] < 1e-13
