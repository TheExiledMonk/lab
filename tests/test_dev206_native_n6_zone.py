import json
from pathlib import Path
import numpy as np
from pbuf.observer.native_n6_zone import zone_tensor
from pbuf.observer.relational_wave_content import directional_tensor
from pbuf.observer.native_zone_em_mapping import mapping_status

def test_fixed_n6_tensor_equals_dev203_operator():
    x=np.arange(18.,dtype=float).reshape(6,3)
    assert np.array_equal(zone_tensor(x),directional_tensor(x))

def test_dev206_preserves_partial_boundary():
    out=Path(__file__).resolve().parents[1]/'runs/dev206_n6_zone_radiative_sector'
    c=json.loads((out/'final_contract.json').read_text())
    assert c['DEV205_MAPPING_RESULT_PRESERVED'] and c['ZONE_TENSOR_EQUIVALENCE_TO_DEV203_CLASSIFIED']
    assert c['OUTGOING_SECTOR_DEFINITION']=='BLOCKED_NO_FREE_NATIVE_PROPAGATION_WINDOW'

def test_mapping_cannot_promote_without_full_native_geometry():
    assert mapping_status('PARTIAL','BIDIRECTIONAL','NOT_DERIVED') == 'PARTIAL'
