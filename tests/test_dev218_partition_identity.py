import json
from pathlib import Path

RUN = Path(__file__).parents[1] / 'runs/dev218_exact_interface_dynamic_polarity'


def test_dev217_artifacts_are_byte_identical_copies():
    partition = json.loads((RUN / 'dev217_partition_identity.json').read_text())
    bonds = json.loads((RUN / 'dev217_interface_identity.json').read_text())
    assert partition['DEV217_PARTITION_BYTE_IDENTICAL']
    assert bonds['DEV217_INTERFACE_BONDS_BYTE_IDENTICAL']
    assert bonds['bond_count'] == 121
