import json
from pathlib import Path

RUN = Path(__file__).parents[1] / 'runs/dev218_exact_interface_dynamic_polarity'


def test_four_frozen_initial_force_classes_are_complete():
    data = json.loads((RUN / 'initial_interface_force_class.json').read_text())
    assert set(data['INITIAL_INTERFACE_FORCE_CLASS']) == {'pp', 'pm', 'mp', 'mm'}
