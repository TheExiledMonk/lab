import json
from pathlib import Path

RUN = Path(__file__).parents[1] / 'runs/dev218_exact_interface_dynamic_polarity'


def test_direct_interface_is_reciprocal_in_every_state():
    data = json.loads((RUN / 'direct_interface_action_reaction.json').read_text())
    assert data['DIRECT_INTERFACE_ACTION_REACTION'] in {'EXACT', 'ROUND_OFF'}
    assert set(data['rows']) == {'pp', 'pm', 'mp', 'mm'}
