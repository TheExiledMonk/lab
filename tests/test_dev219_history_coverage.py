import json
from pathlib import Path
R=Path(__file__).parents[1]/'runs/dev219_magnetic_mechanism_wide_net'
def test_history():
 d=json.loads((R/'historical_magnetic_inventory.json').read_text()); assert d['NATIVE_MAGNETIC_HISTORY_COMPLETE'] and 'DEV218' in d['devs']
def test_failed_route_preservation():
 d=json.loads((R/'native_magnetic_failed_routes.json').read_text()); assert d['NO_RELABELED_FAILED_MECHANISM'] and len(d['rows'])>=5
