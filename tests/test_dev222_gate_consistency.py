import json
from pathlib import Path

from pbuf.registry.validate import validate

ROOT = Path(__file__).parents[1]


def test_blocked_pair_gate_cannot_select_pair_operation():
    registry = json.loads((ROOT / "docs/PBUF_MECHANISM_REGISTRY.json").read_text())
    assert not validate(registry)
    selection = json.loads((ROOT / "runs/dev222_dev221_reconciliation/dev223_test_selection.json").read_text())
    assert selection["DEV223_TEST_SELECTION"] != "PAIR_ORIENTATION_INTERACTION"
