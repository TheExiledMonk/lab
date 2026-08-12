import json
from pathlib import Path


def test_cycle_and_big_bang_claims_remain_unestablished():
    run = Path(__file__).parents[1] / "runs/dev208_native_cosmic_turnaround"
    assert json.loads((run / "native_cyclic_cosmology.json").read_text())["NATIVE_CYCLIC_COSMOLOGY"] == "NOT_DERIVED"
    assert json.loads((run / "big_bang_as_previous_bounce.json").read_text())["BIG_BANG_AS_PREVIOUS_BOUNCE"] == "NOT_DERIVED"
