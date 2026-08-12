import json
from pathlib import Path


def test_turnaround_is_not_inferred_from_a_divergent_local_stress():
    run = Path(__file__).parents[1] / "runs/dev208_native_cosmic_turnaround"
    assert json.loads((run / "finite_extension_stress_behavior.json").read_text())["FINITE_EXTENSION_STRESS_BEHAVIOR"] == "DIVERGENT_AT_BOUND"
    assert json.loads((run / "native_cosmological_turnaround.json").read_text())["NATIVE_COSMOLOGICAL_TURNAROUND"] == "BLOCKED"
