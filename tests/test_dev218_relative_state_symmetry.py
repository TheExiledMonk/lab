import json
from pathlib import Path

RUN = Path(__file__).parents[1] / 'runs/dev218_exact_interface_dynamic_polarity'


def test_reversal_symmetries_are_classified():
    for name, key in [('pp_mm_force_symmetry.json', 'PP_MM_FORCE_SYMMETRY'), ('pm_mp_force_symmetry.json', 'PM_MP_FORCE_SYMMETRY')]:
        assert json.loads((RUN / name).read_text())[key] in {'EXACT', 'ROUND_OFF', 'BROKEN'}
