import json
from pathlib import Path
def test_unloaded_mismatch_is_classified_without_subtraction():
 d=json.loads((Path('runs/dev223_pattern_boundary_interface')/'unloaded_pattern_mismatch.json').read_text())
 assert d['UNLOADED_PATTERN_MISMATCH'] in {'ZERO','LATTICE_BASELINE'} and d['NO_BACKGROUND_FIT']
