import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'runs/dev229_persistent_native_source_derivation'
def load(n): return json.loads((OUT/n).read_text())
def test_repo_first_and_selection_are_frozen():
 c=load('final_contract.json'); assert c['CURRENT_GITHUB_INSPECTED'] and c['CURRENT_HEAD_VERIFIED']; assert c['DEV229_TEST_SELECTION']=='PERSISTENT_NATIVE_SOURCE_DERIVATION_GATE'; assert c['DEV229_TEST_SELECTION_FROZEN']
