import json
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'runs/dev230_native_emission_dependency_split'
def load(name): return json.loads((OUT / name).read_text())
def test_repo_first():
    c = load('final_contract.json')
    for key in ('DEV229_PUBLICATION_VERIFIED', 'CURRENT_GITHUB_INSPECTED', 'CURRENT_HEAD_VERIFIED', 'MECHANISM_REGISTRY_QUERIED', 'DEVELOPMENT_LEDGER_READ', 'HISTORICAL_INDEX_READ', 'DERIVATION_GRAPH_READ'):
        assert c[key]
