from pathlib import Path
from pbuf.foundation.native_neighbor_survivor_registry import load_dev151_survivors, validate_dev151

ROOT=Path(__file__).resolve().parents[1]
def test_registry_is_complete_and_validated():
    rows=load_dev151_survivors()
    assert [r.constitutive_law_id for r in rows]==["C08","C10","C12","C13","C16","C18"]
    assert validate_dev151(ROOT)["validated"]

