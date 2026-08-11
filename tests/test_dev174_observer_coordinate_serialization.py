"""Regression checks for DEV174's provenance-only coordinate package."""
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "runs/dev174_observer_coordinate_serialization001"
D171 = ROOT / "runs/dev171_independent_3d_abell001"

def test_dev174_contract_and_primary_arrays_are_frozen():
    contract = json.loads((OUT / "final_contract.json").read_text())
    assert contract["DEV174_COMPLETE"] is True
    assert contract["PRIMARY_6X6_ARRAY_MODIFIED"] is False
    assert contract["ALL_8_REALIZATIONS_SERIALIZED"] is True
    hashes = json.loads((OUT / "frozen_dev171_hashes.json").read_text())
    for i in range(8):
        name = f"observer_realization_{i:02d}.npy"
        assert hashlib.sha256((D171 / name).read_bytes()).hexdigest() == hashes[name]
        sidecar = json.loads((OUT / f"observer_realization_{i:03d}.coordinate_provenance.json").read_text())
        assert sidecar["observer_shape"] == [6, 6]
        assert len(sidecar["bin_edges_u"]) == len(sidecar["bin_edges_v"]) == 7
        assert np.load(D171 / name).shape == (6, 6)

def test_dev174_retains_required_approximation_and_lineage_semantics():
    contract = json.loads((OUT / "final_contract.json").read_text())
    assert contract["GRID_TO_SKY_RECOVERABILITY"] == "DETERMINISTIC_APPROXIMATION"
    assert contract["FORMAL_FITS_WCS_CREATED"] is False
    assert contract["SOURCE_LINEAGE_SERIALIZED"] is True
    assert contract["RECEIPT_LINEAGE_SERIALIZED"] is True
    assert contract["DEPTH_SEMANTICS"] == "RELATIVE_ENSEMBLE_DEPTH"
    closure = json.loads((OUT / "coordinate_roundtrip_test.json").read_text())
    assert closure["ROUNDTRIP_CLOSURE_STATUS"] == "PASS_WITH_DISCRETIZATION"
