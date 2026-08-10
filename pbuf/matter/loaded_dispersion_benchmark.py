"""External, post-freeze massive-wave/SR morphology comparison only."""
from __future__ import annotations

def benchmark_contract(native_dispersion_established: bool=False):
    return {"lane":"POST_FREEZE_EXTERNAL_COMPARISON","native_dispersion_frozen":True,
            "SR_USED_TO_CONSTRUCT_NATIVE_LAW":False,"SR_DISPERSION_USED_TO_CONSTRUCT_NATIVE_LAW":False,
            "KLEIN_GORDON_USED_TO_CONSTRUCT_NATIVE_LAW":False,
            "comparison":"NOT_COMPARABLE" if not native_dispersion_established else "PENDING_COMPARISON",
            "relativistic_structure_emerges":False}
