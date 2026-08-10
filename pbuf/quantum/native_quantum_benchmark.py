"""Read-only, post-freeze qualitative benchmark against effective QM relations."""
from __future__ import annotations
def compare(frozen_native: dict) -> dict:
    if not frozen_native.get("native_results_frozen"): raise RuntimeError("QM_BENCHMARK_REQUIRES_FROZEN_NATIVE_RESULTS")
    return {"QM_BENCHMARK_READ_ONLY":True,"PLANCK_RELATION_USED_TO_CONSTRUCT_NATIVE_MODE":False,
            "PLANCK_RELATION_USED_AS_POST_FREEZE_BENCHMARK":True,"HBAR_USED_IN_NATIVE_DERIVATION":False,
            "DE_BROGLIE_USED_TO_CONSTRUCT_NATIVE_LAW":False,"E_EQUALS_PC_USED_TO_CONSTRUCT_NATIVE_LAW":False,
            "planck_relation_status":"NOT_UNIQUE: total norm depends on amplitude and support, not k alone",
            "de_broglie_status":"MISSING_NATIVE_MECHANISM: directional norm flux has no unique k relation",
            "polarization_status":"POST_FREEZE_QM_COMPATIBLE"}
