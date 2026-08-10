"""Read-only post-freeze benchmark; never constructs native transitions."""
def compare(native_frozen):
    if not native_frozen: raise ValueError("native transition result must be frozen before benchmark")
    return {"QM_TRANSITION_BENCHMARK_READ_ONLY":True,"PLANCK_RELATION_USED_TO_BUILD_TRANSITIONS":False,
            "PLANCK_RELATION_USED_AS_POST_FREEZE_BENCHMARK":True,"HBAR_USED_TO_FIT_TRANSITION_SCALE":False,
            "planck_relation_status":"NOT_COMPARABLE_MISSING_NATIVE_TRANSITIONS"}
