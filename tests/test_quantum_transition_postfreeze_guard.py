import pytest
from pbuf.quantum.native_transition_benchmark import compare
def test_benchmark_requires_freeze_and_is_read_only():
    with pytest.raises(ValueError): compare(False)
    r=compare(True); assert r['QM_TRANSITION_BENCHMARK_READ_ONLY'] and not r['PLANCK_RELATION_USED_TO_BUILD_TRANSITIONS']
