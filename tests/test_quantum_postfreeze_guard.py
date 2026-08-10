import pytest
from pbuf.quantum.native_quantum_benchmark import compare
def test_qm_comparison_is_postfreeze_only():
    with pytest.raises(RuntimeError): compare({})
    assert compare({'native_results_frozen':True})['QM_BENCHMARK_READ_ONLY'] is True
