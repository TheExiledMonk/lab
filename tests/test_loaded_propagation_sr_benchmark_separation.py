import inspect
from pbuf.matter import native_loaded_propagation, native_mass_loading_state
from pbuf.matter.relativistic_loading_benchmark import benchmark_contract

def test_sr_is_post_freeze_and_not_imported_by_native_modules():
    assert 'relativistic_loading_benchmark' not in inspect.getsource(native_loaded_propagation)
    assert 'relativistic_loading_benchmark' not in inspect.getsource(native_mass_loading_state)
    c=benchmark_contract(); assert not c['SR_USED_TO_CONSTRUCT_PBUF_LAW'] and c['SR_USED_AS_POST_FREEZE_BENCHMARK']

