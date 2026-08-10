import inspect
from pbuf.matter import native_excitation_state, native_excitation_invariants, native_loaded_excitation_modes
from pbuf.matter.loaded_dispersion_benchmark import benchmark_contract

def test_external_dispersion_does_not_construct_native_modules():
    for module in (native_excitation_state,native_excitation_invariants,native_loaded_excitation_modes):
        assert 'loaded_dispersion_benchmark' not in inspect.getsource(module)
    c=benchmark_contract(False)
    assert c['comparison']=='NOT_COMPARABLE' and not c['KLEIN_GORDON_USED_TO_CONSTRUCT_NATIVE_LAW']
