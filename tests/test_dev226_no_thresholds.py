from pathlib import Path

def test_observer_has_no_tolerance_or_threshold():
    source=(Path(__file__).parents[1]/'pbuf/analysis/native_staggered_order.py').read_text()
    assert 'values == 0' in source and 'np.isclose' not in source and 'np.abs(values)' not in source
