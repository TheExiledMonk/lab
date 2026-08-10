from pbuf.matter.native_mass_loading_state import loading_inventory

def test_persistent_and_transient_classification_is_explicit():
    rows=loading_inventory()
    assert all(r['persistence'] in {'PERSISTENT','TRANSIENT','MIXED','UNDETERMINED'} for r in rows)
    assert any(r['persistence']=='PERSISTENT' for r in rows)

