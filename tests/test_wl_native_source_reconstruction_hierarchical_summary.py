from pbuf.wl.native_source_reconstruction_metrics import hierarchical_summary

def test_duplicate_within_family_does_not_reweight_global_summary():
    a={"morphology":"a","lens_family":"l","source_size":1,"source_depth":2,"response_regime":"M","x":0}
    b={**a,"morphology":"b","x":10}
    before=hierarchical_summary([a,b],"x")["median"]
    after=hierarchical_summary([a,a,a,b],"x")["median"]
    assert before==after==5
