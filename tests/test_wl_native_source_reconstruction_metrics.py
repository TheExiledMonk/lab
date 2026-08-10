import numpy as np
from pbuf.wl.native_source_reconstruction_metrics import ambiguity_area, information_gain, score_depth

def test_metric_sanity_and_multivalued():
    assert ambiguity_area([[0,1],[2,3]]) == .25
    assert information_gain(.25,.5) == .5
    p={"primary_depth":2.,"depth_candidates":[2.,3.],"support_interval":[1.9,2.1]}
    assert score_depth(p,3.,1.)["classification"] == "CORRECT_MULTIVALUED"

def test_false_unique():
    p={"primary_depth":2.,"depth_candidates":[2.],"support_interval":[1.9,2.1]}
    assert score_depth(p,3.,1.)["classification"] == "FALSE_UNIQUE_DEPTH"
