from test_dev230_repo_first import load
def test_future_lanes_preserved():
    d=load('future_lane_dependencies.json'); assert all(d[k] for k in ('PERSISTENT_SOURCE_LANE_PRESERVED','INTERSTITIAL_PATTERN_LANE_PRESERVED','COLLECTIVE_X_BODY_LANE_PRESERVED','THREE_CONSTITUENT_LANE_PRESERVED','N6_N27_FUTURE_LANE_PRESERVED'))
