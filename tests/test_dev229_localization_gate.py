from test_dev229_repo_first import load
def test_noncompact_packet_is_not_thresholded():
 x=load('source_localization_representation.json'); assert x['SOURCE_LOCALIZATION_REPRESENTATION']=='NONUNIQUE' and x['NO_THRESHOLD_LOCALIZATION']
