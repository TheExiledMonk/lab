from pbuf.audit.native_provenance_map import provenance
def test_every_mechanism_has_dev_provenance():
    assert all(x["first_introduced"] and x["latest_verifying_dev"]==154 for x in provenance())
