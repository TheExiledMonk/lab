from pbuf.foundation.native_loaded_link_discriminator import rank,decide
from pbuf.foundation.native_loaded_link_response import execute
def test_null_is_clean_not_parameterized():
    rows=[execute(f"T{i:02d}","LOAD03","EX01") for i in range(1,21)]; d=decide(rank(rows))
    assert not d["loaded_link_response_established"] and d["outcome"]=="PBUF_ESTABLISHED_LONGITUDINAL_LINK_STATE_INSUFFICIENT_FOR_TRANSVERSE_RESPONSE"
