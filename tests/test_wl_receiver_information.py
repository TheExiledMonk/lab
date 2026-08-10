import numpy as np
from pbuf.wl.receiver_information import classify_r2, linear_reconstruction_r2, rank_metrics


def test_rank_and_reconstruction():
    x=np.arange(30.,dtype=float);X=np.column_stack((x,x*x))
    m=rank_metrics(X)
    assert m["channel_count"]==2 and m["effective_rank"]>1
    np.testing.assert_allclose(linear_reconstruction_r2(x, x[:,None]),1)


def test_frozen_classifications():
    assert classify_r2(.49)=="STRONGLY_INDEPENDENT"
    assert classify_r2(.5)=="PARTIALLY_INDEPENDENT"
    assert classify_r2(.9)=="MOSTLY_REDUNDANT"
    assert classify_r2(.99)=="REDUNDANT"
