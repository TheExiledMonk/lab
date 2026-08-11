import numpy as np

from pbuf.excitation.native_finite_receipt import NativeReceivedState
from pbuf.excitation.native_observer_adapter import adapt_native_receipt, execute_frozen_observer


def fixture_receipt(n=48):
    t=np.linspace(-2,2,n); z=np.zeros((n,3)); c=np.zeros((n,4))
    source=np.column_stack((np.ones(n),t,np.sin(t)))
    received=source+np.column_stack((np.full(n,4.),.1*np.sin(t),.05*np.cos(t)))
    direction=np.column_stack((np.ones(n),.02*np.sin(t),.01*np.cos(t)))
    direction/=np.linalg.norm(direction,axis=1,keepdims=True)
    return NativeReceivedState(source,received,direction,np.ones(n),np.arange(n),np.arange(n),z,z,direction,c,"BOND_FLUX")


def test_adapter_is_field_preserving_and_observer_has_45_channels():
    r=fixture_receipt(); a=adapt_native_receipt(r)
    np.testing.assert_array_equal(a["received_position_3d"],r.received_positions)
    np.testing.assert_array_equal(a["direction_3d"],r.directions)
    np.testing.assert_array_equal(a["deposition_weight"],r.weights)
    bank,meta=execute_frozen_observer(a,bins=3)
    assert len(bank)==45 and meta["primary_channel"] in bank
    assert any(np.isfinite(v).any() for v in bank.values())
