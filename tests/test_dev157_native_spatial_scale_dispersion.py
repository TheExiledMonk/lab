import numpy as np

from pbuf.excitation.native_dispersion_observer import (analytic_progression_frequency,
    measure_mode_frequency, mode, wavevector)
from pbuf.excitation.native_spatial_spectrum import native_wavevectors, reconstruct, spectrum3d
from pbuf.excitation.native_spatial_support import radial_correlation, support_metrics


def test_3d_spectrum_roundtrip_parseval_and_explicit_dc():
    rng=np.random.default_rng(157); q=rng.normal(size=(9,7,5)); result=spectrum3d(q)
    assert np.allclose(reconstruct(result["transform"]),q)
    assert np.isclose(result["power"].sum(),np.sum(q*q))
    assert np.isclose(result["dc_power"],abs(result["transform"][0,0,0])**2)


def test_native_wave_number_and_dc_contract():
    axes=native_wavevectors((8,9,10))
    assert axes[0][1]==2*np.pi/8 and axes[0][4]==-np.pi


def test_f02_f03_share_exact_measured_dispersion():
    shape=(17,17,17); indices=(2,1,0); exact=analytic_progression_frequency(wavevector(shape,indices))
    for representation in ("F02","F03"):
        row=measure_mode_frequency(representation,shape,indices,steps=18)
        assert np.isclose(row["progression_frequency"],exact,atol=2e-13)


def test_dispersion_reflection_permutation_and_amplitude_covariance():
    shape=(17,17,17)
    rows=[measure_mode_frequency("F03",shape,m,a) for m in ((2,0,0),(0,2,0),(15,0,0)) for a in (1e-4,1e-2)]
    assert np.ptp([x["progression_frequency"] for x in rows]) < 2e-13


def test_two_mode_superposition_remains_linear():
    shape=(17,17,17); q=mode(shape,(1,0,0),1e-3)+mode(shape,(0,2,0),2e-3)
    from pbuf.excitation.native_relational_dynamics import f03_step
    r=np.zeros_like(q)
    for _ in range(8): q,r=f03_step(q,r)
    # No Fourier support appears outside the two real-mode pairs beyond roundoff.
    power=np.abs(np.fft.fftn(q))
    mask=np.ones(shape,bool)
    for idx in ((1,0,0),(16,0,0),(0,2,0),(0,15,0)): mask[idx]=False
    assert power[mask].max() < 1e-11


def test_support_and_signed_correlation_are_threshold_independent():
    q=np.zeros((9,9,9)); q[4,4,4]=1; q[5,4,4]=1
    metrics=support_metrics(q*q,(4,4,4)); correlation=radial_correlation(q)
    assert np.isclose(metrics["rms_radius"],np.sqrt(.5))
    assert metrics["participation_volume"]==2
    assert correlation["correlation"][0]==1
