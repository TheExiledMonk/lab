import numpy as np

from pbuf.excitation.native_relational_dynamics import f03_invariant
from pbuf.source.native_moving_source import evolve_schedule, integer_schedule, release
from pbuf.source.native_source_medium_interaction import equilibrium_residual, stationary_response
from pbuf.source.native_source_state import NativeSourceState


def test_one_cell_source_has_exact_stationary_relational_equilibrium():
    shape=(11,11,11); source=NativeSourceState((5,5,5),.02)
    q=stationary_response(shape,source)
    assert abs(np.mean(q)) < 1e-15
    assert np.max(abs(equilibrium_residual(q,source))) < 1e-13
    assert np.count_nonzero(abs(q)>1e-14) == np.prod(shape)


def test_stationary_response_is_reflection_permutation_and_amplitude_covariant():
    shape=(11,11,11); q=stationary_response(shape,NativeSourceState((5,5,5),.01))
    assert np.allclose(q,np.flip(q,(0,1,2)))
    assert np.allclose(q,np.transpose(q,(2,1,0)))
    q2=stationary_response(shape,NativeSourceState((5,5,5),.02))
    assert np.allclose(q2,2*q)


def test_same_translated_interaction_generates_residual_without_packet():
    shape=(11,11,11); start=(3,5,5); amplitude=.01
    q0=stationary_response(shape,NativeSourceState(start,amplitude))
    schedule=integer_schedule(start,0,3,1)
    result=evolve_schedule(shape,amplitude,schedule,q0=q0)
    assert np.max(abs(result["dynamic_residual"])) > 0
    assert result["states"].shape[0] == len(schedule)


def test_source_removal_uses_exact_conservative_f03_dynamics():
    shape=(11,11,11); q=stationary_response(shape,NativeSourceState((5,5,5),.01)); r=np.zeros_like(q)
    initial=f03_invariant(q,r); states=release(q,r,12)
    # reconstruct retained state: q_n-q_(n-1)
    retained=states[-1]-states[-2]
    assert np.isclose(f03_invariant(states[-1],retained),initial,rtol=1e-12,atol=1e-12)


def test_dev159_runner_writes_all_required_artifacts(tmp_path,monkeypatch):
    from pbuf.labs.foundation import native_source_medium_interaction001 as lab
    monkeypatch.setattr(lab,"RUN",tmp_path); contract=lab.main()
    required=("report.txt source_candidate_inventory.json stationary_source_response.json static_lane_compatibility.json "
      "source_removal_response.json moving_source_response.json moving_source_spectrum.json dispersion_match.json "
      "source_magnitude_scaling.json source_work_audit.json observer_handoff_contract.json downstream_validity_matrix.json "
      "final_source_medium_contract.json").split()
    assert all((tmp_path/name).exists() for name in required)
    assert contract["DEV159_AUDIT_COMPLETE"] and not contract["DEV156_LAWS_MODIFIED"]
