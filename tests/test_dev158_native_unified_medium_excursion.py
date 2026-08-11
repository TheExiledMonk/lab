import numpy as np

from pbuf.excitation.native_bond_state import positive_gradient
from pbuf.excitation.native_content_density import positivity_audit
from pbuf.excitation.native_dynamic_constitutive_audit import bounded_f03_inverse, bounded_f03_step
from pbuf.excitation.native_excursion_bridge import native_bond_excursion, static_small_excursion
from pbuf.excitation.native_relational_dynamics import f02_invariant


def test_static_small_excursion_has_frozen_quadratic_limit():
    e=np.logspace(-7,-3,8); result=static_small_excursion(e)
    assert np.allclose(result["exact"]/result["quadratic"],1,rtol=1e-6)
    assert np.all(result["remainder"]>=0)


def test_common_excursion_is_existing_relational_gradient():
    rng=np.random.default_rng(158); q=rng.normal(size=(7,7,7))
    assert np.array_equal(native_bond_excursion(q),positive_gradient(q))


def test_isolated_bounded_candidate_is_reversible_while_admissible():
    q=np.zeros((9,9,9)); q[4,4,4]=1e-3; r=np.zeros_like(q)
    q1,r1=bounded_f03_step(q,r); q0,r0=bounded_f03_inverse(q1,r1)
    assert np.allclose(q0,q) and np.allclose(r0,r)


def test_exact_local_summand_does_not_force_positivity():
    rng=np.random.default_rng(9); q=rng.normal(size=(7,7,7)); b=-.5*positive_gradient(q)
    result=positivity_audit("F02",q,b)
    assert np.isclose(result["sum"],f02_invariant(q,b))
    assert result["negative_site_count"]>0 and not result["pointwise_positive_proven"]


def test_dev158_runner_writes_required_contract(tmp_path, monkeypatch):
    from pbuf.labs.foundation import native_unified_medium_excursion001 as lab
    monkeypatch.setattr(lab,"RUN",tmp_path); contract=lab.main()
    required=("report.txt static_small_excursion_analysis.json dynamic_invariant_analysis.json "
      "common_excursion_mapping.json static_dynamic_scaling.json equal_invariant_mode_comparison.json "
      "bounded_dynamic_response_results.json local_positive_density_audit.json observer_handoff_contract.json "
      "downstream_validity_matrix.json final_unified_excursion_contract.json").split()
    assert all((tmp_path/name).exists() for name in required)
    assert contract["DEV158_AUDIT_COMPLETE"] and not contract["ARBITRARY_COUPLING_INTRODUCED"]
