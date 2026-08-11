import json
import numpy as np

from pbuf.excitation.native_vector_pair_dynamics import (
    VectorPairState, directed_relations, inverse_step, net_force, pair_forces,
    pair_power_flux, pair_reciprocity_error, positive_relations,
    relation_antisymmetry_error, source_contact_force, step,
)


def test_n6_relations_are_geometric_and_exactly_antisymmetric():
    rng=np.random.default_rng(1); u=rng.normal(0,1e-3,(5,6,7,3))
    r=directed_relations(u)
    assert r.shape==(5,6,7,6,3)
    assert relation_antisymmetry_error(u)==0
    assert np.all(np.linalg.norm(r,axis=-1)>0)


def test_central_force_reciprocity_homogeneous_null_and_label_independence():
    z=np.zeros((5,5,5,3)); assert np.array_equal(net_force(z),z)
    u=np.random.default_rng(2).normal(0,1e-3,z.shape)
    assert pair_reciprocity_error(u)==0
    # There is no loaded-label argument: identical relation input is identical output.
    assert np.array_equal(pair_forces(u),pair_forces(u.copy()))


def test_force_covariance_under_permutation_reflection_translation():
    u=np.random.default_rng(3).normal(0,1e-3,(5,5,5,3)); perm=(1,0,2)
    up=np.transpose(u,perm+(3,))[...,perm]
    expected=np.transpose(net_force(u),perm+(3,))[...,perm]
    assert np.max(np.abs(net_force(up)-expected))<1e-12
    ur=-np.flip(u,axis=(0,1,2))
    assert np.max(np.abs(net_force(ur)+np.flip(net_force(u),axis=(0,1,2))))<1e-12
    assert np.allclose(net_force(np.roll(u,2,0)),np.roll(net_force(u),2,0),atol=1e-15)


def test_kick_drift_exact_reverse_and_no_damping_or_averaging():
    rng=np.random.default_rng(4); u=rng.normal(0,1e-3,(5,5,5,3)); p=rng.normal(0,1e-3,u.shape)
    s=VectorPairState(u,p)
    for _ in range(20): s=step(s,.04)
    for _ in range(20): s=inverse_step(s,.04)
    assert np.max(np.abs(s.displacement-u))<1e-12
    assert np.max(np.abs(s.momentum-p))<1e-12


def test_local_source_contact_and_derived_pair_flux():
    shape=(7,7,7); f=source_contact_force(shape,(3,3,3))
    assert np.count_nonzero(np.linalg.norm(f,axis=-1))==6
    assert np.array_equal(np.sum(f,axis=(0,1,2)),np.zeros(3))
    u=np.random.default_rng(5).normal(0,1e-3,shape+(3,)); p=np.ones_like(u)
    j=pair_power_flux(u,p)
    assert j.shape==shape+(3,) and np.isfinite(j).all()


def test_dev167_runner_artifact_contract(tmp_path,monkeypatch):
    from tools import generate_dev167_pair_dynamics as lab
    monkeypatch.setattr(lab,"OUT",tmp_path)
    contract=lab.main()
    assert contract["DEV167_COMPLETE"] and contract["PRIMARY_TOPOLOGY"]=="N6"
    assert contract["PAIR_RELATION_ANTISYMMETRY_EXACT"]
    assert contract["PAIR_FORCE_RECIPROCITY_EXACT"]
    assert contract["REVERSIBILITY_ESTABLISHED"]
    required={"report.txt","repository_contract.json","historical_attempt_crosscheck.json",
      "pair_state_contract.json","pair_force_contract.json","reciprocity_audit.json",
      "homogeneous_equilibrium.json","symmetry_covariance.json","static_source_relational_state.json",
      "source_removal_response.json","moving_source_response.json","free_packet_response.json",
      "loaded_packet_response.json","loaded_unloaded_comparison.json","transverse_redirection_audit.json",
      "reversibility_audit.json","invariant_audit.json","pair_flux_audit.json","h07_comparison.json",
      "f02_f03_comparison.json","dev157_dispersion_comparison.json","final_contract.json","discussion_handoff.md"}
    assert required <= {p.name for p in tmp_path.iterdir()}
    audit=json.loads((tmp_path/"transverse_redirection_audit.json").read_text())
    assert abs(audit["C0_no_source_transverse"])<1e-12
    assert abs(audit["C1_centered_transverse_delta"])<1e-12
    assert abs(audit["reflection_odd_sum"])<1e-12
