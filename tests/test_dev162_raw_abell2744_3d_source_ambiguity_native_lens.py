import numpy as np

from pbuf.source.projected_source_3d_family import diagnostic_family,project
from pbuf.lens.native_stationary_lens_from_source import stationary_distributed_response,equilibrium_residual
from pbuf.labs.foundation import raw_abell2744_3d_source_ambiguity_native_lens001 as dev162


def test_projection_equivalent_distinct_family():
    image=np.arange(1,43,dtype=float).reshape(6,7);image/=image.sum();rows=diagnostic_family(image)
    assert len(rows)>=5 and {r.family for r in rows}>={"THIN","FINITE_SYMMETRIC","DOUBLE_LAYER","ASYMMETRIC","SPATIALLY_VARIABLE_DEPTH"}
    assert all(np.allclose(project(r.source),image,rtol=0,atol=1e-15) for r in rows)
    assert all(np.isclose(r.source.sum(),1) and not r.physical_truth_claimed for r in rows)


def test_distributed_stationary_response_reaches_equilibrium_and_is_linear():
    image=np.zeros((8,8));image[2,3]=.4;image[5,4]=.6
    source=diagnostic_family(image)[4].source;q=stationary_distributed_response(source)
    assert np.max(np.abs(equilibrium_residual(q,source)))<1e-10
    assert np.allclose(stationary_distributed_response(4*source),4*q,rtol=1e-12,atol=1e-12)


def test_runner_contract_and_artifacts(tmp_path,monkeypatch):
    monkeypatch.setattr(dev162,"OUT",tmp_path);contract=dev162.main()
    assert contract["DEV162_AUDIT_COMPLETE"] and contract["DEV161_SOURCE_REUSED"]
    assert contract["3D_SOURCE_REALIZATION_COUNT"]>=5 and contract["ALL_3D_REALIZATIONS_PROJECTION_EQUIVALENT"]
    assert contract["DEV159_STATIC_SOURCE_INTERACTION_USED"] and not contract["DEV159_DYNAMIC_PROPAGATION_USED"]
    for key in ("KAPPA_USED","GAMMA_USED","EXTERNAL_MASS_MAP_USED","EXTERNAL_DEPTH_INFORMATION_USED","PHYSICAL_DEPTH_SCALE_ASSUMED",
                "LENSING_PROPAGATION_EXECUTED","OBSERVER_EXECUTED","OBSERVER_MODIFIED"):
        assert contract[key] is False
    required={"report.txt","filter_morphology_consistency.json","common_morphology_contract.json","source_3d_realization_inventory.json",
      "projection_equivalence.json","stationary_native_lens_results.json","native_lens_support_metrics.json","transverse_lens_invariance.json",
      "full_3d_lens_degeneracy.json","amplitude_geometry_separation.json","filter_lens_stability.json","lensing_handoff_contract.json",
      "downstream_validity_matrix.json","final_3d_ambiguity_native_lens_contract.json"}
    assert required<={p.name for p in tmp_path.iterdir()}
