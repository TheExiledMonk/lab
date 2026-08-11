import ast
import json
from pathlib import Path

import numpy as np

from pbuf.data.hst_acs_common_frame import CommonFrame, combine_samples
from pbuf.data.hst_acs_raw_source import exposure_families
from pbuf.source.raw_detector_source_bridge import native_2d_constraint
from pbuf.labs.foundation import raw_abell2744_detector_to_native_source001 as dev161


def test_archive_has_116_complete_families():
    families=exposure_families(dev161.ARCHIVE)
    assert len(families)==116
    assert all(x.raw.exists() and x.flt.exists() and x.flc.exists() for x in families)


def test_target_blind_combination_retains_uncertainty_and_masks_dq():
    values=np.array([[2.,100.],[4.,6.]])
    errors=np.ones((2,2));good=np.array([[True,False],[True,True]])
    ra=np.array([[0.1,0.9],[0.1,0.9]]);dec=np.array([[0.1,0.1],[0.9,0.9]])
    out=combine_samples([(values,errors,good,ra,dec)],CommonFrame(0,1,0,1,(2,2)))
    assert out["coverage"].sum()==3
    assert np.isfinite(out["uncertainty"][out["occupied"]]).all()
    assert out["coverage"][0,1]==0


def test_native_constraint_is_relative_2d_not_rho3():
    image=np.array([[0.,1.],[2.,-1.]])
    sigma=np.ones((2,2))
    field=native_2d_constraint([image],[sigma],["F814W"])
    assert field.amplitude.shape==(1,2,2)
    assert np.isclose(field.amplitude.sum(),1)
    assert not field.depth_assigned and not field.physical_mass_scale
    assert np.all(field.amplitude>=0)


def test_runner_dependency_graph_is_fail_closed():
    audit=dev161._dependency_audit()
    assert audit["passed"]
    assert not audit["executable_forbidden_dependencies"]


def test_runner_contract_and_artifacts(tmp_path,monkeypatch):
    monkeypatch.setattr(dev161,"OUT",tmp_path)
    contract=dev161.main()
    assert contract["DEV161_AUDIT_COMPLETE"] is True
    assert contract["RAW_FILE_COUNT"]==contract["FLT_FILE_COUNT"]==contract["FLC_FILE_COUNT"]==116
    assert contract["PRIMARY_PIXEL_PRODUCT"]=="FLT"
    assert contract["COMMON_FRAME_ESTABLISHED"]=="TRUE"
    assert contract["RAW_TO_3D_SOURCE_UNIQUENESS"]=="NON_UNIQUE"
    assert contract["RAW_DATA_CAN_JUSTIFY_HISTORICAL_RHO3"]=="FALSE"
    assert contract["PREFERRED_RAW_NATIVE_ENDPOINT"]=="NATIVE_MULTI_CHANNEL_SOURCE_CONSTRAINT"
    for key in ("PREPROCESSED_LENSING_INPUT_USED","FIVE_CLUSTER_SOURCE_USED","MASS_MODEL_PRIOR_USED",
      "KAPPA_USED","GAMMA_USED","LENSING_TARGET_USED","ARBITRARY_DEPTH_EXTRUSION_USED",
      "NATIVE_LENS_GENERATED","DEV159_PROPAGATION_EXECUTED","OBSERVER_EXECUTED","OBSERVER_MODIFIED"):
        assert contract[key] is False
    required={"report.txt","archive_inventory.json","raw_flt_flc_role_audit.json","detector_metadata_inventory.json",
      "common_frame_contract.json","exposure_combination_audit.json","candidate_source_inventory.json",
      "historical_rho3_contract.json","depth_information_audit.json","source_support_diagnostics.json",
      "dev159_source_interface_contract.json","preferred_native_source_contract.json","downstream_validity_matrix.json",
      "final_raw_source_bridge_contract.json","native_2d_source_constraint.npz"}
    assert required<={p.name for p in tmp_path.iterdir()}
