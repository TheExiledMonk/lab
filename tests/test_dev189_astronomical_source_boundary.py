import json, subprocess, sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'runs/dev189_astronomical_source_boundary'
sys.path.insert(0,str(ROOT))
from pbuf.observer.astronomical_source_boundary import SimilarityBridge, normalized_second_moment

def test_dev189_contract_is_a_clean_blocker():
    c=json.loads((OUT/'final_contract.json').read_text())
    assert c['ASTRONOMICAL_SOURCE_BOUNDARY']=='BLOCKED_INTRINSIC_SOURCE_UNAVAILABLE'
    for k in ('OBSERVED_LENSED_IMAGE_NOT_USED_AS_INCIDENT_SOURCE','NO_GR_DISTANCE_CONVERSION','NO_FITTED_SCALE','NO_BILINEAR_INTERPOLATION','DEV188_TRANSFER_OPERATOR_UNCHANGED'): assert c[k]
    assert c['SPIN2_OBSERVABLE_GATE']=='CLOSED' and c['OBSERVATIONAL_COMPARISON_GATE']=='CLOSED'

def test_similarity_control_and_native_hash():
    r=np.array([[0.,-1.],[1.,0.]]) ; b=SimilarityBridge(2.,r,np.array([1.,2.]))
    p=np.array([[-1.,0.],[0.,0.],[2.,1.]]) ; w=np.array([1.,3.,2.])
    assert not b.injects_anisotropy()
    assert np.allclose(np.linalg.eigvalsh(normalized_second_moment(p,w)),np.linalg.eigvalsh(normalized_second_moment(b.transform(p),w)))
    assert json.loads((OUT/'dev188_native_domain_verification.json').read_text())['DEV188_NATIVE_DOMAIN_HASH_VERIFIED']

def test_dev189_rerun_is_deterministic():
    before=(OUT/'astronomical_source_boundary.json').read_bytes()
    subprocess.check_call([sys.executable,'tools/generate_dev189_astronomical_source_boundary.py'],cwd=ROOT)
    assert (OUT/'astronomical_source_boundary.json').read_bytes()==before
