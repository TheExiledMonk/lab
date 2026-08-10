from pbuf.wl.native_source_reconstruction_sweep import blind_reconstruct, sha256, synthetic_observation, trial_matrix

def test_resume_prediction_sha_is_deterministic():
    rows=trial_matrix(validation=True)
    def run(items):
        out=[]
        for row in items:
            p=blind_reconstruct(synthetic_observation(row),"C4");p.pop("score_surface")
            out.append({"trial_id":row["trial_id"],**p})
        return sha256(sorted(out,key=lambda x:x["trial_id"]))
    assert run(rows)==run(list(reversed(rows)))
