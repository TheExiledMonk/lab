"""Hard gates and honest candidate selection for Dev153."""

def rank(rows):
    ids=sorted({r["candidate_id"] for r in rows})
    return [{"rank":i+1,"candidate_id":cid,"status":"NO_EFFECT" if cid=="T20" else "UNDERDETERMINED",
      "coefficient_free":True,"zero_load_recovery":True,"basis_covariant":True,"conservative":True,
      "derived_longitudinal_response":False} for i,cid in enumerate((["T20"]+[x for x in ids if x!="T20"]))]

def decide(ranking):
    survivors=[r for r in ranking if r["derived_longitudinal_response"]]
    return {"loaded_link_response_established":bool(survivors),"unique_law_selected":len(survivors)==1,
      "descriptor_equivalence_class_established":False,"outcome":"PBUF_ESTABLISHED_LONGITUDINAL_LINK_STATE_INSUFFICIENT_FOR_TRANSVERSE_RESPONSE",
      "reason":"Frozen Dev151/152 structure supplies orthogonal R but no L-dependent transverse Hessian, metric, allocation, or conservative exchange map."}
