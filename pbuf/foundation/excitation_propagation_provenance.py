"""Static and runtime provenance inventory for the established PBUF propagator."""
from __future__ import annotations

from pathlib import Path

STATE_FIELDS = ("variable_name", "module", "function", "producer", "consumer", "shape",
 "scalar_vector_tensor", "native_units_or_dimensionless", "signed", "positive_definite", "normalized",
 "normalization_location", "persistent_across_steps", "stored_in_history", "stored_at_receiver",
 "derived_from_medium", "derived_from_source", "purely_geometric", "candidate_excitation_state",
 "candidate_excitation_magnitude", "discarded", "discard_location")


def _row(name, module, function, producer, consumer, shape, kind, **kw):
    defaults = dict(native_units_or_dimensionless="native", signed=True, positive_definite=False,
        normalized=False, normalization_location=None, persistent_across_steps=False,
        stored_in_history=False, stored_at_receiver=False, derived_from_medium=False,
        derived_from_source=False, purely_geometric=False, candidate_excitation_state=False,
        candidate_excitation_magnitude=False, discarded=False, discard_location=None)
    defaults.update(kw)
    row = dict(variable_name=name, module=module, function=function, producer=producer,
               consumer=consumer, shape=shape, scalar_vector_tensor=kind, **defaults)
    assert set(row) == set(STATE_FIELDS)
    return row


def state_inventory():
    geo = "pbuf.labs.foundation.los_consistent_ray_geometry001"
    return [
      _row("position",geo,"_propagate_g3d","RayLaunch x0/y0 and recurrence","field sampler/receiver","(rays,3)","vector",persistent_across_steps=True,stored_in_history=True,stored_at_receiver=True,derived_from_source=True,purely_geometric=True),
      _row("direction",geo,"_propagate_g3d","fixed (0,0,1) launch and normalized recurrence","position recurrence/receiver","(rays,3)","vector",native_units_or_dimensionless="dimensionless",normalized=True,normalization_location="_propagate_g3d",persistent_across_steps=True,stored_in_history=True,stored_at_receiver=True,purely_geometric=True),
      _row("rx_sample/ry_sample",geo,"_sample","frozen response field","raw direction update","(rays,2)","vector",derived_from_medium=True,stored_in_history=True),
      _row("v_raw",geo,"_propagate_g3d","direction + path_step*response","normalization","(rays,3)","vector",derived_from_medium=True,purely_geometric=True,candidate_excitation_state=True,discarded=True,discard_location="immediately after normalization"),
      _row("|v_raw|",geo,"_propagate_g3d","Euclidean norm of v_raw","normalization divisor","(rays,)","scalar",positive_definite=True,derived_from_medium=True,purely_geometric=True,candidate_excitation_magnitude=True,discarded=True,discard_location="immediately after normalization"),
      _row("path_step",geo,"_propagate_g3d","PropagationConfig","direction and position recurrence","scalar","scalar",native_units_or_dimensionless="native path length",positive_definite=True,persistent_across_steps=True,purely_geometric=True),
    ]


def source_state_inventory():
    return [
      {"name":"initial position","present":True,"source":"RayLaunch.x0/y0","classification":"GEOMETRIC"},
      {"name":"initial direction","present":True,"source":"hard-coded (0,0,1)","classification":"GEOMETRIC"},
      *[{"name":n,"present":False,"source":None,"classification":"ABSENT"} for n in
        ("initial path state","initial response-channel state","source weight","source amplitude-like fields",
         "ray weight","bundle weight","event weight","hidden source scalar")],
    ]


def call_graph():
    return {"nodes":["RayLaunch","_propagate_g3d","_sample","normalize","checkpoint","receiver","ArrivalEvent2D"],
            "edges":[["RayLaunch","_propagate_g3d"],["_propagate_g3d","_sample"],["_sample","normalize"],
                     ["normalize","checkpoint"],["checkpoint","receiver"],["receiver","ArrivalEvent2D"]],
            "dynamic_non_geometric_edge":False}


def candidate_manifest():
    names=("raw trajectory update vector","raw update magnitude","longitudinal update component",
      "transverse update component","fast-channel amplitude","slow-channel amplitude","combined fast/slow state",
      "pair-transfer state","link-state norm","per-axis modal state","source weight","ray weight","bundle weight",
      "accumulated response integral","curvature-integrated state","path-excess state",
      "receiver-preserved interaction state","discarded pre-normalization state",
      "hidden dynamic state elsewhere in propagation stack","no existing dynamic excitation state")
    rows=[]
    for i,name in enumerate(names,1):
        status = ("NUMERICAL_UPDATE_STATE" if i in (1,2,3,4,18) else
                  "STATIC_MEDIUM_STATE" if i in (5,6,7,8,9,10) else
                  "OBSERVATIONAL_WEIGHT_ONLY" if i in (11,12,13) else
                  "DYNAMIC_STATE_GEOMETRIC_ONLY" if i in (14,15,16,17) else
                  "MISSING_PROVENANCE" if i == 19 else "NOT_APPLICABLE")
        rows.append({"id":f"E{i:02d}","name":name,"attempted":True,"status":status,
                     "survives_dynamic_excitation_gate":False})
    return rows


def provenance_contract():
    return {"contract":"PBUF_EXISTING_EXCITATION_PROVENANCE_V1",
      "trajectory_semantics":"GEOMETRIC_TRACER_ONLY","dynamic_excitation_state_found":False,
      "dynamic_excitation_definition":None,"state_location":None,"source_initialized":False,
      "step_persistent":False,"link_based":False,"node_based":False,"trajectory_attached":False,
      "raw_magnitude_exists":True,"raw_magnitude_physical":False,
      "raw_magnitude_definition":"|n + path_step*(rx,ry,0)|","longitudinal_component_exists":True,
      "transverse_component_exists":True,"fast_slow_dynamic":False,"fast_slow_static":True,
      "normalization_discards_physical_state":False,"receiver_preserves_excitation_state":False,
      "new_dynamic_dof_required":True}


def repository_module_inventory(root=None):
    root = Path(root or Path(__file__).resolve().parents[2])
    patterns=("propagat","traject","fast","slow","receiver","arrival","bundle","response")
    rows=[]
    for path in sorted((root/"pbuf").rglob("*.py")):
        text=path.read_text(errors="replace").lower()
        hits=[p for p in patterns if p in text or p in path.name.lower()]
        if hits: rows.append({"module":str(path.relative_to(root)).replace("/",".")[:-3],"families":hits})
    return rows
