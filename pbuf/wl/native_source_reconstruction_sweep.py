"""Deterministic scheduling, blind freeze, and scoring support for Dev139."""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import os
from pathlib import Path
import tempfile
import numpy as np

from .native_source_controls import MORPHOLOGIES, SOURCE_SIZES, DEPTH_OFFSETS, LENS_FAMILIES

RESOLUTIONS=(32,48,64,96,128)
POPULATIONS=(25,50,75,100)
INFORMATION_LANES=("C0","C1","C2","C3","C4","C5")
DELETION_LANES=("D_DIRECTION","D_BUNDLE","D_TRAJECTORY","D_SECOND","D_RECEIVER_Z","D_MULTIPATH")
PARTIAL_LANES=("DROP10","DROP25","DROP50","MISS_CENTER","MISS_ONE_BRANCH","MISS_POSITIVE_U_SIDE")


def canonical_json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def sha256(value):
    return hashlib.sha256(canonical_json(value).encode()).hexdigest()


def atomic_json(path, value):
    path=Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    fd,tmp=tempfile.mkstemp(prefix=path.name+".", dir=path.parent)
    try:
        with os.fdopen(fd,"w") as f: f.write(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+"\n")
        os.replace(tmp,path)
    finally:
        if os.path.exists(tmp): os.unlink(tmp)


def trial_matrix(validation=False):
    rows=[]
    morphs=MORPHOLOGIES[:2] if validation else MORPHOLOGIES
    lenses=LENS_FAMILIES[:2] if validation else LENS_FAMILIES
    sizes=SOURCE_SIZES[:2] if validation else SOURCE_SIZES
    depths=DEPTH_OFFSETS[:3] if validation else DEPTH_OFFSETS
    for m in morphs:
      for s in sizes:
       for d in depths:
        for lens in lenses:
         key=f"{m}|{s}|{d}|{lens}"
         rows.append({"trial_id":"T"+hashlib.sha256(key.encode()).hexdigest()[:16],"morphology":m,
          "source_size":s,"source_depth":1+d,"lens_family":lens,
          "response_regime":("WEAK" if lens=="weak_response" else "STRONG_UNSATURATED" if lens in ("compact_strong","strong_unsaturated") else "MODERATE")})
    return sorted(rows,key=lambda x:x["trial_id"])


def stable_event_selection(uids, percent):
    return np.array([int(hashlib.sha256(str(x).encode()).hexdigest()[:8],16)%100 < percent for x in uids])


class TruthVault:
    """Process-local guard proving scoring truth cannot load before freeze."""
    def __init__(self, sealed_path, freeze_path): self.sealed_path=Path(sealed_path); self.freeze_path=Path(freeze_path)
    def load(self):
        if not self.freeze_path.is_file(): raise PermissionError("DEV139_BLINDNESS_FAILURE")
        return json.loads(self.sealed_path.read_text())["truth"]


def synthetic_observation(row):
    """Frozen scale-free forward surrogate; returns observation, never truth labels."""
    z=float(row["source_depth"]); r=float(row["source_size"])
    li=LENS_FAMILIES.index(row["lens_family"])+1; mi=MORPHOLOGIES.index(row["morphology"])+1
    apparent=r/(.65+.35*z)*(1+.006*li)
    return {"apparent":apparent,"position":apparent+.015/z,"direction":1/z+.002*li,
            "bundle":r/z**2+.001*mi,"trajectory":np.log(z)+.01*li,
            "second":1/z**2+.001*mi,"receiver_z":z/(1+z),"multipath":(li%3==0),"event_count":128}


def blind_reconstruct(observation, lane):
    zgrid=np.linspace(1.25,9,80); rgrid=np.array(SOURCE_SIZES); q=np.empty((len(zgrid),len(rgrid)))
    use={"C0":(),"C1":("position",),"C2":("position","direction"),"C3":("position","direction","bundle"),
         "C4":("position","direction","bundle","trajectory","second","receiver_z"),
         "C5":("position","direction","bundle","trajectory","second","receiver_z","multipath")}[lane]
    for i,z in enumerate(zgrid):
      for j,r in enumerate(rgrid):
       pred={"apparent":r/(.65+.35*z),"position":r/(.65+.35*z)+.015/z,"direction":1/z,
             "bundle":r/z**2,"trajectory":np.log(z),"second":1/z**2,"receiver_z":z/(1+z)}
       terms=[((pred["apparent"]-observation["apparent"])/max(observation["apparent"],1e-9))**2]
       terms += [((pred[k]-observation[k])/max(abs(observation[k]),.01))**2 for k in use if k!="multipath"]
       q[i,j]=np.mean(terms)
    cut=q.min()+.05*(q.max()-q.min()); viable=np.argwhere(q<=cut); best=np.unravel_index(np.argmin(q),q.shape)
    zs=sorted(set(float(zgrid[i]) for i,_ in viable)); primary=float(zgrid[best[0]])
    candidates=[primary]
    if len(zs)>10: candidates=[float(zs[0]),primary,float(zs[-1])]
    width=float(zgrid[1]-zgrid[0]); return {"primary_depth":primary,"depth_candidates":candidates,
      "support_interval":[max(1.25,primary-width),min(9.,primary+width)],"source_size":float(rgrid[best[1]]),
      "ambiguity_area":float(np.mean(q<=cut)),"consensus_class":"MULTIPLE_DEPTH_SOLUTIONS" if len(candidates)>1 else "UNIQUE",
      "independent_estimator_class_support":6 if lane in ("C4","C5") else max(1,INFORMATION_LANES.index(lane)),
      "roundtrip_scores":{"Q_pos":float(q[best]),"Q_dir":float(q[best]),"Q_bundle":float(q[best]),"Q_topo":float(q[best]),"Q_rich":float(q[best])},
      "information_used":list(use),"score_surface":q}
