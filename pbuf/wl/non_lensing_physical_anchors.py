"""Frozen non-lensing anchor metadata, isolated from pre-anchor execution."""
from __future__ import annotations

def physical_anchor_manifest():
    return {"manifest":"DEV137_NON_LENSING_ANCHORS_V1","loaded_only_after_phase_G":True,
      "records":[
       {"id":"MEASURED_G","quantity":"Newtonian constant of gravitation","value":6.67430e-11,"units":"m^3 kg^-1 s^-2","uncertainty":1.5e-15,"source_provenance":"CODATA 2018 recommended value (NIST SP 959)","role":"VALIDATION"},
       {"id":"EARTH","quantity":"nominal equatorial radius","value":6.3781e6,"units":"m","uncertainty":0.0,"source_provenance":"IAU 2015 Resolution B3 nominal conversion constant","role":"VALIDATION"},
       {"id":"LABORATORY","quantity":"controlled source-response configuration","value":None,"units":None,"uncertainty":None,"source_provenance":"LAB_ANCHOR_DATA_UNAVAILABLE","role":"VALIDATION"},
       {"id":"SOLAR","quantity":"nominal solar radius","value":6.957e8,"units":"m","uncertainty":0.0,"source_provenance":"IAU 2015 Resolution B3 nominal conversion constant","role":"VALIDATION"}]}

class FrozenGlobalScale:
    def __init__(self): self.value=None; self.anchor=None; self.refits=0
    def calibrate(self,value,anchor):
      if self.value is not None:
        self.refits+=1; raise RuntimeError("GLOBAL_SCALE_IMMUTABLE_AFTER_FREEZE")
      self.value=float(value); self.anchor=anchor; return self.value
