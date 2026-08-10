"""Wide-net, contamination-safe native physical-scale candidate factory."""
from __future__ import annotations
from dataclasses import asdict, dataclass, field
from typing import Any

RESULT_CLASSES={"SCALE_ESTIMATE","SCALE_RELATION_ONLY","UNDERDETERMINED",
 "MISSING_PHYSICAL_INPUT","MISSING_NATIVE_INPUT","DIMENSIONALLY_INVALID",
 "NON_IDENTIFIABLE","CONFLICTING_INTERNAL_CONSTRAINTS","FORBIDDEN_INPUT_DEPENDENCE","NOT_APPLICABLE"}
FORBIDDEN=("rmax","strength_0p18","historical 0.18","planck length","kappa_target",
           "gamma_target","sawlens","sigma_crit","lcdm","class","camb","hst")

FAMILIES=(
"trajectory curvature","integrated direction change","fast/slow medium spatial-gradient restoration",
"six-neighbor Laplacian spacing restoration","physical source-loading / native cell response",
"physical mass-volume loading","source-surface response","far-field 1/r response",
"measured-G macroscopic response","potential-vs-response-gradient joint closure",
"field-gradient vs trajectory-curvature joint closure","bundle-Jacobian spatial evolution",
"bundle-Hessian spatial evolution","path-excess physicalization","propagation-speed / temporal closure",
"medium-relaxation closure","elastic-energy-density closure","elastic-gradient-energy closure",
"bounded-strain gradient closure","full dimensional-unit restoration",
"dimensional-rank / Buckingham-Pi identifiability","one-anchor universal-scale closure",
"Earth physical anchor","laboratory physical anchor","Solar-system physical anchor",
"cross-anchor universality","density-scaling universality","mass-scaling universality",
"radius-scaling universality","numerical-resolution invariance")
CLASSES=("TRAJECTORY_LOCAL","TRAJECTORY_INTEGRATED","MEDIUM_DYNAMIC","MEDIUM_STATIC",
"SOURCE_SIDE","SOURCE_SIDE","MEDIUM_STATIC","MEDIUM_STATIC","MACROSCOPIC_ANCHOR",
"MEDIUM_STATIC","MEDIUM_DYNAMIC","BUNDLE","BUNDLE","TRAJECTORY_INTEGRATED",
"MEDIUM_DYNAMIC","MEDIUM_DYNAMIC","ENERGY","ENERGY","ENERGY","DIMENSIONAL",
"DIMENSIONAL","MACROSCOPIC_ANCHOR","MACROSCOPIC_ANCHOR","MACROSCOPIC_ANCHOR",
"MACROSCOPIC_ANCHOR","MACROSCOPIC_ANCHOR","SOURCE_SIDE","SOURCE_SIDE","SOURCE_SIDE","DIMENSIONAL")

@dataclass
class NativeScaleEstimate:
    candidate_id:str; family:str; L0_m_per_native:float|None=None
    L0_uncertainty_m_per_native:float|None=None
    native_inputs:list=field(default_factory=list); physical_inputs:list=field(default_factory=list)
    empirical_anchors:list=field(default_factory=list); equations:list=field(default_factory=list)
    derivation_steps:list=field(default_factory=list); dimension_signature:str=""
    independent_of_target:bool=True; independent_of_lcdm:bool=True
    resolution_stability:str="NOT_TESTABLE"; mass_stability:str="NOT_TESTABLE"
    radius_stability:str="NOT_TESTABLE"; density_stability:str="NOT_TESTABLE"
    source_family_stability:str="NOT_TESTABLE"; internal_consistency:str="UNRESOLVED"
    cross_candidate_consistency:str="NOT_TESTABLE"; assumptions:list=field(default_factory=list)
    limitations:list=field(default_factory=list); status:str="UNDERDETERMINED"
    rejection_reason:str|None=None; independence_class:str="DIMENSIONAL"; L0_power:str="other"
    diagnostic_score:float=0.0
    def to_dict(self): return asdict(self)

def validate_candidate(c):
    if c.status not in RESULT_CLASSES: raise ValueError("invalid candidate result class")
    text=" ".join(map(str,asdict(c).values())).lower()
    hit=next((x for x in FORBIDDEN if x in text),None)
    if hit or not c.independent_of_target or not c.independent_of_lcdm:
        c.status="FORBIDDEN_INPUT_DEPENDENCE"; c.L0_m_per_native=None
        c.rejection_reason="FORBIDDEN_LCDM_DEPENDENCE" if hit in {"lcdm","class","camb"} or not c.independent_of_lcdm else "FORBIDDEN_INPUT_DEPENDENCE"
    return c

def candidate_registry(): return {f"S{i:02d}":f for i,f in enumerate(FAMILIES,1)}

def execute_internal_candidates(rank_audit):
    powers=(-1,0,-1,-2,3,3,0,-1,-2,-1,-1,-1,-1,1,"other","other",3,1,0,"other","other",0,0,0,0,0,0,0,0,0)
    results=[]
    reasons={
      1:("MISSING_PHYSICAL_INPUT","independent physical curvature relation absent"),
      2:("MISSING_PHYSICAL_INPUT","independent physical directional-evolution relation absent"),
      3:("SCALE_RELATION_ONLY","identifies U0/L0; response amplitude unit is free"),
      4:("SCALE_RELATION_ONLY","identifies K0_phys*U0/(S0*L0^2); stiffness and source units are free"),
      5:("NON_IDENTIFIABLE","physical-to-native source normalization requires the unknown cell volume"),
      6:("MISSING_NATIVE_INPUT","no independently normalized native mass-volume counterpart"),
      7:("MISSING_PHYSICAL_INPUT","surface response has no established physical bridge"),
      8:("MISSING_PHYSICAL_INPUT","far-field response amplitude unit is unnormalized"),
      9:("NON_IDENTIFIABLE","macroscopic acceleration cannot map to native response without U0 and T0"),
      10:("SCALE_RELATION_ONLY","amplitude-gradient closure identifies U0/L0 but not either factor"),
      11:("MISSING_PHYSICAL_INPUT","neither physical medium gradient nor physical curvature is independently normalized"),
      12:("SCALE_RELATION_ONLY","bundle Jacobian rate scales as L0^-1; physical rate absent"),
      13:("SCALE_RELATION_ONLY","primitive bundle Hessian contributes L0^-1; physical counterpart absent"),
      14:("MISSING_PHYSICAL_INPUT","native path excess exists but no independent physical path excess"),
      15:("UNDERDETERMINED","UNDERDETERMINED_L0_T0_PAIR"),
      16:("MISSING_PHYSICAL_INPUT","independent physical relaxation time and velocity absent"),
      17:("SCALE_RELATION_ONLY","energy normalization remains co-degenerate with L0^3"),
      18:("SCALE_RELATION_ONLY","gradient energy identifies K0_phys*U0^2*L0"),
      19:("SCALE_RELATION_ONLY","dimensionless strain leaves U0/L0 co-degeneracy"),
      20:("SCALE_RELATION_ONLY","complete unit restoration retains multiple free unit scales"),
      21:("NON_IDENTIFIABLE",rank_audit["L0_identifiability"]),
    }
    for i in range(1,22):
        status,reason=reasons[i]
        c=NativeScaleEstimate(f"S{i:02d}",FAMILIES[i-1],status=status,rejection_reason=reason,
          independence_class=CLASSES[i-1],L0_power=f"L0^{powers[i-1]}",
          equations=[reason],derivation_steps=["restore symbolic SI unit factors from primitive operation","test whether all non-L0 dimensions are independently fixed"],
          dimension_signature=f"L0^{powers[i-1]}",limitations=[reason])
        results.append(validate_candidate(c))
    return results

def dependency_graph(results):
    nodes=[{"id":c.candidate_id,"family":c.family,"independence_class":c.independence_class} for c in results]
    edges=[]
    for a in results:
      for b in results:
        if a.candidate_id < b.candidate_id and a.independence_class==b.independence_class:
          edges.append({"source":a.candidate_id,"target":b.candidate_id,"reason":"shared independence class"})
    return {"nodes":nodes,"edges":edges,"effective_independent_support":0}
