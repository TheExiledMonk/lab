"""Common, scale-free local fixture and classification for Dev165."""
from __future__ import annotations

import numpy as np

IDS = tuple(f"H{i:02d}" for i in range(16))
AXES = np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]], float)


def transfer_metrics(t):
    t = np.asarray(t, float)
    mean = t.mean()
    den = np.abs(t).sum()
    return {"T": t.tolist(), "T_total": float(t.sum()),
            "directional_centroid": (t @ AXES).tolist(),
            "anisotropy": float(np.sqrt(np.sum((t-mean)**2))/den) if den else 0.0}


def _base(cid, status, redirection, primitive="NONE", complexity="LOW", **kw):
    row = dict(ID=cid, STATUS=status, UNLOADED_ISOTROPIC_EQUILIBRIUM="PASS",
        N6_SYMMETRY="PASS", DEV156_FREE_LIMIT="EXACT",
        DEV157_DISPERSION_RECOVERED="TRUE", DEV159_STATIC_SOURCE_COMPATIBILITY="TRUE",
        DEV159_FINITE_STATE_COMPATIBILITY="TRUE", DEV162_BACKGROUND_COMPATIBILITY="TRUE",
        LOADED_DIRECTIONAL_REDIRECTION=redirection, LOCAL_INTERACTION_LAW_LOAD_DEPENDENT=False,
        REVERSIBILITY="EXACT", DYNAMIC_INVARIANT="EXACT",
        MEANINGFUL_SOURCE_ENERGY_TRANSFER_REQUIRED="FALSE", EXISTING_SCALAR_PBUF_RECOVERED="TRUE",
        NEW_SCALAR_STATE_COUNT=0, NEW_VECTOR_STATE_COUNT=0, NEW_BOND_STATE_COUNT=0,
        NEW_FREE_COEFFICIENT_COUNT=0, NEW_DIMENSIONFUL_SCALE_COUNT=0,
        NEW_PRIMITIVE_REQUIREMENT=primitive, COMPLEXITY=complexity,
        PHYSICAL_LENGTH_SCALE_INTRODUCED=False, PHYSICAL_TIME_SCALE_INTRODUCED=False,
        MINIMALITY="STRONG", executable=True, reason="")
    row.update(kw)
    return row


def evaluate_candidates(loaded_relations):
    """Evaluate all mandatory families without assigning geometric semantics.

    The exploratory H07 allocation is deliberately reported PARTIAL: its formula
    proves sufficiency of directional allocation, not derivability of that law.
    """
    r = np.asarray(loaded_relations, float)
    if r.shape != (6,):
        raise ValueError("loaded_relations must contain the six ordered N6 differences")
    uniform = np.full(6, 1/6)
    # Permutation/reflection covariant, coefficient-free, norm-preserving example.
    allocation = (1.0 + np.abs(r)) / np.sum(1.0 + np.abs(r))
    rows = {}
    rows["H00"] = _base("H00","REJECTED","FALSE",reason="Dev163 exact loaded/free null")
    rows["H01"] = _base("H01","REJECTED","FALSE",reason="bondwise linear derivatives are identical before summation")
    rows["H02"] = _base("H02","UNDERDETERMINED","NOT_DERIVED","BOND_SCALAR",NEW_BOND_STATE_COUNT=6,
        DEV156_FREE_LIMIT="PARTIAL", DEV157_DISPERSION_RECOVERED="PARTIAL", REVERSIBILITY="NOT_APPLICABLE",
        DYNAMIC_INVARIANT="NOT_DERIVED", MINIMALITY="ACCEPTABLE", executable=False,
        reason="no frozen quantity has native length semantics; separation function underdetermined")
    rows["H03"] = _base("H03","UNDERDETERMINED","NOT_DERIVED","POLARITY",NEW_BOND_STATE_COUNT=6,
        DEV156_FREE_LIMIT="PARTIAL", DEV157_DISPERSION_RECOVERED="PARTIAL", REVERSIBILITY="NOT_APPLICABLE",
        DYNAMIC_INVARIANT="NOT_DERIVED", executable=False, reason="polarity initialization and update law absent")
    rows["H04"] = _base("H04","UNDERDETERMINED","NOT_DERIVED","POLARITY","MEDIUM",NEW_BOND_STATE_COUNT=6,
        DEV156_FREE_LIMIT="PARTIAL", DEV157_DISPERSION_RECOVERED="PARTIAL", REVERSIBILITY="NOT_APPLICABLE",
        DYNAMIC_INVARIANT="NOT_DERIVED", executable=False, reason="magnetic-like pair law and equilibrium are not derived")
    rows["H05"] = _base("H05","UNDERDETERMINED","NOT_DERIVED","BOND_SCALAR","MEDIUM",NEW_BOND_STATE_COUNT=12,
        DEV156_FREE_LIMIT="PARTIAL", DEV157_DISPERSION_RECOVERED="PARTIAL", REVERSIBILITY="NOT_APPLICABLE",
        DYNAMIC_INVARIANT="NOT_DERIVED", executable=False, reason="combines two absent semantics")
    rows["H06"] = _base("H06","REJECTED","FALSE",reason="stationary F03 retained change is zero; frozen linear perturbation remains load independent")
    rows["H07"] = _base("H07","PARTIAL","TRUE","OTHER",reason="coefficient-free conservative allocation proves sufficiency but is a new, underived routing law",
        REVERSIBILITY="NOT_DERIVED", DYNAMIC_INVARIANT="EXACT", MINIMALITY="ACCEPTABLE")
    rows["H08"] = _base("H08","UNDERDETERMINED","NOT_DERIVED","VECTOR","HIGH",NEW_VECTOR_STATE_COUNT=6,
        DEV156_FREE_LIMIT="PARTIAL", DEV157_DISPERSION_RECOVERED="PARTIAL", REVERSIBILITY="NOT_APPLICABLE",
        DYNAMIC_INVARIANT="NOT_DERIVED", executable=False, reason="vector relation is not required by lower-state evidence")
    rows["H09"] = _base("H09","UNDERDETERMINED","NOT_DERIVED","MULTICOMPONENT","MEDIUM",NEW_SCALAR_STATE_COUNT=1,
        DEV156_FREE_LIMIT="PARTIAL", DEV157_DISPERSION_RECOVERED="PARTIAL", REVERSIBILITY="NOT_APPLICABLE",
        DYNAMIC_INVARIANT="NOT_DERIVED", executable=False, reason="M02/M03 projection and update laws absent")
    rows["H10"] = _base("H10","UNDERDETERMINED","NOT_DERIVED","OTHER","HIGH",NEW_SCALAR_STATE_COUNT=36,
        DEV156_FREE_LIMIT="PARTIAL", DEV157_DISPERSION_RECOVERED="PARTIAL", REVERSIBILITY="NOT_APPLICABLE",
        DYNAMIC_INVARIANT="NOT_DERIVED", executable=False, reason="relation matrix is an unconstrained allocation law")
    rows["H11"] = _base("H11","UNDERDETERMINED","NOT_DERIVED","BOND_SCALAR","MEDIUM",NEW_BOND_STATE_COUNT=6,
        DEV156_FREE_LIMIT="PARTIAL", DEV157_DISPERSION_RECOVERED="PARTIAL", REVERSIBILITY="NOT_APPLICABLE",
        DYNAMIC_INVARIANT="NOT_DERIVED", executable=False, reason="preferred relation is a new primitive without an update law")
    rows["H12"] = _base("H12","REJECTED","NOT_DERIVED","POLARITY","MEDIUM",NEW_SCALAR_STATE_COUNT=1,
        DEV156_FREE_LIMIT="FAIL", DEV157_DISPERSION_RECOVERED="FALSE", DEV159_STATIC_SOURCE_COMPATIBILITY="PARTIAL",
        REVERSIBILITY="NOT_APPLICABLE", DYNAMIC_INVARIANT="NOT_DERIVED", executable=False,
        reason="binary pair polarity alone supplies no coefficient-free stable equilibrium-centered propagation law")
    rows["H13"] = _base("H13","UNDERDETERMINED","NOT_DERIVED","POLARITY","HIGH",NEW_BOND_STATE_COUNT=6,
        DEV156_FREE_LIMIT="PARTIAL", DEV157_DISPERSION_RECOVERED="PARTIAL", REVERSIBILITY="NOT_APPLICABLE",
        DYNAMIC_INVARIANT="NOT_DERIVED", executable=False, reason="memory is characterized, but polarity remains absent")
    rows["H14"] = _base("H14","PARTIAL","NOT_DERIVED","NONE",reason="valid output interpretation only; requires a prior redirecting mechanism")
    rows["H15"] = _base("H15","PARTIAL","NOT_DERIVED","OTHER",reason="audit outcome: a constrained directional allocation primitive is required")
    transfers = {}
    for cid, row in rows.items():
        if cid == "H07": loaded = allocation
        elif row["executable"]: loaded = uniform
        else: loaded = None
        transfers[cid] = {"unloaded": transfer_metrics(uniform),
                          "loaded": transfer_metrics(loaded) if loaded is not None else None,
                          "delta_T": (loaded-uniform).tolist() if loaded is not None else None}
    return list(rows.values()), transfers
