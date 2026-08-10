"""Exact symbolic unit restoration for the frozen PBUF weak-lensing path.

Dimensions are integer exponents of (M, L, T).  The module intentionally
does not assign numerical values to native unit scales.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

BASE_DIMENSIONS = ("M", "L", "T")


def _rref(rows):
    a = [[Fraction(x) for x in row] for row in rows]
    if not a:
        return a, []
    pivots, r = [], 0
    for c in range(len(a[0])):
        p = next((i for i in range(r, len(a)) if a[i][c]), None)
        if p is None:
            continue
        a[r], a[p] = a[p], a[r]
        q = a[r][c]; a[r] = [x/q for x in a[r]]
        for i in range(len(a)):
            if i != r and a[i][c]:
                q = a[i][c]; a[i] = [x-q*y for x, y in zip(a[i], a[r])]
        pivots.append(c); r += 1
        if r == len(a): break
    return a, pivots


def null_space(matrix, ncols=None):
    """Return an exact rational basis for the right null space."""
    n = ncols if ncols is not None else (len(matrix[0]) if matrix else 0)
    rr, pivots = _rref(matrix)
    free = [j for j in range(n) if j not in pivots]
    out = []
    for f in free:
        v = [Fraction(0) for _ in range(n)]; v[f] = Fraction(1)
        for i, p in enumerate(pivots): v[p] = -rr[i][f]
        out.append(v)
    return out


@dataclass(frozen=True)
class DimensionalSystem:
    scales: tuple[str, ...]
    equations: tuple[tuple[int, ...], ...]

    def audit(self, length_scale="L0"):
        basis = null_space(self.equations, len(self.scales))
        rank = len(self.scales)-len(basis)
        li = self.scales.index(length_scale)
        lfree = any(v[li] for v in basis)
        if not lfree:
            ident = "UNIQUELY_IDENTIFIABLE_INTERNALLY"
        elif len(basis) == 1 and all(j == li or not basis[0][j] for j in range(len(basis[0]))):
            ident = "IDENTIFIABLE_WITH_ONE_EXTERNAL_ANCHOR"
        elif len(basis) == 1:
            ident = "CO_DEGENERATE_WITH_OTHER_UNIT_SCALE"
        elif lfree:
            ident = "NON_IDENTIFIABLE_FROM_CURRENT_PHYSICS"
        else:
            ident = "IDENTIFIABLE_WITH_TWO_OR_MORE_ANCHORS"
        combos=[]
        for v in basis:
            terms=[f"{s}^{v[i]}" for i,s in enumerate(self.scales) if v[i]]
            combos.append(" * ".join(terms) or "1")
        return {"unknown_unit_scales_count":len(self.scales),
                "independent_dimensional_equations_count":rank,"matrix_rank":rank,
                "nullity":len(basis),"L0_identifiability":ident,
                "remaining_free_combinations":combos,
                "nullspace_basis":[[str(x) for x in v] for v in basis]}


def frozen_medium_unit_contract():
    """Symbolic contract traced to current primitive operations."""
    return {
      "base_dimensions": list(BASE_DIMENSIONS),
      "unit_scales":{"coordinate":"L0","mass":"M0","time":"T0",
                     "response":"U0","stiffness":"K0_phys","source":"S0"},
      "primitives":{
        "c_state":{"dimensions":"U0","source":"pbuf.wl.native_response:construct_c_state"},
        "fast_channel":{"dimensions":"U0","source":"pbuf.wl.backends.cpu:pair transfer"},
        "slow_channel":{"dimensions":"U0","source":"pbuf.wl.backends.cpu:pair transfer"},
        "displacement_response":{"dimensions":"U0","source":"native response array"},
        "strain":{"dimensions":"U0/L0","source":"neighbor difference / native spacing"},
        "bounded_strain_stress":{"dimensions":"K0_phys*U0/L0","source":"bounded-strain constitutive law"},
        "source_term":{"dimensions":"S0","source":"native source loading"},
        "A_ij":{"dimensions":"U0","source":"0.03*delta_u_fast + 0.003*delta_u_slow"},
        "trajectory_coordinate":{"dimensions":"L0","source":"pbuf.wl.backends.cpu:positions"},
        "path_length":{"dimensions":"L0","source":"Euclidean segment sum"},
        "direction":{"dimensions":"1","source":"normalized displacement"},
        "curvature":{"dimensions":"L0^-1","source":"delta(direction)/delta(path)"},
        "bundle_Jacobian":{"dimensions":"1","source":"endpoint derivative wrt launch coordinate"},
        "bundle_Hessian":{"dimensions":"L0^-1","source":"second endpoint derivative"}},
      "frozen_coefficients":{"fast":0.03,"slow":0.003,"physical_dimensions":"symbolic transfer coefficient"}}


def default_dimensional_system():
    # Current implementation fixes ratios/operations but supplies no SI equation.
    # Columns L0,U0,K0_phys,S0,T0,M0; equations encode only K*U/L ~ S.
    return DimensionalSystem(("L0","U0","K0_phys","S0","T0","M0"),
                             ((-1,1,1,-1,0,0),))


def restored_operators():
    return [
      {"operator":"gradient","native_equation":"Delta u / Delta x_n","physical_equation":"(U0/L0) Delta_n u","L0_power":-1,"unknowns":["U0"]},
      {"operator":"six_neighbor_laplacian","native_equation":"sum_neighbors(u_i-u_j)","physical_equation":"(U0/L0^2) Laplacian_n u","L0_power":-2,"unknowns":["U0","K0_phys","S0"]},
      {"operator":"cell_volume","native_equation":"dV_n","physical_equation":"L0^3 dV_n","L0_power":3,"unknowns":[]},
      {"operator":"trajectory_curvature","native_equation":"dtheta/ds_n","physical_equation":"kappa_n/L0","L0_power":-1,"unknowns":[]},
      {"operator":"path_excess","native_equation":"L_path,n-L_straight,n","physical_equation":"L0 DeltaL_n","L0_power":1,"unknowns":[]},
      {"operator":"bundle_rate","native_equation":"dJ/ds_n","physical_equation":"(1/L0)dJ/ds_n","L0_power":-1,"unknowns":[]},
      {"operator":"gradient_energy_3d","native_equation":"sum |grad_n u|^2 dV_n","physical_equation":"U0^2 L0 sum |grad_n u|^2 dV_n","L0_power":1,"unknowns":["U0","K0_phys"]}]


def spatial_only_native_basis():
    """Dev141 dimensional system: time is deliberately not a native basis axis."""
    system = DimensionalSystem(("L0", "U0", "K0_phys", "S0"),
                               ((-1, 1, 1, -1),))
    out = system.audit()
    out.update({"base_dimensions": ["M", "L"], "native_time_dimension": False,
                "unknowns": list(system.scales)})
    return out


def emergent_time_mapping(length_scale_m_per_native, path_native, *, c=299_792_458.0):
    """Convert an established physical path to elapsed seconds; accept no T0."""
    if length_scale_m_per_native is None:
        return {"available": False, "EMERGENT_TIME_PER_NATIVE_LENGTH": None,
                "elapsed_seconds": None, "relation": "s_physical / c"}
    if length_scale_m_per_native <= 0 or c <= 0:
        raise ValueError("length scale and c must be positive")
    e0 = float(length_scale_m_per_native / c)
    return {"available": True, "EMERGENT_TIME_PER_NATIVE_LENGTH": e0,
            "elapsed_seconds": float(path_native * e0), "relation": "L0*s_native/c",
            "native_T0": False}


def dev137_time_ontology_reconciliation():
    """Compare the frozen Dev137 model with the corrected Dev141 ontology."""
    old = default_dimensional_system().audit()
    corrected = spatial_only_native_basis()
    return {"old_unknown_count": old["unknown_unit_scales_count"],
            "corrected_unknown_count": corrected["unknown_unit_scales_count"],
            "old_rank": old["matrix_rank"], "corrected_rank": corrected["matrix_rank"],
            "old_nullity": old["nullity"], "corrected_nullity": corrected["nullity"],
            "removed_equations": [],
            "removed_unknowns": ["T0", "M0"],
            "T0_degeneracy_classification": "ONTOLOGY_ARTIFACT_REMOVED",
            "remaining_degeneracies": corrected["remaining_free_combinations"],
            "apparent_old_degeneracy_ontological": True}


def reject_fundamental_time_dimension(*, T0=None, solver_iterations=None):
    if T0 is not None or solver_iterations is not None:
        raise ValueError("REJECT_FUNDAMENTAL_TIME_DIMENSION")
    return "REJECT_FUNDAMENTAL_TIME_DIMENSION"
