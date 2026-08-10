"""Real code dependencies only; conceptual relations are kept separate."""
def graph():
    edges=[("rho3","u_slow/u_fast initialization"),("u_slow/u_fast","c_state"),("c_state","bounded equilibrium source"),("bounded equilibrium source","u accumulated"),("u accumulated","epsilon"),("epsilon","sigma/W/K_tangent"),("packet controls","X"),("X_n","X_n+1"),("local frames","R_ij"),("R_ij","X_n+1"),("X","rho_X"),("rho_X","rbar"),("rbar history","packet trajectory"),("X1","lambda_n")]
    conceptual=[("c_state","L"),("u accumulated","excitation progression"),("packet trajectory","weak-lensing ray"),("rbar history","direction"),("direction","curvature"),("lambda_n","k_n")]
    return {"nodes":sorted({x for e in edges+conceptual for x in e}),"edges":[{"from":a,"to":b,"kind":"CODE_DEPENDENCY"} for a,b in edges],"conceptual_only":[{"from":a,"to":b,"kind":"CONCEPTUAL_RELATION_ONLY"} for a,b in conceptual]}
