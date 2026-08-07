# TIME-001 medium-evolution dependency graph

```mermaid
flowchart TD
    O["PBUF ontological premise:<br/>3D medium"]
    S["Complete instantaneous Cauchy state<br/>X in Q_phys = configurations / gauge"]
    M1["Missing state specification:<br/>fields, momenta/memory, gauge quotient"]
    M2["Missing admissible-history law:<br/>orientation, causality, well-posedness"]
    H["Oriented history class<br/>gamma: I -> Q_phys modulo reparametrization"]
    L["Arbitrary sequence label lambda<br/>(not physical time)"]
    M3["Missing relational clock functional<br/>and monotonicity domain"]
    C["Relational clock readings"]
    M4["GEOMETRY-001 missing map:<br/>medium clocks/rulers -> one Lorentzian metric"]
    G["Effective Lorentzian metric g"]
    P["Effective proper time and causal cones"]
    V["V11 relativistic regime:<br/>GR, Lorentz invariance, H(a), propagation"]
    CG["Optional coarse graining / causal limits / dissipation"]
    NR["Non-unique or operationally inaccessible past"]

    O --> S
    M1 --> S
    S --> H
    M2 --> H
    H --> L
    H --> C
    M3 --> C
    C --> G
    M4 --> G
    G --> P
    P --> V
    H --> CG
    CG --> NR
```

Solid arrows show required dependencies, not derived field equations. `lambda` is a representation of the history and has no edge to proper time. Past non-reconstructibility is conditional on the additional branch through coarse graining, causal limits, or non-invertibility; it does not follow directly from emergent ordering.

