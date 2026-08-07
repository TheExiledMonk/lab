# STATE-001 ontology-to-emergent-time dependency graph

```mermaid
flowchart TD
    V["V11 authority:<br/>effective GR, Lorentz invariance,<br/>standard clocks and dynamics"]
    O["PBUF premise:<br/>fundamental 3D medium"]
    M1["Missing complete-state specification:<br/>fields, matter, constraints, gauge, boundaries"]
    S["Candidate physical state<br/>X in Q_phys"]
    M2["Missing global-state conditions:<br/>Cauchy surface / admissible slicing"]
    GS["Conditional global instantaneous state"]
    T["TIME-001 conditional ordered history<br/>gamma modulo reparametrization"]
    M3["TIME-001 missing temporal identification:<br/>clock functional + metric proper time"]
    ET["Effective relativistic temporal order/duration"]
    M4["STATE-001 missing actuality-and-slicing principle"]
    A1["Single-actual-state ontology"]
    A2["Growing-history ontology"]
    A3["Whole-history ontology"]
    D["Missing dynamics / boundary conditions"]
    DU["Determinism, branching, reversal,<br/>arrow and reconstructibility"]
    R["Coordinate/spacetime histories<br/>(representations, not actuality rules)"]

    O --> S
    M1 --> S
    S --> GS
    M2 --> GS
    GS --> T
    T --> R
    T --> ET
    M3 --> ET
    M4 --> A1
    M4 --> A2
    M4 --> A3
    T --> A1
    T --> A2
    T --> A3
    D --> DU
    T --> DU
    ET --> V
    A1 -. "must pass covariance test" .-> V
    A2 -. "must pass covariance test" .-> V
    A3 -. "must pass covariance test" .-> V
```

The three actuality branches are alternatives, not simultaneous conclusions. The graph separates actuality from dynamics: deterministic evolution would not by itself select any actuality branch. Dashed edges are required compatibility tests with V11, not derived equivalences.

