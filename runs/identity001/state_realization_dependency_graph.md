# IDENTITY-001 state-realization dependency graph

```mermaid
flowchart TD
    I["Identity premise:<br/>medium = physical universe"]
    C["CONTINUUM-STATE-001:<br/>complete-state single-occupancy"]
    Q["Admissible physical-state space<br/>Q_phys / representation"]
    R["One complete realized state<br/>at each admissible state position"]
    N["No simultaneous incompatible<br/>complete configurations there"]
    AR["Anti-reification boundary:<br/>alternatives are not actual merely<br/>because they are represented"]
    H["Ordered history representation"]
    A["Missing relativistic<br/>actuality principle"]
    S["Single-actual-state ontology"]
    G["Growing-history ontology"]
    W["Whole-history ontology"]
    D["Missing dynamics:<br/>successors, probabilities,<br/>reversibility, arrow"]
    T["TIME-001 missing link:<br/>relational clock -> proper time"]
    V["V11 compatibility gate:<br/>Lorentz invariance, one metric,<br/>refoliation-equivalent observables"]

    I --> R
    C --> R
    R --> N
    I --> AR
    Q --> AR
    R --> H
    D --> H
    H --> S
    H --> G
    H --> W
    A --> S
    A --> G
    A --> W
    T --> V
    S -. "must pass" .-> V
    G -. "must pass" .-> V
    W -. "must pass" .-> V
```

Solid arrows indicate logical dependence or required input. Dashed arrows are compatibility tests, not derivations. Identity and continuum mechanics converge on same-position realization uniqueness. The graph branches only when actuality is assigned across an ordered history.

