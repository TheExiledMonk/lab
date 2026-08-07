# GEOMETRY-001 derived dependency graph

```mermaid
flowchart TD
    V11["V11: standard GR + Quantum Engine/alpha hierarchy"]
    MP["Missing constitutive principle:<br/>field content, dimensions, gauge, locality,<br/>normalization, G: chi -> g"]
    CHI["microscopic medium chi"]
    MET["one Lorentzian metric g_mu nu = G[chi]"]
    SM["standard matter action S_m[g,Psi]"]
    T["T_m^mu nu = 2/sqrt(-g) delta S_m/delta g_mu nu"]
    DG["metric response delta G/delta chi"]
    J["medium source J_A = -1/2 T : delta G/delta chi"]
    SS["missing covariant elastic action S_sigma"]
    CEQ["covariant chi equation"]
    CP["missing normalized C and projection P"]
    US["u = C[chi,g]; s = P[J]"]
    MB["conditional MB-001 balance<br/>K u - div(G_MB grad u) = s"]
    GEO["standard geodesics / null rays from g"]
    WL["frozen weak-lensing implementation<br/>(unchanged; no derived edge yet)"]

    V11 --> MP
    MP --> CHI
    CHI --> MET
    MET --> SM
    SM --> T
    MET --> DG
    T --> J
    DG --> J
    SS --> CEQ
    J --> CEQ
    CEQ --> CP
    CP --> US
    US --> MB
    MET --> GEO
    GEO -. "future validation only" .-> WL
```

Solid arrows are authoritative identities, necessary dependencies, or exact conditional derivations. The dashed arrow is deliberately not an implementation authorization. V11 provides the background elastic history but supplies no documented edge from an alpha quantity to `chi`, `G`, `S_sigma`, `C`, or `P`.
