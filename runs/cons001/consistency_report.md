# PBUF CONS-001 — Top-down consistency constraint on fundamental coupling

## Result: g_dev is currently indeterminate

Replacing the fixed coupling premise by symbolic `g_dev` after ERR-001 correction produces no finite bound and no preferred value. `g_dev` now directly normalizes the microscopic matter source, so the former inverse-rescaling argument is withdrawn. The common admissible region nevertheless remains unbounded because no completed sector supplies a value-selecting identity or independent cross-sector equation. Imposing the premise that the coupling is nonzero merely removes `g_dev=0` and does not bound its magnitude.

This is Outcome **D: existing theory is insufficient**. It is not evidence that a particular g_dev is wrong. No observation was fitted, no accepted numerical value was used in the analysis, and the frozen weak-lensing implementation was neither imported nor changed.

## Dependency graph

`g_dev -> microscopic source -> conditional coarse source -> u -> missing n(u) -> photon observables`

In parallel, `g_dev -> equal-component vertex -> bright/dark structure and normalized ratios`; g_dev cancels from the normalized ratios.

| Edge | From | To | Type | Established in | Consequence |
|---|---|---|---|---|---|
| C001-D01 | g_dev | microscopic source J | direct linear coupling | corrected CORE-001/FND-002/FND-003 | g_dev directly normalizes the matter-aligned microscopic source; no independent coupling multiplier remains. |
| C001-D03 | microscopic source J | coarse source s(rho) | conditional coarse graining | CORE-001; MB-001 | The map is conditional and its normalization is not independently fixed. |
| C001-D04 | coarse source s(rho) | continuum deformation u | conditional constitutive response | WL-003/MB-001 | Stiffness and closure coefficients introduce further independent scales. |
| C001-D05 | u | photon response n(u) | missing map | PHOTON-001 | beta=(dn/du)\|_0 is undetermined and is not identified with g_dev. |
| C001-D06 | n(u) | deflection and phase | conditional optical propagation | PHOTON-001 | Photon observables depend on n or beta and u, not on g_dev separately. |
| C001-D07 | g_dev | equal-component vertex vector g_dev*(1,1,1) | explicit premise | FND-004/FND-005 | Absolute amplitude scales with g_dev; normalized component-counting ratios cancel it. |
| C001-D08 | equal-component vertex | bright/dark and coherent ratios | linear algebra | FND-004/FND-005 | The ratios sqrt(3) and 3 are g_dev-independent for nonzero common coupling. |

## Constraint matrix

| Sector | G_dev dependence | Functions of g_dev | G_dev-independent outputs | Consistency condition | Admissible region | Strength | Interpretation |
|---|---|---|---|---|---|---|---|
| Foundational ontology (FND-001--FND-003) | directly normalizes the corrected microscopic matter vertex | microscopic source amplitude | three-component count and rotation-representation audit | a separate principle would be required to derive a numerical value | all finite g_dev; nonzero if the stipulated coupling is imposed | no magnitude bound | former inverse-rescaling conclusion withdrawn; numerical value remains a premise |
| CORE-001 microscopic model | linear source term -epsilon_* g_dev eta e.q | source and conditional equilibrium response amplitude | homogeneous stability, mode count, normalized coarse-graining form | stiffness positivity constrains kappa coefficients, not g_dev | all finite g_dev; nonzero for matter loading | excludes zero only if nonzero matter coupling is imposed | direct coupling is well defined but not numerically selected |
| FND-004 consequences | component vertices and unnormalized bright amplitude | \|g_vec\|=sqrt(3)\|g_dev\|; quadratic weight proportional to g_dev^2 | two dark modes and normalized multiplicity ratios | common equal nonzero coupling for bright/dark interpretation | g_dev != 0 for the stated coupling premise; no magnitude bound | excludes zero only if nonzero coupling is imposed | premise, not a derived consistency constraint |
| FND-005 experimental consequences | absolute calibrated source response | amplitudes proportional to g_dev; powers proportional to g_dev^2 | normalized coherent/single-channel ratios and component counting | access map and independent source calibration are absent | unbounded; separate g_dev unobservable | none on magnitude | proposed tests can constrain equality/counting before magnitude |
| Weak-lensing laboratory (WL-001) | none in frozen implementation | none | all archived fields, trajectories, residuals and RMSE | empirical scalar interface contains no g_dev mapping | all g_dev | none | cannot constrain a parameter absent from the code |
| Constitutive studies (WL-002/WL-002A/WL-003) | none in the implemented candidate scalar laws | only the upstream source amplitude once a physical closure is supplied | catalogue rankings and frozen reproducibility comparisons | physical micro--macro closure remains missing/conditional | all g_dev | none | empirical laws supply no equation relating their behavior to g_dev |
| Micro--macro closure (MB-001) | direct upstream source dependence should propagate through a completed closure | conditional coarse source | closure-gap finding and frozen Version-D reproduction | quantitative coarse-graining/response law absent | all finite g_dev; no closure-derived bound | none on magnitude | missing closure prevents a cross-sector consistency equation |
| Elasticity, rigidity, stiffness and thermal assumptions | no established equation ties g_dev to elastic, thermal, damping, or stiffness coefficients | none established | positivity/stability and symmetry conditions | K, G, kappa, thermal and response scales remain independent inputs | all g_dev | none | a relation would be new physics and is forbidden here |
| Photon coupling (PHOTON-001) | none established; optical beta or full n(u) is independent and missing | none without an added g_dev-to-electromagnetic map | conditional Fermat/ray equations and symmetry null tests | photon action/effective metric must independently fix n(u) | all g_dev | none | identifying beta with g_dev would be an unsupported assumption |

## Observable classification

| Observable | Dependence | Reason |
|---|---|---|
| normalized component count / two dark modes | independent of g_dev (assuming a nonzero equal access vertex) | depends on dimension and direction of the equal-coupling vector, not its magnitude |
| coherent-to-single and coherent-to-incoherent ratios | independent of g_dev | common powers of g_dev cancel |
| absolute microscopic/coarse response amplitude | microscopic source depends directly on g_dev; coarse response also depends on established stiffness/closure structure | the auxiliary rescaling degeneracy is gone, but the downstream closure is incomplete |
| continuum profile and weak-lensing residuals | independent in existing implementations | the frozen laboratory has no g_dev input |
| photon deflection and optical phase | conditional on beta and u; no established separate g_dev dependence | PHOTON-001 leaves n(u) unspecified |
| stability, stiffness positivity and symmetry gates | independent of g_dev in supplied equations | g_dev appears as a source scale rather than in the quadratic stability operator |

## Consistency overlap and dominant constraints

The intersection is `all finite g_dev; if the stipulated fundamental coupling must be nonzero, all finite g_dev except zero`. No completed sector supplies a closed interval, and none selects a preferred point. The strongest restrictions are not numerical bounds but identifiability gates:

1. No supplied symmetry or consistency identity selects a numerical g_dev.
2. The micro--macro response is not completed into an independent cross-sector constraint.
3. The photon response n(u), especially beta=(dn/du)|_0, is missing and is not linked to g_dev.

Consequently no sector dominates by numerical strength. The prior foundational/CORE rescaling obstruction no longer exists. Instead, the lack of a value-selecting theoretical identity, together with incomplete micro--macro and photon maps, prevents formation of independent simultaneous constraints. Apparent disagreement in a downstream absolute response could reveal missing closure/optical physics or a bad g_dev; the current framework cannot yet distinguish those explanations.

## Recommendation

Regard `g_dev` as **currently indeterminate**, not free in the sense of a completed predictive theory and not bounded or strongly constrained. The corrected microscopic source makes it a direct parameter rather than a factorization convention. A future consistency study can become informative only after a PBUF principle predicts or constrains its value, a quantitative micro--macro response propagates the source to `u`, and a photon/electromagnetic action fixes `n(u)` or `beta` and states whether it depends on `g_dev`. Those additions must be derived or independently specified before repeating the symbolic overlap; fitting them jointly to lensing is forbidden here.

Automated completion checks pass: **True**.
