# PBUF ERR-001 — Removal of the non-PBUF auxiliary coupling

## Corrective result

The independent exploratory parameter conventionally written as lambda has been removed from CORE-001 and every affected downstream derivation. It was not the cosmological constant and is not part of PBUF. No replacement coefficient, hidden normalization factor, or new free fundamental parameter was introduced.

The corrected matter vertex contains `g_dev` directly. The former inverse-rescaling degeneracy and every conclusion relying on it are withdrawn. All seven milestones were regenerated and their validations pass; the frozen weak-lensing benchmark was not imported, executed, or modified by this correction.

## Milestone audit

| Milestone | Affected locations before correction | Correction | Status |
|---|---|---|---|
| CORE-001 | microscopic energy, continuum source, traceability and model mapping | Removed the auxiliary multiplier; matter couples directly through g_dev. | corrected and revalidated |
| FND-002 | A02 assumption audit and irreducible-postulate discussion | Withdrew the effective-product substitution and source-rescaling degeneracy; retained g_dev as a direct, underived PBUF premise. | corrected and revalidated |
| FND-003 | T08/T09, coupling derivation, postulate P5 and recommendation | Replaced the product-identifiability proof by direct-source dependence; explicitly withdrew the former inverse-rescaling conclusion. | corrected and revalidated |
| FND-004 | P10 identifiability consequence | Absolute calibrated response is now g_dev-sensitive; normalized component ratios still cancel g_dev. | corrected and revalidated |
| FND-005 | P08, ontology comparison and generated catalogue | Removed the nuisance-normalization claim; distinguished absolute vertex sensitivity from g_dev-independent ratios. | corrected and revalidated |
| PHOTON-001 | no auxiliary coupling occurrence | No coupling equation required modification; the missing optical response remains independent of g_dev unless PBUF derives a link. | unchanged and revalidated |
| CONS-001 | dependency graph, three sector rows, overlap logic, report and recommendation | Reran the audit with direct g_dev loading and withdrew rescaling degeneracy as a reason for indeterminacy. | revised and revalidated |

The DOCX specifications for all seven milestones were also inspected. They contained no instance of the auxiliary product and required no binary edits.

## Corrected equations

| ID | Milestone | Corrected equation | Effect |
|---|---|---|---|
| ERR-001-E01 | CORE-001 | F=epsilon_* sum_i[kappa_0\|q_i\|^2/2+kappa_1 sum_<ij>\|q_j-q_i\|^2/2-g_dev eta_i e.q_i] | g_dev is the sole fundamental matter-state coupling. |
| ERR-001-E02 | CORE-001 | s(rho)=epsilon_* g_dev (rho/rho_*)/a^d | The coarse source follows directly from the corrected vertex. |
| ERR-001-E03 | CORE-001/MB-001 | K u-div(G grad u)=s(rho) | Form unchanged; corrected source E02 replaces the former exploratory source. |
| ERR-001-E04 | FND-004/FND-005 | g_vec=g_dev(1,1,1); \|g_vec\|=sqrt(3)\|g_dev\| | Absolute vertex is g_dev-sensitive; two dark directions are unchanged. |
| ERR-001-E05 | FND-004/FND-005 | \|sum_i g_dev\|^2 / sum_i\|g_dev\|^2=3 | The normalized coherent ratio remains g_dev-independent. |
| ERR-001-E06 | PHOTON-001 | n(u)=1+beta u+O(u^2) | Unchanged: beta remains a missing optical response and is not assumed to equal g_dev or a new free fundamental coupling. |

## Change log

| Item | Before ERR-001 | After ERR-001 | Classification |
|---|---|---|---|
| microscopic source normalization | g_dev multiplied by an auxiliary exploratory coupling | g_dev appears directly | revised |
| coarse source | proportional to a two-factor coupling product | proportional directly to g_dev | revised |
| inverse-rescaling degeneracy | claimed to make g_dev separately unidentifiable | withdrawn as an artefact of the auxiliary parameter | withdrawn |
| numerical derivation of g_dev | not derived | still not derived by any supplied symmetry or consistency identity | unchanged |
| component multiplicities and normalized ratios | independent of coupling magnitude | independent of coupling magnitude | unchanged |
| CORE stability and propagation length | controlled by stiffness coefficients | controlled by stiffness coefficients | unchanged |
| micro--macro closure status | incomplete | incomplete, but no longer obscured by source-rescaling degeneracy | unchanged boundary |
| photon response | n(u) and beta missing | n(u) and beta remain missing; no g_dev link is invented | unchanged |
| CONS-001 final classification | indeterminate partly because of inverse rescaling | indeterminate solely because no value-selecting or closed cross-sector constraints exist | revised rationale |

## Revised dependency graph

`g_dev -> microscopic source -> conditional coarse source -> u -> missing n(u) -> optical observables`

The parallel component branch is `g_dev -> equal-component vertex -> bright/dark structure and normalized ratios`, with the magnitude canceling only in normalized ratios.

## Revised CONS-001 conclusion

`g_dev` remains currently indeterminate, but for revised reasons. The previous rescaling argument was wholly an artefact of the non-PBUF auxiliary parameter and has been withdrawn. Direct microscopic dependence now makes g_dev operationally meaningful in a completed calibrated model. Nevertheless, the present theory provides no value-selecting symmetry or identity, no closed micro--macro consistency equation, and no photon-response relation tied to g_dev. Thus there is still no finite interval or preferred value from top-down consistency alone.

Required future theory is limited to genuine PBUF derivations: a principle that derives or bounds g_dev, a quantitative closure without a new fundamental coupling, and a photon/electromagnetic action stating whether its response follows from g_dev.

## Absence confirmation

No independent lambda-like coupling or equivalent surrogate remains in the corrected CORE-001 through CONS-001 sources or regenerated artifacts. Ordinary Python anonymous-function syntax is a programming-language construct, not a physical parameter. PHOTON-001 uses no such physical coupling; its ray functional is expressed without introducing one.

Automated ERR-001 checks pass: **True**.
