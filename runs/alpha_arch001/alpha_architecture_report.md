# PBUF ALPHA-ARCH-001 — Restoration of the V11 alpha architecture

## Outcome

The V11 hierarchy is restored as four physically distinct quantities:

- `alpha_QM`: Quantum Engine elastic amplitude.
- `alpha_resolved`: resolved elastic amplitude entering the cosmological pipeline.
- `alpha_T(a)`: scale-factor-dependent thermal LUT elastic amplitude.
- `alpha_EM`: electromagnetic fine-structure constant.

The exploratory `alpha_*` and the generic `alpha` propagated through CORE/FND/CONS/ERR do not have the physical role of any one of those quantities. They are the coefficient of an invented matter-state vertex. They are therefore renamed `g_dev` and classified as one new post-V11 construct. Its retained numerical premise `g_dev=1/137` does not identify it with `alpha_EM` and has no V11 derivation.

The frozen weak-lensing implementation was not changed.

## Authoritative relationships

The V11 flow is:

`Quantum Engine -> alpha_QM -> metadata/LUT -> {alpha_resolved, alpha_T(a)} -> cosmological pipeline`

V11 equation (4) implements `alpha_resolved ~= 3 alpha_EM`. This is a relationship between distinct quantities, not a license to substitute one symbol for the other. V11 describes the factor three as a motivating dimensional-counting consistency argument and defers a first-principles QFT derivation.

V11 equations (5) and (8) use `alpha_T(a)` in the thermal evolution. Equations (9) and (16) use `alpha_resolved` in flat-today closure and baryon normalization. The generic `alpha` in equation (15) is resolved by its surrounding cosmological context to `alpha_resolved`; it is not `alpha_EM` and not the exploratory matter vertex.

## Development audit

| Document | Finding | Corrected classification |
|---|---|---|
| CORE-001 | `alpha_*=1/137` loads the microscopic matter state and normalizes `s(rho)`. | `g_dev`, a post-V11 exploratory vertex coefficient; no V11 mapping. |
| MB-001 | Names `alpha_T` only as a possible closure input. | `alpha_T(a)`; V11 thermal LUT amplitude. No identification with the vertex. |
| FND-001 | No source, report, DOCX, or run directory is present. | Explicit coverage gap; no occurrence can be inferred. |
| FND-002 | Carries `alpha_*=1/137` as an underived assumption. | `g_dev`; post-V11 premise. |
| FND-003 | Uses the same coefficient in a three-component ontology. | `g_dev`; three-dimensionality does not map it to a V11 alpha. |
| FND-004 | Uses the coefficient for equal component loading and amplitude/power scaling. | `g_dev`; conditional post-V11 linear algebra. |
| FND-005 | Propagates the coefficient into proposed absolute observables; normalized ratios cancel it. | `g_dev`; no V11 or observational identification. |
| PHOTON-001 | Notes that the upstream coefficient is not linked to `n(u)` or `beta`. | `g_dev`; the negative mapping result is preserved. |
| CONS-001 | Replaced the fixed coefficient with generic `alpha`, effectively presenting it as framework-wide. | `g_dev`; its consistency result applies only to the exploratory vertex. |
| ERR-001 | Made the coefficient direct and formed `g_dev(1,1,1)`; its quadratic ratio gives 3. | `g_dev`; this is not V11 equation (4), which is amplitude-level. |

The occurrence-level equation and downstream mapping is in `traceability_matrix.csv`. The corrected symbol rules are in `corrected_notation_map.csv`.

## Dependency architecture

```text
regulators + field content
          |
          v
     Quantum Engine
          |
       alpha_QM
          |
    metadata / thermal LUT ---------------------+
          |                                      |
  alpha_resolved                           alpha_T(a)
     |          |                                |
 eq. (9)     eq. (16)                       eqs. (5)-(8)
     |          |                                |
 closure     Omega_b0                       thermal sigma

alpha_EM -- V11 eq. (4), implemented relation --> alpha_resolved

g_dev -- post-V11 matter vertex --> source --> [missing closure] --> u
   |                                                        |
   +--> equal-component ratios                    [missing n(u)/beta]

No edge identifies g_dev with alpha_QM, alpha_resolved, alpha_T(a), or alpha_EM.
```

The machine-readable edge list is `dependency_graph.csv`.

## Genuine deviations from V11

The genuine deviations are architectural, not merely typographical: a new matter-state coefficient was given an alpha name and the numerical value `1/137`; CONS-001 generalized it into a framework-wide `alpha`; the V11 production/resolution/thermal hierarchy was omitted; and ERR-001's quadratic component-counting factor could be confused with V11's amplitude-level equation (4). Each deviation and its resolution is recorded in `genuine_deviations.csv`.

## Applied corrections and recommendations

The active Python generators and regenerated CORE-001, FND-002 through FND-005, PHOTON-001, CONS-001, and ERR-001 outputs now use `g_dev`. MB-001 now writes `alpha_T(a)` explicitly. The obsolete generated `runs/cons001/alpha_dependency_graph.csv` was removed and replaced by `g_dev_dependency_graph.csv`. No weak-lensing source or archived weak-lensing output changed.

The DOCX mission specifications are retained as historical inputs. Their two ambiguous directives are explicitly mapped here: FND-002's `alpha_*` and every CONS-001 generic `alpha` mean the post-V11 `g_dev`. When those specifications are reissued, replace those symbols and add the sentence: “`g_dev` is not any V11 alpha quantity; `g_dev=1/137` is an independent exploratory premise.”

Do not create an FND-001 reconstruction from references in later documents. Supply the actual FND-001 artifact, then run the same occurrence audit and append its rows to the traceability matrix.

## Completion statement

Every alpha occurrence in every available required active exploratory artifact is either authoritative (`alpha_T(a)` in MB-001) or has been removed in favor of `g_dev`. Every alpha occurrence in the retained historical mission DOCX files is mapped by document and role. FND-001 is explicitly marked unavailable rather than guessed. Thus there is no unidentified alpha in the available corpus and no implicit physical identification based on numerical equality.
