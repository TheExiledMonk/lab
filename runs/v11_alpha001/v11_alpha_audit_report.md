# PBUF V11-ALPHA-001 — Authoritative alpha audit

## Result

V11 first introduces `alpha_resolved` in section 2.2 as the single curvature parameter resolved from microphysics metadata. Section 2.3 identifies the upstream Quantum Engine output as `alpha_QM`. Equation (4), section 2.3.1, states `alpha_resolved ~= 3 alpha_EM = 3/137.036 ~= 0.0219` and defines `alpha_EM` as the electromagnetic fine-structure constant.

V11 does **not** present the factor three as a first-principles derivation. It says the relation is consistent with one contribution per spatial dimension and explicitly calls this a “motivating consistency argument”; a QFT derivation is deferred. The exact classification is therefore: equation (4) is a numerical observation/implemented relationship, while the spatial-dimension explanation is a motivating premise/interpretation.

## Complete alpha traceability

| Record | PDF page | Section | Equation | Statement | Classification | Assumptions | Downstream |
|---|---|---|---|---|---|---|---|
| V11-2.2-A1 | 2 | 2.2 | prose | thermal table supplies alpha_T(a) | definition | active thermal LUT | background, distance, and growth calculations |
| V11-2.2-A2 | 2 | 2.2 | prose | pipeline resolves alpha_resolved from microphysics metadata | definition | metadata with table-metadata fallback | single curvature/elastic parameter used by pipeline |
| V11-2.3-A1 | 2 | 2.3 | prose | Quantum Engine produces elastic amplitude alpha_QM | definition | specified regulators and field content | thermal LUT, then cosmological sector |
| V11-E04 | 2 | 2.3.1 | (4) | alpha_resolved ~= 3 alpha_EM = 3/137.036 ~= 0.0219 | empirical/numerical observation | current implementation; inherited quantum-microphysics metadata | fixed elastic amplitude and density normalization |
| V11-2.3.1-G1 | 3 | 2.3.1 | restatement of (4) | alpha_resolved ~= 3 alpha_EM; one contribution per spatial dimension | premise/motivating consistency argument | three-dimensional physical space and additive equal dimensional contributions | geometric interpretation only; first-principles QFT derivation deferred |
| V11-2.3.1-B1 | 3 | 2.3.1 | unnumbered / later (16) | Omega_b0 = 2 alpha_resolved | pipeline identity | two transverse electromagnetic polarizations; fixed normalization | present-day baryon density |
| V11-E05 | 3 | 2.3.2 | (5) | k_max(a)=epsilon_0,T(a)-alpha_T(a) | definition | thermal LUT fields | S(a), Omega_sigma_raw(a) |
| V11-E08 | 3 | 2.3.2 | (8) | Omega_sigma_raw(a)=alpha_T(a)(1-decay(a))S(a) | definition | activation and saturation definitions (6)-(7) | rescaling (10)-(11) |
| V11-E09 | 3 | 2.3.2 | (9) | Omega_sigma_target=1-Omega_m0-Omega_r0-alpha_resolved | pipeline normalization identity | flat_today normalization mode | sigma_rescale (10), Omega_sigma(a) (11) |
| V11-E15 | 4 | 2.3.3 | (15) | Omega_sigma(a)=alpha(1-exp(-a/R_max))S(a) | model definition | alpha denotes elastic amplitude; S defined by (7) | E(a), H(a), distances and growth via (12)-(14) |
| V11-E16 | 4 | 2.3.4 | (16) | Omega_b0=2 alpha_resolved | pipeline identity | fixed polarization-counting normalization | Omega_m0 through (17) |
| V11-T2-A1 | 7 | 5.1 Table 2 | Table 2 note | Omega_m0 and Omega_b0 derived from alpha_resolved; reported Omega_k0 corresponds to alpha_resolved | definition/reporting identity | V11 enforced identities (16)-(17) | reported PBUF parameters and likelihood predictions |
| V11-7-A1 | 11 | 7 | prose | alpha fully captures effective large-scale geometric response; no independent curvature parameter | model interpretation | effective-curvature usage in V11 pipeline | interpretation of background results |

Repeated prose references in sections 2.3.1–2.3.4 and Table 2 restate these same definitions and identities: `alpha_resolved` is fixed before likelihood evaluation, is not fitted to late-time data, determines `Omega_b0` through (16), enters flat-today closure through (9), and is reported as `Omega_k0` although it is not independent spatial curvature.

## V11-only dependency graph

`regulators + field content -> Quantum Engine -> alpha_QM, epsilon_0(T) -> metadata/LUT -> alpha_resolved, alpha_T(a)`

`alpha_EM --[numerical relation (4); dimensional-count motivation]--> alpha_resolved`

`alpha_T, epsilon_0,T -> (5)-(8) -> Omega_sigma_raw -> (9)-(11) -> Omega_sigma(a)`

`alpha_resolved -> (9) Omega_sigma_target; alpha_resolved -> (16) Omega_b0 -> (17) Omega_m0`

`Omega_sigma + matter/radiation -> (12)-(14) E(a), H(a) -> distances, q(a), growth, f sigma8`

## Geometric interpretation and ontology

“Geometric” in V11 means dimensional counting: one contribution for each of three spatial dimensions. V11 does not specify three microscopic state components, an SO(3) representation, or an equal-component vector whose norm produces the factor. Accordingly, the later three-dimensional ontology is compatible with the V11 motivation but is not its demonstrated origin.

## Later-development comparison and genuine deviations

| Document | Deviation from V11 | Assessment |
|---|---|---|
| CORE-001 | Recasts alpha_*=1/137 as a microscopic matter-state vertex and stipulates q in R^3; V11 instead uses alpha_resolved ~=3 alpha_EM ~=0.0219 as elastic/curvature metadata. | substantive symbol/value/role drift |
| FND-003 | Correctly labels the three-component spatial-vector bridge as requiring post-V11 representation premises; this qualification is stronger than V11's motivating dimensional-count argument. | clarification, not contradiction |
| ERR-001 | Its equal-component identities use g=alpha_*(1,1,1), giving a squared norm factor 3, but V11's factor 3 multiplies alpha_EM at amplitude level. ERR-001 does not establish these are the same construction. | unsupported retrospective identification if conflated |
| CORE/FND/ERR chain | Does not preserve the V11 Quantum Engine -> metadata/LUT -> cosmology chain, equations (5)-(17), or the distinction among alpha_QM, alpha_T(a), alpha_resolved, and generic model alpha. | missing V11 traceability |

## Recommendation

Choose **B) clarify wording**, **C) strengthen the mathematical derivation**, and **D) revise later development documents**. Preserve equation (4) as the V11 implemented numerical relationship, but do not promote its dimensional interpretation beyond V11's own “motivating consistency argument” classification. Later documents should restore the distinctions among `alpha_QM`, `alpha_T(a)`, `alpha_resolved`, generic `alpha`, and `alpha_EM`, and must not substitute the quadratic equal-component factor three for V11's amplitude-level relationship without a new derivation.

## Provenance and completion

Authoritative source: `/home/fabian/pbuf-test/docs/Planck-Bound_Unified_Framework_v11_preprint.pdf`  
SHA-256: `d4eece722f4329824eddd57850f89890dac1d88cad7d999319065d1e048585c2`

All 13 pages were read before later-development comparison. Completion checks pass: **True**.
