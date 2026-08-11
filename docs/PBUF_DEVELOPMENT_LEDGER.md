# PBUF Development Ledger

`LEDGER_EPOCH=POST_RESTORED_N6_CURRENT_STATE`
`PRE_LEDGER_HISTORY_CANONICAL=false`

GitHub code, this ledger, and current tests define canonical state. Conversation memory does not.

## PRE_LEDGER_HISTORY

`status=HISTORICAL_ONLY`

> Earlier development includes exploratory, superseded, and lower-dimensional surrogate implementations. Historical results must not override the current restored N6 code unless explicitly revalidated or preserved.

## LEDGER ENTRY 000 — CURRENT RESTORED N6 BASELINE

- Date: 2026-08-11
- Branch: `dev-doc-112-fullscale-vulkan-observer-validation`
- Start/final SHA: `f56e6a8bea77572ba80772e81806a68742c864d1`
- Repository at epoch capture: 52 untracked intended Dev156–165 files; tracking branch otherwise synchronized
- Tests: canonical Dev156–165 suite, 38 passed
- Active modules: native N6 F02/F03 dynamics and dispersion; source interaction; RAW bridge; 3D source family/static lens; wide-net mechanism audit; existing historical geometric observer
- Frozen: loaded scalar perturbation cancels to free F03; scalar static state does not derive geometry; no Dev165 candidate survives promotion
- Partial: H07, H14, H15; non-unique 3D source family
- Rejected: loaded scalar F03 coupling, scalar-derived geometry, scalar bond/memory-only redirection, minimal binary magnetic-like pair
- Underdetermined: richer relational state and loaded propagation mechanism families
- RAW: 116 RAW/FLT/FLC, one effective F814W channel, relative projected source constraint
- Observer: present; blocked by upstream finite-loaded receipt physics
- Current root blockers: native loaded propagation mechanism; absolute normalization; source-depth uniqueness

## LEDGER ENTRY 001 — DEV166 COMPLETE MISSING-PIECE AUDIT

- Date: 2026-08-11
- Start commit: `f56e6a8bea77572ba80772e81806a68742c864d1`
- Final commit: `7d919178b29e27832e231b740b9b6839f1805c13`
- Branch: `dev-doc-112-fullscale-vulkan-observer-validation`
- Push confirmed: `true`
- Question: what exists, what is missing, and which missing pieces are aliases or downstream symptoms?
- Result: complete dependency, semantic, routing/flux, variable-c, observer, circularity, and missing-piece audit; no new physics
- Status: `COMPLETE`
- Root blocker changed: yes—the immediate allocation requirement is retained, but its deeper identity is `MULTIPLE_CANDIDATES / UNDERDETERMINED`
- Downstream blocked: finite loaded propagation, received native state, finite-state observer reconnection
- Next allowed action: discussion/research review only
- `NEXT_DEV_AUTHORIZED=false`

## Mechanism ledger

| Mechanism | Status | Established by | Current evidence | Reopen only if | Dependents |
|---|---|---|---|---|---|
| Linear loaded F03 changes disturbance | REJECTED | Dev163 | background cancels; perturbation is free F03 | governing dynamic operator changes independently | loaded propagation |
| Scalar static state defines geometry | REJECTED | Dev164 | scalar differences are not lengths/orientations | independent native geometry semantics/state exists | geometric trajectory |
| Existing scalar bond state alone redirects | REJECTED | Dev165 H01 | no loaded redirection | governing loaded update changes independently | loaded response |
| Frozen F02/F03 memory alone redirects | REJECTED | Dev165 H06 | no loaded redirection | independently derived load-dependent memory law | loaded response |
| Minimal binary magnetic-like pair | REJECTED | Dev165 H12 | equilibrium/propagation gates fail | richer state and reversible law independently established | magnetic-like program |
| H07 directional allocation | PARTIAL | Dev165 | redirects, coefficient-free, conservative, permutation-covariant; not derived | independent derivation from native state | loaded response |
| H14 emergent geometry as output | PARTIAL | Dev165 | output effect only under nonderived routing | governing mechanism derived | observer path |
| H15 new primitive required | PARTIAL | Dev165 | missing capability localized, identity unknown | candidate semantics/law independently supplied | loaded response |
| State-dependent propagation speed | UNDERDETERMINED | Dev166 audit | conceptual lead only; `c_state` is not speed | native state-to-speed law independently derived | flux/allocation research |
| Richer magnetic-like interaction | UNDERDETERMINED | Dev165/166 | H12 does not exclude richer families | native semantics and reversible law supplied | loaded response |

## Missing-piece ledger

| Missing item | Layer | Status | Root/downstream | Blocked by | Blocks | Needed now |
|---|---|---|---|---|---|---|
| Loaded directional response | Native dynamics | MISSING | ROOT | deeper mechanism underdetermined | finite loaded propagation | yes |
| Finite loaded propagation | Native dynamics | BLOCKED | DOWNSTREAM | loaded directional response | received 3D state | yes |
| Received finite native 3D state | Receipt | BLOCKED | DOWNSTREAM | finite loaded propagation | observer reconnection | yes |
| Observer adapter | Interface | ADAPTER_REQUIRED | DOWNSTREAM | receipt contract absent | finite-state observer use | yes |
| Absolute normalization | Physical scale | UNRESOLVED | INDEPENDENT ROOT | — | physical comparison | no |
| Unique 3D depth | Source | NON_UNIQUE | INDEPENDENT ROOT | projected detector data | unique physical source | no |

## Permanent rules

Every future Dev begins by reading this ledger and current GitHub implementation, checking aliases and reopen conditions, and only then considering implementation. A Dev is incomplete until tests pass, ledger is updated, a commit SHA is recorded, push succeeds, and the remote SHA is verified. Rejected mechanisms are not retested unless their explicit reopen condition changes. Directional allocation is not assumed fundamental until state-dependent propagation and conservative flux equivalence are resolved.
