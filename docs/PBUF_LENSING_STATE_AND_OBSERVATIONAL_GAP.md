# PBUF lensing state and observational-data gap

## Current result in one sentence

The native PBUF weak-lensing chain is frozen through a deterministic, finite-footprint sky representation for eight Abell 2744 source realizations, but it has **not** been compared with the intended Harvey–Massey (2024) pyRRG-JWST weak-lensing reconstruction because the released numerical map/catalogue with deterministic astrometry could not be recovered.

This is a morphology-only programme.  No PBUF output has been interpreted as physical convergence, mass, or an absolute weak-lensing amplitude.

## What is established on the native side

The native chain is closed in this order:

1. **Source ensemble (Dev171).** An independently acquired Abell 2744 spectroscopic catalogue supplies the frozen source membership and an eight-realization 3D phase-space/depth ensemble.  Spectroscopic redshift is not treated as direct geometric depth.
2. **Native propagation (Dev167).** The frozen distance-bound vector-pair law and propagation procedure provide the native disturbance evolution.  They were not changed after the source ensemble was frozen.
3. **Finite receipt (Dev168).** Finite received states retain the 3D position, direction, source lineage, flux/content information, and progression state needed by the downstream observer adapter.
4. **Observer output (Dev171).** The frozen observer channel bank produces eight 6×6 primary-channel arrays, stored as `observer_realization_00.npy` through `observer_realization_07.npy`.
5. **Coordinate provenance (Dev174).** The previously omitted serialized coordinate information was recovered without changing any 6×6 value.  Each realization now has receipt/source lineage, screen basis, dynamic extent, bin edges, native-bin footprints, and deterministic approximate sky footprints.

Dev175 independently rechecked the Dev174 coordinate-package manifest and the hashes of all eight frozen arrays.  They match the frozen package.  The resulting native coordinate status is:

| Item | Current status |
|---|---|
| Native excitation | established; not reopened |
| Pair law / propagation / receipt | frozen; not modified |
| Source ensemble | frozen; not modified |
| Observer channels and decoder | frozen; not retuned |
| 6×6 scientific arrays | byte-verified unchanged |
| Observer coordinate lineage | serialized / closed |
| Native grid → sky mapping | deterministic cell footprint serialized |
| Formal WCS for PBUF output | not claimed |

The coordinate footprints are finite cells, not fitted point registrations.  Their limited 6×6 resolution and the native discretization remain part of the result.

## Intended observation

The fixed external target is the Abell 2744 weak-lensing reconstruction in Harvey & Massey (2024), using pyRRG-JWST and UNCOVER DR1 data.  The paper reports shear measurements in F115W, F150W, and F200W and a 64×44 convergence map made with 12.8-arcsec pixels plus the documented smoothing.  It states that code, catalogues, and maps were publicly released through the pyRRG project.

The intended asset may be any one of the following, provided its provenance is unambiguous:

- the released convergence raster with sky coordinates;
- the exact released shear catalogue, including sky positions and shear quantities, from which only the paper's documented reconstruction could be reproduced; or
- an exact author-archived copy of either product.

The essential requirement is deterministic observational sample-to-sky mapping: a FITS WCS, explicit RA/DEC coordinates, an explicit sky grid, or a published reference origin with deterministic pixel geometry.  A manually registered image or a figure-derived map is not an acceptable substitute.

## What was searched and what is missing

Dev175 searched the current pyRRG repository, its reachable commit history and deleted files, all listed branches and pull-request refs, tags, GitHub releases, potential LFS history, repository-linked storage, the journal record/supplement route, Durham and EPFL deposits, arXiv ancillary material, and public author-controlled repositories.

One important finding is historical Abell 2744 training data in the pyRRG Git history.  Those files are legacy star/galaxy-classifier data and are not provenance-tied to the 2024 JWST analysis, its exact shear catalogue, or its convergence map.  They are explicitly classified as `UNVERIFIED_CANDIDATE` and were not used.

No qualifying 2024 product was recovered.  Consequently:

| External item | Status |
|---|---|
| Exact released convergence raster | unavailable |
| Exact released shear catalogue | unavailable |
| Author-archived numerical copy | unavailable |
| Immutable asset hash | unavailable because no asset was retrieved |
| Deterministic observational astrometry | insufficient |
| Observational uncertainty product | unavailable |
| Observational asset blocker | external unavailable |

The detailed query record and candidate classification are in `runs/dev175_pyrrg_recovery_blind_wl001/observational_asset_search_log.json` and `observational_asset_candidates.json`.

## Why no comparison was performed

The comparison is deliberately blind: both the PBUF prediction and the observational product must be frozen independently before either is compared.  With no valid observational product, there is no honest common sky grid and no footprint coverage calculation.  Opening the PBUF science arrays or computing a correlation in this state would defeat the gate rather than test the prediction.

The following were therefore not performed:

- projection of observations into native finite footprints;
- correlations for any of the eight realizations;
- centroid, peak, shape, radial, multipole, or topology comparisons;
- source-only and spatial-null controls; and
- any sign selection or morphology result.

Their status is `NOT_RUN_BLOCKED_BY_T09_T12`, not a failed morphology metric.

## Frozen comparison protocol if the authentic asset is recovered

Dev172 already predeclared the comparison and Dev175 preserved it unchanged.  Once a qualifying asset is recovered, it must first be hashed and frozen with its astrometry and processing.  Only then may the frozen PBUF arrays be accessed.

The common comparison grid must be derived only from Dev174 sky footprints and the observational astrometry.  Observational values should be area-integrated or averaged inside each native footprint; the 6×6 PBUF array must not be upsampled, interpolated, smoothed, or position-renormalized.  Coverage must be predeclared from observational geometry alone, never from PBUF values.

The frozen primary metric is zero-lag Pearson morphology correlation over the common mask after zero-mean/unit-RMS normalization.  The full required analysis is across all eight realizations, with ensemble mean, median, minimum, maximum, and spread, plus centroid, peak morphology, principal-axis shape, radial profile, low-order multipoles, threshold topology, a projected-source-only control, and predeclared spatial permutation/fixed-rotation nulls.  No best-realization promotion is allowed.

No translation, rotation, scale, mirror, centroid alignment, peak alignment, smoothing choice, sign choice, mask choice, source retuning, or native-physics change may be made to improve agreement.  Absolute amplitude remains off limits.

## Current outcome and responsible next step

The current result is **Outcome F**:

```text
OBSERVATIONAL_ASSET_BLOCKER=EXTERNAL_UNAVAILABLE
BLIND_WL_MORPHOLOGY_STATUS=NOT_EVALUATED
```

The native coordinate chain remains closed; the missing edge is external released data, not another theoretical or native-model issue.  The appropriate next action is to obtain an authentic author/repository archival copy of the 2024 numerical product, or to predeclare a separate independent weak-lensing dataset in a new development effort.  Neither route authorizes changing the frozen PBUF source, propagation, receipt, observer, footprints, registration, smoothing, mask, or metrics.

## Evidence locations

- Dev171 frozen native outputs: `runs/dev171_independent_3d_abell001/`
- Dev172 frozen comparison contract: `runs/dev172_blind_wl_morphology001/primary_metric_contract.json`
- Dev174 serialized footprint package: `runs/dev174_observer_coordinate_serialization001/`
- Dev175 recovery record and final contract: `runs/dev175_pyrrg_recovery_blind_wl001/`
- Canonical ledger and history: `docs/PBUF_DEVELOPMENT_LEDGER.*`, `docs/PBUF_HISTORICAL_ATTEMPT_INDEX.*`
