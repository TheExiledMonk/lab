# Cluster baryonic source input inventory 001

## Purpose

Step back from weak-lensing prediction and inspect the physical baryonic source chain one link at a time for the five current Frontier Fields clusters.

This is a **fact-finding inventory**, not a lensing run, not an SI-density reconstruction, and not a calibration fit.

The preceding baryonic-density normalization audit established that the current F160W source is photometrically calibrated but remains a morphology proxy rather than an absolute baryonic mass-density field. This lab separates the unresolved chain into explicit links so each can be closed independently.

## Directly admitted observed input

The STScI Hubble Frontier Fields archive directly lists the cluster redshifts used here:

- Abell 2744: `z = 0.308`
- MACSJ0416.1-2403: `z = 0.396`
- MACSJ1149.5+2223: `z = 0.543`
- Abell S1063: `z = 0.348`
- Abell 370: `z = 0.375`

Provenance: `https://archive.stsci.edu/prepds/frontier/`

Redshift is treated only as observed source metadata. The lab does **not** convert redshift to luminosity distance or angular-diameter distance. Such a conversion requires a cosmological geometry and therefore must be independently supplied or explicitly audited later rather than silently importing LCDM into the PBUF source normalization.

## Audited links

For each cluster the lab inventories:

1. observed cluster redshift;
2. HST detector-to-flux calibration;
3. redshift-to-physical-distance geometry;
4. stellar baryonic mass conversion/catalog;
5. hot/diffuse intracluster gas baryons;
6. physical pixel area;
7. surface-density to volume-density deprojection;
8. preservation of absolute amplitude into the native `rho2/rho3` source.

Each link is classified as one of:

```text
OBSERVED_AVAILABLE
PIPELINE_AVAILABLE_BUT_AMPLITUDE_ERASED
ASTROPHYSICAL_CONVERSION_REQUIRED
COSMOLOGY_DEPENDENT_GEOMETRY
INDEPENDENT_EXTERNAL_DATA_REQUIRED
NOT_YET_CLOSED
```

## Optional future source manifests

The lab checks for per-cluster JSON manifests under:

```text
pbuf/data/baryonic_source_inputs/<cluster_id>.json
```

No manifest is required for this audit. Missing manifests are expected until independently sourced baryonic data are deliberately added.

Recognized future fields include independently supplied stellar/gas maps, distances, physical pixel area, LOS/deprojection information, or an already absolute baryonic density/surface-density map. The lab does not manufacture any of them.

## Hard guardrails

- no observed kappa pixel values;
- no shear or lensing morphology;
- no lensing amplitude fit;
- no historical `strength=0.18`;
- no solving baryonic mass from `G` or the native transfer;
- no fitted stellar M/L;
- no fitted gas fraction;
- no fitted LOS depth;
- no automatic `z -> distance` conversion;
- no Quantum Engine;
- no Planck-scale input;
- gravity is not fundamental in PBUF;
- stdout only; no run directory is created.

## Run

Checkout PR branch:

```text
foundation/cluster-baryonic-source-input-inventory001
```

From repository root run exactly:

```bash
PYTHONPATH=. python pbuf/labs/foundation/cluster_baryonic_source_input_inventory001.py
```

## Runner contract

The runner is an **executor only**.

Do not modify the lab, source files, redshifts, status classifications, manifests, thresholds, or production modules.

Do not add a cosmology, distances, stellar M/L values, gas fractions, deprojection depths, or baryonic maps during the run.

Do not repair a failure. If execution exits nonzero, return the complete raw failure unchanged.

Return exactly:

1. current HEAD SHA and branch name;
2. process exit code;
3. complete raw stdout and stderr;
4. `git status --short` after the run;
5. confirmation that the preservation stash was not altered.

Do not delete, clean, move, modify, or commit historical untracked `runs/...` directories.

Do not pop, apply, drop, rewrite, or otherwise alter the preservation stash.
