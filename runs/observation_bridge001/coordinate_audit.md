# Coordinate Audit

## Published

- **system:** Equatorial RA / Dec on WCS (TAN projection)
- **projection:** gnomonic (TAN)
- **pixel_scale:** 6.25 - 11.36 arcsec / pixel (cluster dependent)
- **origin:** Field centre at CRVAL1, CRVAL2; CRPIX = (N/2 + 0.5, N/2 + 0.5)
- **handedness:** RA increases to the EAST; Dec increases to the NORTH; pixel (0, 0) is at the BOTTOM-LEFT of the array
- **orientation:** CD matrix is diagonal (no rotation); CD1_1 < 0 (RA flips east-to-west in pixel index)

## Version_A

- **system:** Dimensionless Cartesian grid
- **projection:** Identity (no projection)
- **pixel_scale:** 16 / 128 = 0.125 dimensionless units / pixel (matter); 16 / 64 = 0.25 dimensionless units / pixel (observables)
- **origin:** Pipeline centre at (0, 0)
- **handedness:** x increases to the right; y increases up; pixel (0, 0) is at the BOTTOM-LEFT of the array
- **orientation:** Identity (no rotation)

## Mismatch

- Published WCS is astronomical (RA/Dec); Version A is dimensionless Cartesian. The two cannot be aligned without imposing an external angular scale.
- WEAK-LENSING-OBSERVATION-001 mapped pixel index (0, N-1) -> (-extent, +extent) and discarded all angular information. The cluster's true angular diameter (1500 arcsec for four of the five benchmarks; 1200 arcsec for MACS1149) was not used.
- WCS handedness was not preserved: the published CD1_1 < 0 sign flip was dropped.

