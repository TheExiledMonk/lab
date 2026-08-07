# Matter Input Audit

**Claim:** WEAK-LENSING-OBSERVATION-001 used ρ = max(κ, 0) / max(κ) as the matter-density input to the frozen Version A constitutive law.

**Verdict:** `approximation`

**Justification:**

- In standard weak lensing, κ = Σ / Σ_crit where Σ is the projected surface mass density (baryonic + dark) and Σ_crit is the critical density that depends on the cosmology and the lens/source redshift pair. Hence the published κ is a *mass map in units of the critical density*, not the bare matter density ρ.
- Treating max(κ, 0) as the bare matter density ρ ignores the Σ_crit normalisation and the cosmology that fixes it. It also ignores the baryonic / dark partition; the reconstructed κ integrates both.
- The non-negativity clamp (max(κ, 0)) and the peak-normalisation (/ max(κ)) are defensive choices that suppress the negative tails of the reconstruction (which arise from noise and from mass-deficit regions). They are necessary for Version A's constitutive law to receive a non-negative input, but they are not a physical identification.
- An alternative matter input (e.g., X-ray gas density, stellar mass map, or Σ = κ Σ_crit evaluated at a chosen cosmology) would carry explicit physical units. None is supplied by the benchmark.

**Required for full justification:** An external Σ_crit(z_l, z_s, cosmology), or an explicit baryonic / dark matter partition, is not provided in the benchmark and would have to be supplied by an external cosmological module. Under the frozen laboratory no such module is admitted.
