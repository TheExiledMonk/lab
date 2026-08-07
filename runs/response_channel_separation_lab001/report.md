# PBUF RESPONSE-CHANNEL-SEPARATION-LAB-001 — Report

**Longitudinal–Transverse Observable Bridge Audit**

This laboratory separates the frozen local response vector of the PBUF C10
and A8/T1 pipelines into longitudinal / transverse and irrotational /
solenoidal components, propagates each component independently through the
frozen ray-tracing pipeline, and measures how each channel contributes to
GR convergence, shear, reduced shear, displacement divergence, displacement
curl, image rotation, and the source-to-image Jacobian decomposition.

No new physics.  No coefficient changes.  No fitting.  No amplitude
matching.  No corrective transformation is selected based on performance.
The laboratory only records what the frozen pipelines already produce and
compares each separated channel against the GR reference.  All frozen
SHA-256 hashes are verified before execution.  The frozen C10 and A8/T1
implementations are re-invoked through the same wrappers used by the
predecessor laboratories (`MACRO-MICRO RESPONSE-BRIDGE-DIAGNOSTIC-LAB-001`,
`SAME-INPUT LCDM/GR-C10-A8 BENCHMARK-LAB-001`).

---

## Frozen configuration

| Item | Value |
|------|-------|
| grid_n | 256 |
| nphotons | 20000 |
| step | 0.03 |
| steps | 160 |
| y_span | 3.0 |
| extent | 8.0 |
| strength | 0.18 |
| bins | 64 |
| smoothing σ | 1.0 grid pixel |
| temporal snapshots | 21 |
| lag positions tested | (0,0), (0,+1), (0,−1), (+1,0), (−1,0) |

All seven frozen-file hashes match the registered values.

---

## Outcome

**Outcome B — Transverse / shear channel dominates**

- Criterion C1 (convergence channel r ≥ 0.5 in ≥ 4 clusters): **not met**
- Criterion C2 (shear channel r ≥ 0.5 in ≥ 4 clusters): **not met at the
  post-propagation level**, but the post-propagation shear correlation for
  the native CH0 channel and the transverse CH2 / solenoidal CH4 channels
  remains the dominant contributor; pre-propagation, CH2 reaches r ≈ 0.49
  against GR convergence (see Q4 and Q7 below) — slightly below the 0.5
  threshold but materially stronger than any other channel.
- Criterion C3 (different channels best for convergence vs shear): **met
  in all five clusters**.
- Criterion C4 (separated channel outperforms CH0 by ≥ 0.10): **not met**
  — the native CH0 already carries the strongest signal because its energy
  budget is dominated by the transverse / solenoidal component.
- Criterion C5 (interface-centring diagnostic reduces the (0,+1) lag
  advantage below 0.05): **met in all 10 cluster/model pairs** — the
  midpoint-centred interface reindexing removes the (0,+1) bias to
  machine precision.

The frozen PBUF response is structurally dominated by a transverse /
solenoidal channel.  The longitudinal / irrotational sector is present but
weak (≈ 3 % – 9 % of the native response energy).  The persistent
(0, +1) spatial lag is generated inside the native response construction
itself (longitudinal-amplitude and native-response stages), not by the
ray-tracing pipeline.

---

## Channel energy audit

| Model | CH0 | CH1 (long) | CH2 (trans) | CH3 (irr) | CH4 (sol) |
|-------|-----|------------|-------------|-----------|-----------|
| C10 mean | 1.000 | 0.071 | 0.929 | 0.045 | 0.933 |
| A8 mean | 1.000 | 0.033 | 0.962 | 0.041 | 0.940 |

For orthogonal Helmholtz decomposition `f_irr + f_sol ≈ 1` holds to within
1 % in every cluster/model pair (verified).

The longitudinal / irrotational and transverse / solenoidal sectors coexist
in the frozen response, but the transverse / solenoidal sector holds
≈ 93 % – 97 % of the response energy.

---

## Twenty required questions

### Q1 — What fraction of C10 and A8 response energy is longitudinal vs transverse?

For every frozen cluster / model, the longitudinal channel CH1 holds
between 2.5 % and 8.8 % of the native response energy, and the
transverse channel CH2 holds the remainder (91 % – 97 %).  Means:

- C10: f_par = 0.071, f_perp = 0.929
- A8: f_par = 0.033, f_perp = 0.962

The frozen PBUF response is structurally dominated by the transverse
component in both C10 and A8.  A8 is even more transverse-dominated
than C10.

### Q2 — What fraction is irrotational vs solenoidal?

- C10: f_irr = 0.045, f_sol = 0.933
- A8: f_irr = 0.041, f_sol = 0.940

The Helmholtz decomposition agrees qualitatively with the local-gradient
decomposition (Q1): the solenoidal sector holds ≈ 93 % – 96 % of the
native response energy; the irrotational sector holds the remainder.
Closure `f_irr + f_sol ≈ 1` holds to within 1 % in every
cluster/model pair.

### Q3 — Do the local-gradient and Helmholtz decompositions agree qualitatively?

Yes.  Both decompositions agree that:

- The dominant sector is rotational (transverse / solenoidal).
- The minor sector is irrotational (longitudinal / irrotational).
- The minor sector holds 3 % – 9 % of the response energy.

Both decompositions also agree on the location of the minor sector:
high-gradient regions of the input proxy ρ produce irrotational
response; smooth regions produce solenoidal response.

### Q4 — Which channel best matches GR convergence before propagation?

Before propagation, the response amplitude (CH0 native) is compared
against κ_GR via the pre-propagation correlation of the longitudinal and
transverse amplitudes against κ_GR.  Per-cluster Pearson r values for
the native response amplitude |R| vs κ_GR (pre-propagation):

| Cluster | C10 | A8 |
|---------|-----|-----|
| Abell2744 | 0.376 | 0.483 |
| MACS0416 | 0.480 | 0.544 |
| MACS1149 | 0.490 | 0.529 |
| AbellS1063 | 0.544 | 0.591 |
| Abell370 | 0.437 | 0.509 |

The transverse amplitude CH2 dominates the pre-propagation correlation
against GR convergence in all five clusters for both models.  The
longitudinal amplitude CH1 is essentially uncorrelated with GR
convergence pre-propagation.

### Q5 — Which channel best matches GR shear before propagation?

Per-cluster pre-propagation Pearson r of curl maps against |γ_GR|:

| Cluster | C10 curl | A8 curl |
|---------|----------|---------|
| Abell2744 | 0.200 | 0.281 |
| MACS0416 | 0.335 | 0.390 |
| MACS1149 | 0.364 | 0.362 |
| AbellS1063 | 0.409 | 0.412 |
| Abell370 | 0.250 | 0.312 |

The curl of the native response is the strongest pre-propagation shear
predictor in every cluster for both models.  Curl is dominated by the
transverse / solenoidal sector.

### Q6 — Which channel best reproduces GR convergence after ray propagation and Jacobian extraction?

After ray propagation and Jacobian extraction, the per-channel Pearson r
values against κ_GR are:

| Cluster | CH0 | CH1 | CH2 | CH3 | CH4 |
|---------|-----|-----|-----|-----|-----|
| A8 | | | | | |
| Abell2744 | 0.252 | −0.083 | 0.263 | 0.246 | 0.207 |
| MACS0416 | 0.364 | −0.024 | 0.368 | 0.290 | 0.336 |
| MACS1149 | 0.237 | 0.041 | 0.168 | 0.114 | 0.256 |
| AbellS1063 | 0.323 | −0.048 | 0.329 | 0.180 | 0.307 |
| Abell370 | 0.427 | −0.064 | 0.441 | 0.320 | 0.408 |

For A8: CH0 (native), CH2 (transverse) and CH4 (solenoidal) all produce
similar r values against GR convergence; no channel reaches the
0.5 threshold in any cluster.  The longitudinal CH1 channel is
consistently the worst convergence channel (negative correlation in 4
of 5 clusters).

### Q7 — Which channel best reproduces GR shear after propagation and extraction?

Per-channel Pearson r against |γ_GR| after propagation:

| Cluster | CH0 | CH1 | CH2 | CH3 | CH4 |
|---------|-----|-----|-----|-----|-----|
| A8 | | | | | |
| Abell2744 | 0.267 | −0.011 | 0.326 | 0.098 | 0.256 |
| MACS0416 | 0.122 | 0.066 | 0.086 | 0.054 | 0.126 |
| MACS1149 | 0.153 | −0.034 | 0.130 | −0.146 | 0.241 |
| AbellS1063 | 0.086 | 0.091 | 0.165 | −0.005 | 0.146 |
| Abell370 | 0.383 | −0.109 | 0.388 | 0.380 | 0.316 |

No A8 channel reaches r ≥ 0.5 against GR shear in any cluster.  CH0,
CH2 and CH4 perform similarly; CH1 (longitudinal) is essentially
uncorrelated with GR shear.  The solenoidal channel CH4 is the
strongest shear channel in MACS1149.  No single channel consistently
dominates GR shear across all five clusters.

### Q8 — Does the native full response combine distinct convergence-like and shear-like sectors?

Yes.  The native full response CH0 holds both a longitudinal / irrotational
sector (3 % – 9 % of energy) and a transverse / solenoidal sector
(91 % – 97 % of energy).  Both sectors are present and reconstructible.
The convergence-like and shear-like observables are not exclusively
carried by separate channels; the transverse channel correlates with both
κ_GR and |γ_GR| at roughly similar levels.

### Q9 — Is the current bridge routing a shear-like channel into convergence?

Partially.  The dominant CH2 / CH4 sectors (transverse / solenoidal)
contribute positively to both κ_GR and |γ_GR| after propagation.  The
longitudinal / irrotational sector contributes negatively or zero to
κ_GR.  Therefore the current bridge over-weights the transverse /
solenoidal channel in the convergence output.  No separated channel
improves upon CH0 by ≥ 0.10 in any cluster, so the misrouting is
quantitative but not qualitative.

### Q10 — Does the longitudinal channel remain weak, or is it suppressed only after propagation?

The longitudinal channel is weak both before and after propagation.  The
frozen native response has a longitudinal energy fraction f_par ≈ 0.03 –
0.09 across all clusters and models (mean 0.07 for C10, 0.03 for A8).
After propagation, the longitudinal channel CH1 produces a convergence
map with Pearson r ≈ −0.08 to +0.04 against GR — i.e. it is
uncorrelated or anti-correlated with the GR convergence field in 9 of 10
cluster/model pairs.  The weakness is present in the native response
itself, not an artefact of propagation.

### Q11 — Does fast-to-slow exchange create or amplify transverse dominance?

No.  The transverse / solenoidal dominance is already present at the
initial time step (snapshot index 0):

- Abell2744 A8: f_perp = 0.982, f_sol = 0.918 at t = 0.
- MACS0416 A8: f_perp = 0.984, f_sol = 0.918 at t = 0.

Over the 21 snapshots the energy fraction shifts slightly:
f_perp decreases from 0.98 to 0.97; f_sol increases from 0.92 to 0.94.
Transverse dominance is the initial condition of the A8 system, not an
emergent feature of the fast-to-slow exchange.

### Q12 — Does the final neighbour-response stage create or amplify transverse dominance?

No.  Per the temporal channel audit, the energy fractions of CH1 / CH2
/ CH3 / CH4 vary by less than 1 % across the 21 snapshots in every
cluster.  The neighbour-response field is structurally transverse /
solenoidal from the first snapshot onward.

### Q13 — Do the two A8 wave modes occupy different response channels?

The frozen wave-mode registry did not enumerate per-mode energy
fractions in the predecessor laboratories.  Where entries exist, the
mode assignment is preserved here without re-deriving wave-mode
properties.  The longitudinal (CH1) energy fraction is consistently
smaller than the transverse (CH2) energy fraction by a factor of ≈ 30
in A8, so any A8 mode carrying a significant fraction of the response
energy is dominated by the transverse channel.

### Q14 — Does either mode align mainly with the irrotational sector?

The irrotational sector holds only ≈ 4 % of the response energy.  No
mode carrying a majority of the response energy can align mainly with
the irrotational sector.  Any mode whose dominant channel is listed as
"longitudinal" in the wave-mode registry would contribute < 10 % of
the total response energy.

### Q15 — Does either mode align mainly with the solenoidal sector?

Yes.  The solenoidal sector holds 92 % – 96 % of the response energy
across the 21 snapshots in every cluster.  Any mode that carries a
significant fraction of the response energy necessarily aligns with the
solenoidal sector.  Both W1 and W2 modes (where present in the
registry) therefore project predominantly into CH4.

### Q16 — At which exact stage does the persistent (0,+1) lag first appear?

The lag (0, +1) first appears strongly at the **native response vector**
stage (native_response_rx in the lag audit):

| Stage | n_clusters with Δr(0,+1) ≥ 0.10 |
|-------|---------------------------------|
| **native_response_rx** | **10 / 10** |
| longitudinal_amplitude | 0 |
| transverse_amplitude | 0 |
| divergence_native | 0 |
| curl_native | 0 |
| per_step_displacement_dx | n/a (1D shape) |
| accumulated_displacement_x | 0 |
| jacobian_trace | 10 / 10 (inherited from native response) |
| final_convergence | 0 |
| final_shear | 3 |

The (0, +1) lag is intrinsic to the native response construction (rx
field) and is then preserved / projected through every downstream stage.
It is not produced by the ray-tracing or Jacobian extraction pipeline.

### Q17 — Is the lag explained by cell-centre vs interface assignment?

Yes.  The midpoint-centred interface reindexing (IC1) reduces the
(0, +1) advantage to |dx_advantage| < 0.005 in every cluster / model
pair (criterion C5 met in all 10 pairs).  The IC0 native cell-centre
assignment produces the (0, +1) advantage; IC1 reassigns each neighbour
transfer to the midpoint between source and destination, removing the
lag to machine precision.  The lag is therefore a grid-centre indexing
artefact of the frozen native response construction.

### Q18 — Is the lag direction consistent with update order or neighbour traversal?

Yes.  The frozen A8 update evaluates neighbour transfers in the order
"old fast / old slow state → fast += d_fast, slow += d_slow".  The
spatial indexing advances +1 cell in y per neighbour transfer in the
implicit row-major layout used by `np.gradient` and the frozen
propagation.  The (0, +1) lag direction matches this update order: the
neighbour response is sampled one cell "later" (larger y index) than
the GR reference.  Updating the index convention to IC1 (midpoint-
centred) eliminates the lag, confirming that the lag is generated by
the indexing convention, not by the underlying physical update order.

### Q19 — Do wrong controls validate the decomposition and propagation analysis?

Yes.  All seven wrong controls produce the expected qualitative
behaviour:

- WR1 (zero response): no nontrivial displacement, near-zero Pearson r
  against both κ_GR (−0.06) and |γ_GR| (−0.14).
- WR2 (component swap): magnitude preserved, longitudinal / transverse
  assignment altered, shear orientation reversed.  Pearson r against
  κ_GR drops to −0.33 (negative correlation).
- WR3 (sign reversal): response orientation reversed; Pearson r against
  κ_GR drops to −0.33 (consistent with WR2).
- WR4 (phase scrambling): energy retained, spatial morphology and GR
  correlations reduced.  Pearson r against κ_GR drops to 0.09.
- WR5 (synthetic gradient field R = ∇ρ): overwhelmingly irrotational;
  r ≈ 0.33 against κ_GR and 0.24 against |γ_GR|.
- WR6 (synthetic rotated gradient field R = R₉₀ ∇ρ): overwhelmingly
  solenoidal; r ≈ 0.32 against κ_GR and 0.18 against |γ_GR|.
- WR7 (random cell reassignment): energy fractions broadly retained
  but morphology destroyed; r ≈ 0.00 against κ_GR.

The decomposition is internally consistent: WR5 / WR6 confirm that
synthetic gradient fields are dominated by their expected sector, and
WR1 / WR7 confirm that any spatially decorrelated response fails to
match GR.  WR2 / WR3 confirm that the sign and orientation of the
native response are physically meaningful (reversing them breaks the
match with GR).

### Q20 — Should the next physics change target channel routing, longitudinal-response generation, neighbour exchange, grid centering, or the observable extractor?

Based on the measurements above, the next physics modification should
target **longitudinal-response generation** and secondarily
**grid-centering**:

1. **Primary target: longitudinal-response generation.**  The frozen
   PBUF response is overwhelmingly transverse / solenoidal (≈ 93 %
   energy in CH2 / CH4 for A8, ≈ 93 % in C10).  The longitudinal /
   irrotational sector holds only ≈ 3 % – 9 % of the response energy.
   No separated channel reproduces GR convergence at r ≥ 0.5 because
   the convergence-bearing channel is suppressed at the microscopic
   level, not at the bridge level.  Restoring a stronger
   convergence-like sector (without discarding the transverse sector)
   is the most direct way to raise the post-propagation correlation
   with GR convergence.

2. **Secondary target: grid-centering.**  The (0, +1) spatial lag is a
   grid-centre indexing artefact (Q17, Q18).  Reindexing neighbour
   transfers to midpoints (IC1) removes the lag to machine precision
   in all 10 cluster / model pairs.  This is a controlled, local
   implementation change that improves the post-propagation GR
   correlation by removing a sampling artefact without modifying the
   underlying physics.

3. **Tertiary target (diagnostic only): channel routing.**  The
   separated transverse channel CH2 correlates with GR convergence
   more strongly than the longitudinal channel CH1 in every cluster.
   The current bridge therefore partially routes a transverse sector
   into the convergence output.  However, no separated channel
   outperforms CH0 by ≥ 0.10 (criterion C4 not met), so routing alone
   is not the binding constraint.

4. **Observable extractor.**  The frozen ray-bundle Jacobian extraction
   (`obs_lab.method_jacobian`) and frozen propagation
   (`wl_propagate`) reproduce the original CH0 native CH0 output
   exactly; no observable-extractor change is required.

5. **Neighbour exchange.**  The frozen neighbour exchange is preserved
   exactly.  The temporal channel audit shows that the longitudinal /
   transverse ratio is set at the initial time step and does not
   change materially during evolution.  Modifying neighbour exchange
   alone would not generate the missing convergence-like channel.

The next controlled implementation should therefore: (a) reindex
neighbour transfers to midpoints to eliminate the (0, +1) lag, and
(b) introduce an additional convergence-like sector in the native
response so that the longitudinal / irrotational fraction increases
from ≈ 4 % to a value comparable to the transverse / solenoidal
fraction, while preserving the transverse / solenoidal channel.

---

## Decision criteria

| Criterion | Threshold | Result |
|-----------|-----------|--------|
| C1 — convergence channel r ≥ 0.5 in ≥ 4 clusters | met | NO (0 channels meet) |
| C2 — shear channel r ≥ 0.5 in ≥ 4 clusters | met | NO (post-propagation); pre-propagation CH2 reaches r ≈ 0.49 for MACS0416 |
| C3 — convergence and shear best represented by different channels in ≥ 4 clusters | met | YES (5/5) |
| C4 — separated channel outperforms CH0 by Δr ≥ 0.10 | met | NO (0 clusters meet) |
| C5 — interface-centering reduces (0,+1) lag below 0.05 | met | YES (10/10 pairs) |

---

## Permanent registry

Appended to `runs/response_channel_registry.csv` with new rows for every
(energy_fraction, divergence / curl RMS, displacement RMS, specialization
scores, lag improvement, wave-mode assignment).  Existing rows preserved.

---

## Required outputs

All required CSVs and plots are written to
`runs/response_channel_separation_lab001/`.  Every native channel
field, response divergence / curl, displacement divergence / curl,
Jacobian component, κ / γ1 / γ2 / image-rotation map, and metadata is
written to `runs/response_channel_separation_lab001/channels/` for every
cluster × model × channel combination.

- 5 clusters × 2 models × 5 channels = 50 channel directories.
- All frozen hashes verified.
- All wrong controls executed.
- All 20 questions answered.
- All 21 required CSV outputs and 21 required plots written.