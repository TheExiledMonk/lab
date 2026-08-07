# WL-001 clean-room experiment

The completed PBUF EM-TRANSPORT-001 audit is recorded in
`runs/em_transport001/em_transport001_report.md`. It tests whether the
neighbour-to-neighbour propagation law required for weak-lensing
wavefront evolution is already implied by the V11 electromagnetic
microscopic structure (`alpha_resolved ~ 3 alpha_EM`, the CORE-001
`q in R^3` with `g_dev = 1/137`). It finds Outcome B: the
microscopic structure supplies static neighbour coupling only;
under the CORE-001 overdamped local evolution the coarse field obeys
the time-independent Helmholtz equation, so no wavefront follows.
The exact missing local physical principle is the kinetic sector
already identified by INERTIA-001 (positive momentum density or an
equivalent symplectic structure), and it cannot be derived from
`alpha_EM`, `alpha_resolved`, or `g_dev` alone. Regenerate the
artifacts with:

```bash
python em_transport001.py
```

The completed PBUF TRANSPORT-RESEARCH-001 comparative audit is
recorded in `runs/transport_research001/transport_research001_report.md`.
It analyses six required local wave-transport systems
(water surface, elastic solid, acoustic, Maxwell EM, spin wave,
plasma / MHD) using a common eight-mechanism framework and finds
two structural families. Family A (mechanical: S1, S2, S3, S5, S6)
requires a restoring mechanism plus inertial resistance; Family B
(Maxwell: S4) requires neither and uses mutual curl coupling of
`E` and `B`. PBUF matches Family A on the spatial half (local state,
neighbour coupling, restoring, steering) but is missing the
resistance / inertia slot; Family B is structurally incompatible
with CORE-001. The closest match is the elastic-solid pattern (S2),
whose completion requires the same kinetic sector flagged by
INERTIA-001. Outcome B. Regenerate the artifacts with:

```bash
python transport_research001.py
```

`PBUF LENS-LOAD-001` is implemented as an executable identifiability gate in
`lens_load001.py`, with its decision report in
`runs/lens_load001/lens_load001_report.md`. Run `python lens_load001.py` to
audit the archived Lens 001 inputs and domain. The frozen corpus does not select
the medium-to-metric map needed to infer native placement from weak lensing, so
the gate deliberately does not fabricate a reconstructed load. An independently
known placement may be mechanically reconstructed and forward-checked with the
explicitly conditional `--placement` interface.

The completed PBUF INVERSE-SOURCE-001 methodology is recorded in
`runs/inverse_source001/inverse_source_report.md`. It proves that a known
placement and boundary work uniquely determine the required generalized
native load, while weak-lensing observations alone do not identify that
placement because the frozen metric map is unselected and the optical inverse
has nontrivial null directions. It defines the set-valued inverse problem,
identifiability theorem, generalized-load target, dimensionless comparison
variables, multi-system universality test, interpretation classes, unchanged
forward-validation pipeline, and readiness decision without fitting data,
adding a coupling, changing the constitutive law, or modifying V11.

The completed PBUF MATTER-MEDIUM-INTERACTION-001 physical-origin audit is
recorded in
`runs/matter_medium_interaction001/matter_medium_interaction_report.md`. It
proves that matter must be an organization or distinction of the complete
one-medium state rather than a second substance, and that native elastic
loading can only mean internal generalized placement-work transfer into the
single frozen load channel. It also proves that the frozen ontology does not
select deformation, energy, waves, stress, or boundary action as the universal
driver and therefore leaves one normalized, gauge-basic universal
matter–placement work postulate necessary before the native source projection
can be uniquely derived. It introduces no new ontology, constitutive law,
constant, metric map, fit, or V11 change.

The completed PBUF EVOLUTION-LAW-001 structural derivation is recorded in
`runs/evolution_law001/native_law_of_successive_state_evolution.md`. It proves
that the strongest frozen-compatible evolution object is a gauge-basic
selection of oriented unparameterized admissible histories, with a set-valued
state continuation correspondence as its projection. It separates ontic
completeness from Markovity, determinism, probability, and reversal, states the
conditional conservation and wave gates, and identifies the exact missing
history-selection role without adding state variables, momentum, coordinate
time, or a proposed kinetic law.

The completed PBUF INERTIA-001 physical-origin audit is recorded in
`runs/inertia001/inertia_origin_report.md`. It proves that the frozen rate-free
elastic energy cannot determine inertia, while stable reversible wave support
conditionally requires a positive momentum-carrying kinetic sector. It defines
inertia intrinsically as a calibrated history/tangent-to-cotangent response,
derives the emergent-time, gauge, positivity, conservation, and V11 causal-cone
gates on any future `K_tau`, and leaves its normalization and final operator
explicitly open. The accompanying CSV files provide the continuum survey and
machine-readable minimum requirements; `validation.json` records completion.

The completed PBUF EQUILIBRIUM-001 audit is recorded in
`runs/equilibrium001/equilibrium_report.md`. It reduces established equilibrium
principles to static constrained stored-energy variation and stationary native
history action, proves that neither selects the nonlinear remainder or D1–D3
endpoint class, and recommends a parametrically constitutive governing-equation
family. Regenerate the report and machine-readable audits with:

```bash
python equilibrium001.py
```

The completed PBUF ENERGY-SEARCH-001 classification is recorded in
`runs/energy_search001/energy_search_report.md`. It reduces named formula
families to one invariant fixed-reference-2-jet function class, audits all
frozen constraints, proves the native eliminations, and identifies the exact
remaining freedom without ranking a law or changing weak lensing. Regenerate
the report and machine-readable catalogues with:

```bash
python energy_search001.py
```

The completed PBUF DURATION-001 derivation is recorded in
`runs/duration001/emergent_duration_derivation.md`. It derives physical duration
as a positive additive reparametrization-invariant line functional of
propagation-bearing medium evolution, distinguishes the order label `s`, clock
duration `tau`, and effective coordinate time `t`, defines physical clocks and
the static limit, and gives the conditional V11 proper-time matching relation.
It also records the irreducible calibration family: the accepted inputs do not
select `tau=ell/c`, a microscopic standard clock, or the missing normalized
medium-to-metric map.

Regenerate and validate the artifacts with:

```bash
python duration001_derivation.py --output runs/duration001
```

The completed PBUF DYNAMICS-001 derivation is recorded in
`runs/dynamics001/native_evolution_principle.md`. It defines the native action
on oriented histories modulo monotone reparametrization, derives the degree-one
homogeneity requirement, tangent/cotangent evolution variables and canonical
pairing, and proves why an ordinary quadratic `T-W` action is not native before
a relational clock is derived. It supplies a dependency graph, equation
traceability, machine-readable action catalogue, and validation record without
deriving field equations or changing V11.

The canonical frozen PBUF ontology is recorded in
`runs/foundation001/foundational_ontology.md`. FOUNDATION-001 establishes FP-1
through FP-6 as explicit framework postulates for subsequent development and
provides a downstream citation contract. The postulates are design inputs, not
derived theorems, and may be reopened only by an explicitly authorized ontology
review, a demonstrated contradiction, or a failed prediction. FOUNDATION-001
is distinct from the unavailable historical `FND-001` artifact and does not
reconstruct its missing alpha traceability.

The completed PBUF IDENTITY-001 ontology audit is recorded in
`runs/identity001/identity_ontology_report.md`. It finds Outcome B, with Outcome
C retained at the cross-history boundary: identifying the spacetime medium with
the physical universe removes a medium/universe duality and prevents automatic
reification of state spaces, while identity plus the inherited continuum result
gives one complete realized state per admissible state position. It does not
select single-state, growing-history, or whole-history actuality. The report,
comparison matrix, dependency graph, consequence register, and recommendation
introduce no physics or change to V11.

`pbuf_experiment.py` runs the single-lens forward pipeline described in
`docs/`: baryonic mass → Version A deformation → gradient → photon paths →
predicted image, alongside isolated baryonic-GR and LCDM reference outputs.

Run it with:

```bash
python pbuf_experiment.py --output runs/wl001
```

The output directory contains CSV fields, PNG visualizations, `run.json`
provenance/checksums, and `execution_log.csv`. Version A is the `version_a`
function in `constitutive_equations.py`; selecting another equation changes only
that independently replaceable component.

Run the WL-002 catalogue (Versions A-D) with:

```bash
python wl002_catalogue.py --output runs/wl002
```

This produces a complete artifact directory per equation plus an RMSE table,
comparison plot, machine-readable catalogue, findings, and recommendation.
The WL-002A discovery record is `runs/wl002/discovery.json`; it ranks the
experimentally evaluated scalar laws, records topology/trajectory/stability/cost
diagnostics, and documents why diffusion and tensor-elastic candidates cannot be
admitted without constants or interfaces absent from the frozen laboratory.

Run the WL-003 physical-derivation audit with:

```bash
python wl003_derivation.py --output runs/wl003
```

It writes a conditional conservation-form derivation, symbol-level traceability,
the exact missing PBUF micro--macro closure, and a frozen-laboratory
reproducibility comparison against archived Version D. Because the supplied PBUF
materials do not quantitatively define that closure, the audit records Outcome B
and deliberately leaves Version D unchanged and labelled empirical.

Run the MB-001 micro--macro closure audit with:

```bash
python mb001_closure.py --output runs/mb001
```

MB-001 audits every supplied PBUF relation, emits equation-level traceability,
and compares an unchanged frozen-laboratory rerun with archived WL-003 Version D,
including deformation, gradients, photon paths, topology, residuals, and RMSE.
Because the supplied theory names possible microscopic inputs but defines no
quantitative coarse-graining or response law, it records Outcome C and precisely
specifies the physical law needed before the constitutive model can be updated.

Run the CORE-001 microscopic-state formalization with:

```bash
python core001_definition.py --output runs/core001
```

CORE-001 defines a three-component dimensionless lattice state, a normalized
rotationally symmetric coarse-graining operator, and the conditional mapping to
the MB-001 continuum coefficients. It emits a report, machine-readable model,
traceability matrix, and executable consistency checks. The three degrees of
freedom and the `1/137` coupling scale remain explicitly labelled hypotheses;
the frozen weak-lensing laboratory is not changed.

Run the FND-002 microscopic-state justification audit with:

```bash
python fnd002_justification.py --output runs/fnd002
```

FND-002 assigns one of the four required labels to every CORE-001 assumption,
tests alternative microscopic realizations, consolidates the irreducible
postulates, and recommends a regulator-independent effective model and an
identifiability-focused FND-003. It is theory-only and does not alter the frozen
weak-lensing laboratory.

Run the FND-003 three-dimensional microscopic-state audit with:

```bash
python fnd003_three_dimensional_justification.py --output runs/fnd003
```

FND-003 distinguishes the three-spatial-dimensions ontology from the additional
faithful-vector assumptions needed to obtain exactly three microscopic
components. It audits rotation-covariant mappings to the scalar continuum field,
shows that post-V11 `g_dev` directly normalizes the corrected CORE-001 source, emits the required
traceability and remaining-postulate artifacts, and records Outcome C without
importing or changing the frozen weak-lensing laboratory.

Run the FND-004 consequence audit with:

```bash
python fnd004_consequences.py --output runs/fnd004
```

FND-004 accepts the three-component ontology and `1/137` coupling as explicit
axioms, derives their exact counting and equal-coupling consequences, and keeps
continuum, propagation, and weak-lensing claims conditional on their additional
requirements. It emits axiom-level traceability, alternative-ontology and unique
prediction catalogues, a derivation report, and executable completion checks
without fitting data or changing the frozen weak-lensing laboratory.

Run the FND-005 experimental-consequence audit with:

```bash
python fnd005_experimental_consequences.py --output runs/fnd005
```

FND-005 emits a no-fit prediction catalogue, observable-to-axiom traceability,
scalar/three-component/generic-N comparison, and a ranked experimental test
programme. It recommends component-resolved bright/dark response tomography as
the first falsification experiment and explicitly separates ontology-fixed
multiplicities from signatures requiring representation, dynamics, or photon
coupling. The frozen weak-lensing laboratory is not imported or changed.

Run the PHOTON-001 microscopic-to-photon coupling derivation with:

```bash
python photon001_derivation.py --output runs/photon001
```

PHOTON-001 derives the general conditional isotropic geometrical-optics action,
ray equation, small-deformation limit, optical phase, and effective-metric
representation. It identifies the missing photon response `n(u)` (especially
`beta=(dn/du)|_0`), audits compatibility with the current scalar WL interface,
and emits a no-fit PHOTON-002 validation specification. It neither imports nor
modifies or executes the frozen weak-lensing laboratory.

Run the CONS-001 top-down coupling-consistency audit with:

```bash
python cons001_consistency.py --output runs/cons001
```

CONS-001 replaces the fixed microscopic coupling premise by symbolic `g_dev`,
traces it through every completed theory sector, and emits a dependency graph,
sector constraint matrix, observable classification, and consistency-overlap
analysis. It finds that `g_dev` is currently unconstrained by cross-sector
consistency even though it directly normalizes the microscopic source, because
the micro--macro and photon maps remain independently unspecified. It performs no observational fit and does not
import, modify, or execute the frozen weak-lensing laboratory.

Run the ERR-001 corrective audit with:

```bash
python err001_correction.py --output runs/err001
```

ERR-001 removes the non-PBUF auxiliary coupling from CORE-001 and all affected
downstream derivations, regenerates and revalidates CORE-001 through CONS-001,
and records corrected equations, a milestone audit, change log, revised
dependency graph, and revised CONS-001 conclusion. The obsolete inverse-scaling
argument is withdrawn; `g_dev` remains without a top-down numerical bound only
because no completed sector supplies a value-selecting or closed cross-sector
constraint.

The PBUF NATURE-001 natural-analogy audit is recorded in
`runs/nature001/nature001_report.md`. It separates cross-material mechanical
behaviors from material-specific mechanisms, supplies the required behavior
matrix and analogy audit, and recommends structural constraints on a future
stored-energy functional. It introduces no equation, parameter, fit, or change
to the frozen weak-lensing laboratory.

Run the V11-ALPHA-001 source-bounded audit with:

```bash
python v11_alpha_audit.py --output runs/v11_alpha001
```

When the authoritative preprint is available, pass it explicitly, for example:

```bash
python v11_alpha_audit.py \
  --v11-preprint docs/Planck-Bound_Unified_Framework_v11_preprint.pdf \
  --output runs/v11_alpha001
```

It enforces the V11-AUDIT-ERR-001 authoritative-source gate. Without an
explicit `--v11-preprint` it stops with `BLOCKED – PRIMARY SOURCE NOT
SUPPLIED`, emits no scientific traceability or later-development comparison,
and does not treat the milestone brief as V11 evidence. A separate equation
set and official amendments can be supplied with `--equation-set` and repeated
`--errata` arguments.

The completed ALPHA-ARCH-001 restoration is in
`runs/alpha_arch001/alpha_architecture_report.md`. It preserves the four-part
V11 hierarchy, renames the unrelated exploratory matter-state coefficient to
`g_dev`, and supplies the notation map, occurrence traceability, dependency
graph, genuine-deviation register, and validation record. FND-001 is explicitly
recorded as absent from the supplied repository; the frozen weak-lensing
implementation is unchanged.

Run the ARCH-001 `g_dev` disposition audit with:

```bash
python arch001_gdev_audit.py --output runs/arch001
```

ARCH-001 inventories every repository occurrence, audits each distinct
equation, and tests mappings to the authoritative V11 alpha hierarchy. Its
decision supersedes the earlier convenience-based retention: no V11 mapping is
derived, normalized component ratios eliminate the common magnitude exactly,
and the absolute matter-vertex branch must be removed or redesigned rather
than retained as fundamental. No replacement parameter is introduced and the
frozen weak-lensing implementation remains unchanged.

Run the MATTER-001 matter-action derivation audit with:

```bash
python matter001_derivation.py --output runs/matter001
```

MATTER-001 derives the stress-energy tensor as the universal matter operator
from V11's retained standard-GR architecture and obtains the parameter-free
medium source as a functional derivative of the matter action. It identifies
the exact remaining gap as the absent local elastic action and normalized
medium-to-metric/coarse-graining map. It introduces no replacement coupling,
does not fit data, and does not import, execute, or modify weak-lensing code.

The completed GEOMETRY-001 metric-emergence audit is in
`runs/geometry001/metric_emergence_derivation_report.md`. It derives the exact
functional chain from an arbitrary covariant medium-to-metric map to the
stress-energy source, classifies the admissible metric families, and identifies
the missing normalized constitutive identification of microscopic clocks and
rulers with one Lorentzian metric. Its decision is Outcome D with mathematical
non-uniqueness C. It introduces no coupling, changes no V11 alpha definition,
performs no fit, and does not import, execute, or modify weak-lensing code.

The completed SOURCE-PROJECTION-001 audit is in
`runs/source_projection001/source_projection_report.md`. It proves that the
frozen ontology fixes the native load's codomain, virtual-work role, units, and
symmetry constraints but not the physical matter-to-load projection. It also
proves the point-local isotropic tensor obstruction, classifies the exact
remaining operator freedom, and identifies a normalized universal
matter–medium interaction/virtual-work principle as the minimal missing
ingredient. It introduces no field, coupling, fit, metric map, or V11 change.
