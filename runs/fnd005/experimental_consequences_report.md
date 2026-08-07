# PBUF FND-005 — Experimental consequences of the microscopic ontology

## Result

The minimal ontology has **no nontrivial positive laboratory observable without an access map**: component count is latent until components can be prepared or read out. Once that requirement is stated, its cleanest falsifiable prediction is a three-dimensional response space. With the equal common coupling added, source-response tomography must contain one bright combination and exactly two source-dark combinations; a calibrated coherent-to-single-channel power ratio must equal 3. These counting relations do not require a constitutive propagation law or fitted parameter.

The often-associated one-longitudinal/two-transverse spectrum is less fundamental. It additionally assumes that the three components are a spatial vector and that isotropic dynamics exists. Branch speeds, gaps, damping, conservation laws, and weak-lensing profiles are not fixed by A1–A3.

No parameter was fitted, no new ontology was adopted, and the frozen weak-lensing laboratory was neither imported nor changed.

## Axioms and experimental assumptions

- **A1:** The microscopic state has exactly three independent components.
- **A2:** A common matter load couples equally to those components with bare scale g_dev=1/137.
- **A3:** Macroscopic observables arise by coarse graining the microscopic state.

- **E1:** The three components can be independently prepared or read out with a calibrated linear map.
- **E2:** The common loading and readout normalization are fixed independently; linear response applies.
- **E3:** The components form a spatial vector in an isotropic, parity-even reference state.
- **E4:** A stable local quadratic long-wave theory and a specified time kinetic law apply.
- **E5:** Symmetry breaking is a controlled perturbation rather than an unknown apparatus effect.
- **E6:** A photon/effective-metric coupling maps microscopic or coarse fields to light.

E1–E6 are validation conditions, not consequences silently added to the ontology. A failed conditional prediction falsifies PBUF only when its listed conditions have been independently established.

## Prediction catalogue and observable-to-axiom traceability

| ID | Category | Prediction | Axioms | Extra assumptions | Ontology only | Observable | Comparison class | Falsifier/boundary |
|---|---|---|---|---|---|---|---|---|
| F005-P01 | mode counting | Component-resolved state tomography has rank three. | A1 | E1 | False | Three independently preparable/readable response directions. | unique to PBUF against scalar or N!=3 | Resolved rank is not three after sensitivity and constraints are accounted for. |
| F005-P02 | coherent versus incoherent loading | One normalized common source excites one bright combination and leaves exactly two orthogonal combinations dark at the source vertex. | A1,A2 | E1,E2 | False | A rank-one common-source response with two null source channels. | unique multiplicity: scalar has 0 dark; generic N has N-1 | The calibrated common vertex has nullity other than two or unequal entries. |
| F005-P03 | coherent versus incoherent loading | The coherent bright-mode amplitude is sqrt(3) times one equally normalized component and its quadratic weight is three times one component; equivalently, the squared unnormalized coherent sum is three times the incoherent sum of three equal component powers. | A1,A2 | E2 | False | No-fit ratio R_amp=sqrt(3), R_single-power=3, or R_coherent/incoherent=3 under the stated normalization. | distinguishable by multiplicity; generic N gives sqrt(N), N and scalar gives 1, 1 | Any normalized ratio differs beyond preregistered uncertainty. |
| F005-P04 | longitudinal/transverse mode structure | An isotropic spatial-vector realization has one longitudinal and two transverse polarizations. | A1 | E3,E4 | False | One L eigenvector and a two-dimensional T eigenspace for nonzero k. | unique versus scalar; not implied for a generic internal N-state | A stable isotropic vector sector lacks the 1+2 multiplicity. |
| F005-P05 | degeneracy relations | The two transverse modes are exactly degenerate before symmetry breaking. | A1 | E3,E4 | False | Two equal transverse poles/rates at fixed \|k\|. | distinguishable from scalar; generic N is representation-dependent | Reproducible splitting remains in an isotropic limit. |
| F005-P06 | propagation branches | If inertial dynamics applies, three polarizations occur as one L branch plus a double T branch; overdamped dynamics gives the same multiplicities as relaxation poles, not waves. | A1,A3 | E3,E4 | False | A 1+2 pole pattern, without predicted speeds or gaps. | conditional PBUF realization, not unique from A1-A3 | Only conditional after E3-E4 are independently established. |
| F005-P07 | symmetry breaking signatures | Weak anisotropy may split the transverse doublet; restoring isotropy must restore degeneracy. | A1 | E3,E4,E5 | False | Reversible transverse splitting correlated with a controlled anisotropy. | generic vector signature, not unique to exactly three components | Nonzero intercept in the controlled isotropic limit, subject to systematics. |
| F005-P08 | coupling identifiability | g_dev directly normalizes the common matter vertex, while normalized multiplicity ratios cancel its magnitude. | A2 | E1,E2 for an absolute measurement | False | A calibrated absolute vertex is g_dev-sensitive; bright/dark counts and normalized ratios are not. | direct coupling is distinguishable only with absolute calibration | A calibrated equal vertex inconsistent with the stipulated g_dev falsifies A2; ratios alone cannot measure its magnitude. |
| F005-P09 | conservation consequences | No conservation law follows from component count, equal loading, or coarse graining alone. | A1,A2,A3 | none | True | A conserved charge requires an action and continuous symmetry not present in the axioms. | identical boundary for all compared ontologies | Not applicable: this prevents an unsupported prediction. |
| F005-P10 | weak-lensing equivalence | If only one scalar projection reaches photons, scalar, three-component, and generic N ontologies can be observationally identical. | A3 | E6 and decoupled hidden components | False | No component-count discriminator in scalar lensing alone. | identical across ontologies | Not a PBUF signature; observation of extra coupled polarizations would break equivalence. |

## Comparison against alternative ontologies

| Observable | PBUF three-component | Scalar | Generic N | Verdict |
|---|---|---|---|---|
| resolved response-space rank | 3 | 1 | N | unique to PBUF only relative to N!=3; requires E1 |
| dark channels under one equal common source | 2 | 0 | N-1 | unique multiplicity; requires E1-E2 |
| normalized coherent amplitude/power | sqrt(3) / 3 | 1 / 1 | sqrt(N) / N | unique multiplicity; requires E2 |
| longitudinal/transverse multiplicity | 1 L + 2 T | one scalar | undefined without representation | conditional on spatial-vector realization E3-E4 |
| transverse degeneracy | twofold | absent | representation-dependent | symmetry signature, conditional rather than ontology-only |
| branch speeds, gaps, damping | undetermined | undetermined | undetermined | identical lack of prediction without constitutive dynamics |
| scalar weak-lensing profile | undetermined | undetermined | undetermined | identical without closure and photon coupling |
| separate value g_dev=1/137 | direct vertex premise; measurable only with calibrated response/readout | model-dependent | model-dependent | no auxiliary coupling degeneracy; numerical value remains postulated |
| conserved microscopic charge | not implied | not implied | not implied | identical without action and continuous symmetry |

Scalar coarse observations can hide any number of decoupled components. Consequently, matching a scalar continuum or weak-lensing curve cannot determine microscopic component count. Conversely, finding `N=3` is not logically exclusive to the PBUF name; it supports the stated three-component ontology against the comparison classes, while the equal-coupling and representation tests supply additional discrimination.

## Ranked experimentally testable signatures

| Rank | Signature | Discrimination | Dependencies | Readiness | Rationale |
|---|---|---|---|---|---|
| 1 | Bright/dark source-response tomography: one equal bright channel plus exactly two dark channels | N-1 dark channels directly counts N; PBUF predicts 2 | A1,A2,E1,E2 | Simulation-ready once component source/readout operators exist; laboratory platform not yet specified | No propagation coefficients or parameter fit are needed; it tests counting and equal coupling together. |
| 2 | Coherent-to-single-channel power ratio of 3 | Scalar=1; generic N=N | A1,A2,E2 | Simulation-ready and potentially experimental with calibrated channels | Dimensionless no-fit ratio, but vulnerable to channel normalization and cross-talk. |
| 3 | One longitudinal plus a degenerate transverse doublet | Separates spatial-vector PBUF from scalar; generic N depends on representation | A1,E3,E4 | Requires a dynamical realization and directional spectroscopy | Sharp multiplicity test, but not forced by the minimal ontology. |
| 4 | Transverse splitting vanishes with controlled anisotropy | Tests vector symmetry, not exactly N=3 by itself | A1,E3,E4,E5 | Requires tunable anisotropy and systematic-error control | Useful corroboration after the vector representation is established. |
| 5 | Additional non-scalar imprint in photon observables | Could break scalar equivalence but is not presently predicted | A3,E6 plus a closure | Not test-ready | No photon coupling, profile, or amplitude follows from the supplied ontology. |

## Recommended first falsification experiment

Run a preregistered **component-resolved source-response tomography** in the first microscopic simulation (or physical platform) that supplies calibrated preparation and readout operators:

1. Determine the response matrix using independent small-amplitude probes, with rank and detection thresholds fixed before examining the result.
2. Test whether its resolved state space has rank 3.
3. Apply the normalized common load `g=g_dev(1,1,1)` and rotate the measured basis to `q_B=(q1+q2+q3)/sqrt(3)` plus two orthogonal directions.
4. Without refitting, test one bright source vertex, two null vertices, equal component entries, amplitude ratio `sqrt(3)`, and power ratio `3` relative to one channel.
5. Repeat across scale and resolution. Count a failure only where all three modes remain within sensitivity and calibration/cross-talk controls exclude a hidden or constrained channel.

This experiment is first because it tests the exact multiplicities inherited from A1–A2 while avoiding unknown stiffnesses, dispersions, photon coupling, and weak-lensing assumptions. Outcomes rank 1, 2, or greater than 3 favor the corresponding scalar or generic-N alternatives over PBUF; unequal common-source entries falsify A2. Until a microscopic simulator or material system defines E1, this is a fully specified protocol class rather than an immediately executable bench experiment.

## Constitutive and phenomenological boundary

Ontology/counting fixes the numbers `3`, `2`, and—after calibrated equal linear coupling—`sqrt(3)` and `3`. It does not fix pole locations, speeds, damping, spectral weights after propagation, nonlinear response, symmetry-breaking magnitude, conserved charges, lensing deflection, or experimental platform. Those quantities must not be fitted and then reported as consequences of the ontology.

## Completion checks

All required categories, deliverables, comparisons, traceability fields, and the first falsification recommendation are present. All automated checks pass: **True**.
