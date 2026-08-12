# PBUF Historical Attempt Index

This is an anti-circularity index, not a canonical-result ledger. Results may
be historical, invalidated, superseded, or preserved. The development ledger
alone records canonical status.

| Mechanism | Historical dev or PR | Dimension | Result | Canonical status | Why noncanonical (if applicable) | Current follow-up | Do not retest unless |
|---|---|---|---|---|---|---|---|
| Scalar gradient propagation | Pre-ledger | surrogate | steering demonstrated | historical only | assumes response from scalar gradient | vector-pair dynamics | independently derived law requires it |
| Tangent-stiffness propagation speed | Pre-ledger | scalar/N6 diagnostic | explored | historical only | changes law with loading | none | load-dependent law is independently derived |
| Scalar induced geometry | Pre-ledger | scalar field | explored | rejected by Dev164 | scalar does not provide bond orientation/embedding | Dev167 explicit relations | independent geometry state exists |
| Mass-loading speed coupling | Dev145 | native scalar | unresolved | unresolved | no state-to-speed law | none | such a law is independently derived |
| Shared-state cross-coupling | Dev151–153 | N6 scalar/link | insufficient | rejected | longitudinal link state did not derive transverse response | richer relation state | independent multicomponent law exists |
| Loaded transverse link response | Dev151–153 | N6 link | insufficient | rejected | renamed cross-coupling | richer relation state | independent multicomponent law exists |
| 1D frame transport | Dev151–153 | 1D surrogate | explored | noncanonical | not full N6 | full N6 only | independently justified as diagnostic |
| Scalar loaded F03 | Dev163 | N6 scalar | loaded perturbation cancels | rejected/frozen | unchanged operator gives free F03 | vector-pair branch | governing operator independently changes |
| Scalar geometry | Dev164 | N6 scalar | lengths/directions underived | rejected/frozen | scalar differences are not geometry | explicit vector relations | independent relational state exists |
| H07 allocation | Dev165 | N6 allocation | partial redirection | partial, nongoverning | allocation not derived | structural comparison only | derived dynamics produces it |
| Minimal binary magnetic pair | Dev165 H12 | N6 binary pair | gates failed | rejected | binary polarity insufficient | richer magnetic-like space open | richer reversible law exists |
| DEV169 invariant-drift forensic audit | Dev170 | N6 accounting | persistent source potential omitted from medium-only invariant | resolved bookkeeping | not a pair-law failure | source-work invariant semantics | source force or invariant definition independently changes |
| DEV169 depth-sensitivity provenance audit | Dev170 | native→observer | mixed native direction/receipt variation and observer amplification | resolved mixed | does not establish unique source depth | source-depth uniqueness | new independent 3D source constraint |
| DEV169 deposition-semantics audit | Dev170 | observer interface | loaded summary was stale signed-channel value; arrays reconcile exactly | resolved | label/provenance issue | serialized-array semantics | serializer or summary path changes |
| DEV169 nonfinite-channel audit | Dev170 | observer bank | knn_density convergence undefined with duplicate launch points / zero kNN radius | expected degeneracy | no implementation defect | retain unchanged bank | input geometry or kNN implementation changes |
| Arbitrary 2D→3D depth lane selection | Dev162/169 | diagnostic source family | depth sensitivity demonstrated | diagnostic only / must not be promoted | arbitrary guessed geometric depth is not an independent reconstruction | Dev171 frozen independent catalog ensemble | a new independent data set materially changes depth constraints |
| Spectroscopic redshift as direct geometric depth | Dev171 | source reconstruction | forbidden | forbidden | peculiar velocity is degenerate with geometric position within the cluster | phase-space-only component inference | independently justified dynamical reconstruction establishes such a relation |
| Post-freeze WL morphology comparison | Dev172 | frozen observer vs. pyRRG-JWST WL map | blocked before comparison | blocked / no result | public code repository did not provide an A2744 map/catalogue asset with WCS; native observer grid has no serialized WCS | freeze released WL asset and independently justify a native-grid-to-sky WCS bridge | new independent observational dataset or independently derived physical normalization/observable bridge changes the comparison basis; do not retest by changing source depth, source weights, native pair law, packet law, receipt, observer channel, registration, or smoothing |
| Source-to-observer coordinate lineage recovery | Dev173 | source→observer interface | catalog-to-native footprint deterministic; serialized 6×6 provenance dropped | resolved serialization boundary | do not solve lost provenance by source-depth retuning, packet-law changes, observer-channel changes, manual WCS/registration, or reopening native excitation | serialize receipt/adapter/grid-footprint provenance | new independent coordinate provenance or external WCS asset |
| Observer coordinate-provenance serialization | Dev174 | source→receipt→observer serialization | all eight frozen 6×6 outputs have native and deterministic sky footprints | closed native coordinate package | finite-cell / 6×6 discretization remains explicit; no formal WCS claimed | obtain independent astrometric observational asset | new observational asset or independently justified physical normalization |
| Published pyRRG asset recovery + blind morphology | Dev175 | frozen observer vs. Harvey–Massey 2024 pyRRG-JWST | exhaustive public release-provenance search did not recover exact numerical map/catalogue | blocked / external unavailable | no authentic product with deterministic astrometry; comparison not run | contact authors or predeclare a separate independent WL dataset | authentic released asset or independently predeclared observational dataset; never retune PBUF/registration/smoothing/mask/metric |
| Direct receipt spin-2 candidate matrix (P1--P7) | Dev176 | frozen DEV168 receipt / DEV174 footprints | P1--P6 implemented without observational values; sparse per-footprint support and native-control replays remain incomplete | partial / not promoted | no valid JWST pyRRG PSF-corrected shape catalogue is frozen; structural gates incomplete | complete native controls and acquire/freeze authentic direct shape catalogue | a frozen receipt with materially denser support or a new independently frozen observational shape asset; never introduce sign/rotation/registration/retuning |
| Full native received-state channel audit | Dev177 | individual DEV168 3D receipt records | preserves geometry, displacement, direction, momentum, flux, content and lineage; audits complementarity and intrinsic `J3/G3` before collapse | canonical observer-input boundary | does not define what an image-forming observer measures | derive a deterministic physical observer mapping from the retained state | `FULL_NATIVE_RECEIPT_AUDIT_COMPLETE` and an independently justified observer mapping; do not invent a new 2D extractor merely from a diagnostic ranking |
| Vulkan/KDE and source-plane coverage infrastructure | historical PRs c620cd9 / b54caa8 and DEV178 | historical A8/M10 ray plane; current DEV167/168 receipt diagnostics | exact Vulkan KDE, CPU/GPU parity, viewer patterns, deterministic coverage reference recovered | infrastructure reusable; old physics historical only | historical five-lens A8/M10 propagation, processed kappa proxies, 0.18 strength and legacy observer fusion remain noncanonical | apply Vulkan/viewer only to frozen DEV167/168/177 receipt artifacts; use historical 25% geometry solely as a sampling reference | a physically defined current-native launch representation exists; do not reactivate processed five-lens physics or introduce 2D observer extractors before distribution and observer-mapping closure |
| DEV179 native sub-cell source representation | DEV179, corrected by DEV180 | DEV167 integer source contact / DEV168 receipt lineage | no physical off-node source-to-medium mapping is present | valid narrow closure: `OUTCOME_D` | absence from the current implementation was incorrectly generalized to absence from PBUF history; historical C25 rays were launches, not matter locations | recover and reconcile historical PBUF source-medium coupling before any new law is derived; determine whether the desired density is packet/launch, source, or grid resolution | an independently derived current-native bridge is required only for the relevant unresolved density/source representation; never reopen by bilinear/trilinear/spline/Gaussian/nearest-node interpolation or by receipt/J3/viewer/observational performance |

New 2D observer extractors are allowed only if `FULL_NATIVE_RECEIPT_AUDIT_COMPLETE` and `PHYSICAL_OBSERVER_MAPPING_REQUIRES_IT` are both true.

Historical continuous G3D launch density cannot be transferred directly to DEV167/168 without a current-native source-position representation. Sub-cell coordinates do not become physical information merely because interpolation can assign values; high-density launches remain blocked until source-to-medium representation closure.

DEV182 reconciliation: executable PR16 C25 is `267×267=71,289`, not the stale nominal `266×266=70,756`; executable PR106 is `534×534=285,156`. PR107 establishes that actual launch-coordinate arrays are primary receipt lineage. Current DEV167 packet launches are finite discrete node supports under exact-reset replay, not continuous historical ray coordinates.

DEV183 closure: longitudinal packet translations change the fixed launch-to-receipt experiment and are not density samples. The current same-experiment domain is the complete 121-state transverse periodic integer-translation set; every scaled replay requires additive, lossless packet/launch/realization receipt provenance. Continuous packet interpolation and output-selected launches remain forbidden.

DEV184 result: the frozen 1/30/60/121 nested exact-reset packet ladder was completed across all eight DEV171 realizations. The preserved DEV177 full-state rank and depth-rank diagnostics are stable by C50 in all eight lanes; C100 is retained as the canonical finite-domain production coverage. This is density sufficiency evidence only, not an observer/channel selection result.

DEV185 result: the complete frozen C100 receipt distribution is structurally sufficient across all eight realizations (full raw-feature rank 14; DEV177 depth increment 4 in every realization). This authorizes only the next native mode/channel sufficiency audit; no channel, statistical component, DFT mode, or observer is promoted.

DEV186 result: only the exact storage-derived content columns W01, W03, and W04 are reconstructable from retained momentum, W02, and flux; the reduced representation preserves all eight rank/depth gates. PCA and DFT remain diagnostic bases, and the raw native basis retains physical priority. Physical observer derivation—not observational comparison—is now authorized.

DEV187 result: DEV168 positive receipt-crossing weight is an additive native detector measure in exact native receipt cells, but the C100 121-launch collection is a sampling of propagation transfer response rather than literal brightness multiplicity. A source-image pushforward remains required; DEV176 P1–P7 and J3/G3 are not promoted, and shape/spin-2/observational gates remain closed.

## Registry integration (DEV181)

The mechanism registry at `docs/PBUF_MECHANISM_REGISTRY.json` is now the primary historical lookup layer. This Historical Attempt Index remains a concise anti-circularity index and must not be expanded into a duplicate of the full registry. Every row above maps to one or more target/attempt records in that registry.

DEV188 anti-circularity rule: C100 equal launch coverage is transfer sampling, not uniform physical brightness unless external source weights specify it. The native pushforward is a finite position-dependent operator, not a convolution kernel.

DEV189 anti-circularity rule: the observed lensed galaxy image cannot be reused as the incident source for its own PBUF prediction. WCS is observed-image astrometry, not an unlensed source-plane map; an angular-to-native scale may not be fit to weak-lensing morphology.

DEV190 anti-circularity rule: full multichannel source rank does not establish a physical distortion metric. SVD/PCA/DFT and K^T K remain diagnostic without a native output-space metric; transfer geometry requires explicit launch/receipt coordinates and never source fitting or observational selection.

DEV191 anti-circularity rules: conventional weak-lensing terminology does not establish a PBUF bridge; PSF-corrected moments and calibrated shear are different layers; representation matching is not physical equivalence; ensemble claims must declare intrinsic-source assumptions.

DEV192 anti-circularity rules: a new lens cannot repair an observable-definition failure on the development lens; cross-lens validation begins only after observer freeze; intrinsic-orientation cancellation must be derived and cannot be silently assumed.

DEV193 anti-circularity rules: linear combination of independently sampled transfer columns does not prove physical simultaneous-source superposition in a nonlinear native medium; a centroid Jacobian does not automatically transport full image morphology.

DEV194 anti-circularity rule: additive detector columns cannot be rebranded as independent physical events merely because they were computed by exact reset; stationarity, event injection, and post-event memory must be established separately.

DEV195 rule: a nonzero global post-event state cannot alone establish persistent local memory; local matched-background force balance, outward transport, and periodic recurrence must be separated.

DEV196 rule: additive DEV182 packet initialization is valid on an evolved VectorPairState, but exact floating-state support must be classified before interpreting a residual as a disjoint-event failure.

DEV197 rule: exact support overlap is diagnostic geometry, while sequential-event significance is a channelwise comparison of B-after-A against B-fresh and deterministic replay.

DEV198 rule: residual local native force and B influence are separately measured over the fixed pre-recurrence sequence; exact-zero denominators stay undefined.

DEV199 rule: four-state force residuals are exact DEV167 bond identities; receipt differences are path-integrated local consequences, not evidence for a new nonlocal rule.

DEV200 rule: canonical excitation propagation is an evolving N6 bond-state plus momentum pattern; scalar force reductions and packet labels are not dynamically sufficient.

DEV201 rule: the canonical DEV167 N6 central-force tangent is not the rejected DEV151 longitudinal-link representation; nevertheless its unloaded spectrum contains no robust propagating transverse sector.

DEV202 rule: inspect self-generated DEV195 bond loads with the full DEV167 central-force tangent; do not repeat DEV151 scalar longitudinal-link failure or impose prestress.
