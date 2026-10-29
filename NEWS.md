# BreastSubtypeR (development)

## Bug fixes

* `RORgroup()` keeps both proliferation panel subsets as matrices and no
  longer fails when exactly one 11-gene or 18-gene panel gene is represented
  (#206). Represented genes are still averaged with `mean(..., na.rm = TRUE)`
  without any gene-count threshold; absent genes are never synthesized. A
  sample with no observed panel value, or a nonfinite panel estimate, now
  receives a missing proliferation score instead of `NaN` or `Inf`, so
  ROR-P, original ROR-PC and the additional ROR-PC score and groups are
  missing for that sample rather than the literal text `"NaN"` or an
  apparently valid `high`/`low` category. Every exported numerical ROR score
  is now finite or `NA`: a `NaN` or infinite result, including overflow in
  the final `100 * (score + 0.35) / 0.85` scaling or an infinite tumor-size
  covariate, is reported as `NA` with a missing group for ROR-S, ROR-P,
  ROR-C, original ROR-PC and the additional subtype-specific ROR-PC score.
  Finite complete-input scores, subtype calls, calibration, correlations,
  coefficients, cutoffs, the exact-equality limitation (#172), the reference
  standardization question (#198), the correlation choice (#156) and the
  independent NODE >= 4 high-risk rule are unchanged. Incomplete panels are
  documented as a compatibility limitation, not a validated score.

# BreastSubtypeR 1.5.11

## Mapping

* `Mapping()` now defaults to `impute = FALSE`; primary examples and Shiny
  explicitly opt out. Duplicate resolution precedes optional imputation and
  complete required outputs bypass donor preparation.
* Optional nonraw imputation uses the fixed-k=10 pure-R estimator over the full
  eligible gene-level donor input, targeting only original holes in represented
  genes selected for existing NC/SSP outputs. Exact ties retain donor row order;
  estimates use one to ten available original contributors per coordinate among
  those ten selected profiles.
  Unsupported/partial targets and nonfinite newly projected estimates error
  with internal diagnostics. Absent rows are not created and observed/nontarget
  values, RNG state, output structure and expression scales are preserved.
* Optional KNN no longer uses legacy row/column missingness cutoffs. Unusable
  background donors are excluded with diagnostics, required genes remain
  protected, and diagnostic capacity is checked before estimation.
* Reject missing raw counts with either setting and never impute raw-derived
  expression. Full-count normalization and library-size handling are unchanged.
  Remove unused `impute` imports/dependency after replacing the legacy backend.

* Infinity rejection applies when optional KNN runs, not to disabled or complete
  output bypasses. This revision differs from the original published pipeline;
  incomplete-input estimates and downstream predictions may change even with
  explicit `impute = TRUE`.

## Bug fixes

* AIMS and the SSPBC subtype models now return missing calls and probabilities
  for samples with zero evaluable model rules (#212, PR #215). Observed FALSE
  comparisons still count as evidence; supported samples and model definitions
  are unchanged.

# BreastSubtypeR 1.5.10

## Bug fixes

* Preserve all Shiny reference-median labels under C and UTF-8 locales
  (#210, PR #211). Construct the four non-ASCII option names at runtime;
  all option values, ordering and defaults remain unchanged. Regression
  tests compare the rendered UI and HTTP option labels across locales.

* Load Shiny server and UI regression-test sources explicitly as UTF-8
  (PRs #209 and #211). This fixes locale-dependent test-source parsing
  without weakening the citation assertions or changing prediction logic.

* `BS_Multi()` restores AUTO `ssBC` and `ssBC.v2` subcohort results by
  character `PatientID` labels after checking complete, unique identities
  (#202). Factor levels no longer cause eligible calls to be lost during
  restoration. Intended missing calls, valid character-ID results, AUTO
  eligibility, prediction algorithms and ROR formulas are unchanged.

* `BS_cIHC.itr()` now converts `PatientID` values to character before
  selecting expression columns in each iteration (#201). Factor labels,
  including unused levels, select the intended ER-balanced centering samples
  and values. Missing, empty, duplicated and unmatched identifiers produce
  informative errors. Valid character-ID results, both subsampling policies,
  iteration counts, RNG behavior and numerical ROR formulas are unchanged.

* `BS_cIHC()` now converts `PatientID` values to character before
  name-based indexing (#200). Factor labels, including those with unused
  levels, select the intended ER-balanced centering samples and expression
  values. Missing, empty, duplicated and unmatched identifiers produce
  informative errors. Valid character-ID results, balancing rules, RNG
  behavior and numerical ROR formulas are unchanged.

* `BS_PCAPAM50()` now converts `PatientID` values to character before
  name-based indexing (#197). Factor labels, including those with unused
  levels, select the intended expression samples and centering values.
  Missing, empty, duplicated and unmatched identifiers produce informative
  errors. Valid character-ID results, sampling order, output order, RNG
  behavior and numerical ROR formulas are unchanged.

* `BS_PCAPAM50()` now orders complete intermediate rows by character
  `PatientID` before second-stage sampling, matching the PCAPAM50 1.0.4
  first-step merge under the caller's collation (#195). Selected sets,
  centering medians, calls and derived scores may change. First-step
  calculations and final input-order restoration are preserved. This does
  not resolve the separate ROR expression-standardization investigation.

* `BS_PCAPAM50()` now restores all sample-indexed outputs to input order,
  including `score.ROR` and the prediction, distance and expression components
  of `outList`. Internal calibration and sampling, and per-sample values for
  the same input, seed and environment, are unchanged.

* Report an informative error from `Mapping(method = "max")` when every
  probe of a multi-probe Entrez ID has an undefined row total (`NA`, `NaN`
  or `Inf - Inf`), naming the Entrez ID and the number of candidate probes,
  instead of failing with an unrelated dimnames error. This changes only the
  diagnostic: no probe is selected, imputed, averaged or dropped for such a
  group, and selection and results for every group with a defined total are
  unchanged. Follows the singleton preservation of #178; the duplicate
  resolution review of #177 remains open.

## Interpretation and documentation

* Declare `S4Vectors` and `htmltools` in `Suggests` for their direct use in
  tests and the optional Shiny application.

* The factor-identifier corrections (#197, #200, #201 and #202) can change
  affected sample selections, predictions, derived scores or restored calls.
  The PCA-PAM50 intermediate-order correction (#195, PR #196) can change
  results even with valid character identifiers. Unchanged ROR formulas do
  not imply unchanged scores when their inputs change. Output-axis
  restoration (PR #193) preserves per-sample values but can affect downstream
  positional joins. These fixes do not resolve the separate ROR scientific
  questions or establish that manuscript results are unaffected.

* Clarify current version and installation-channel wording in the README
  and vignette. Earlier release notes are retained.

# BreastSubtypeR 1.5.9

## Bug fixes

* `RORgroup()` now requires finite nonnegative integer-valued `NODE` counts
  for `ROR-PC Group (Subtype + Clinic + Prolif.Subtype)` (#187). Negative,
  fractional or infinite numeric values previously left the unscaled score
  as numeric text in the group column or were silently classified; they now
  give one warning per call and a missing group. The column is consistently
  character, including all-missing batches. `NA`/`NaN` `NODE`, all numerical
  ROR scores, the original Parker groups and the valid-count rules from #157
  are unchanged. The exact-equality limitation in #172 remains deferred.
* Preserve singleton gene rows under `Mapping(method = "max")` without
  ranking them (#178). In the nonraw `impute = FALSE` path, remaining `NA`
  and `NaN` measurements are retained instead of causing a dimension error.
  Selection among multiple probes and all other mapping methods are unchanged.
* Reject mixed AUTO/manual selections and duplicate method names in
  `BS_Multi()` before accessing input data or executing classifiers (#188).
  AUTO must be the only entry; manual selections must contain at least two
  distinct methods. Valid method selections and their execution order are
  unchanged.

# BreastSubtypeR 1.5.8

## Bug fixes

* Correct PAM50 `Confidence` gene alignment, Spearman-test handling and
  two-decimal rounding; undefined values remain missing. Earlier code could
  lose centroid gene names when extracting a data.frame column, fail to pair
  genes and return false zeros. This defect is verified in the Bioconductor
  3.23 version 1.4.1 and 3.24 version 1.5.7 source archives. GitHub revisions
  sharing a version number can differ: the Confidence fix was already in
  GitHub commit 6f68e09 (version 1.5.7). Users of affected builds should
  recompute Confidence summaries and downstream uses of those summaries.
  A legitimately computed Confidence can still be zero. This finding alone
  does not invalidate subtype calls or numerical ROR scores, which the fix
  does not change. Confidence is not a calibrated subtype probability.
* Correct the additional NODE-dependent
  `ROR-PC Group (Subtype + Clinic + Prolif.Subtype)` boundaries (#154):
  NODE = 0 uses <= 40 / > 40 through 60 / > 60, and NODE = 1--3 uses
  <= 15 / > 15 through 40 / > 40; at least four positive nodes remain high
  risk with nonmissing tumor size. Scores are not rounded for grouping.
* Represent an undefined numerical `NaN` score as an actual missing category
  in the additional NODE-dependent ROR-PC group for NODE = 0 or 1--3,
  instead of the literal text `"NaN"`. Numerical scores, finite/infinite
  categories, clinical missing-value handling and the independent NODE >= 4
  high-risk rule are unchanged. Original Parker grouping is unaffected.

## Reporting and presentation

* Resolve Shiny preprocessing through the installed package namespace when
  launched with `BreastSubtypeR::iBreastSubtypeR()`. Empty intentional Shiny
  aborts clear stale results without blank error callouts; meaningful
  validation and ordinary error messages remain visible.
* Remove the unused ggrepel import without changing plot geometry. Clarify
  NODE error wording, OSLO sample-subset provenance and `.txt` TSV exports.
* Add descriptive `InputCoverage()` and Shiny coverage reporting for PAM50,
  AIMS and SSPBC, using bundled-model gene/rule denominators. Reporting does
  not change method eligibility, imputation, predictions or scientific
  reliability, and does not assess every ROR component.
* Refresh the five visualization functions' typography, spacing, legends
  and annotations while retaining the established subtype colors. Display
  the already calculated heatmap row dendrogram. These are presentation
  changes: scaling, clustering, displayed order, calls, counts, PCA values,
  boxplot statistics and raw entropy remain unchanged.
* Update the local Shiny presentation, public OSLO2-EMIT0 README screenshot,
  installation/channel guidance, input and export documentation, and
  Editor's Choice links. Clarify ssBC.v2 reference-quantile provenance.
* Keep strict package-error/warning CI gates while excluding account checks
  only from the GitHub CI BiocCheck invocation. Mailing-list subscription,
  Support Site registration and watched-tag registration require separate
  maintainer verification; skipped checks are not reported as successful.

## Scientific clarification and known limitation

* Clarify research Spearman correlations versus commercial Prosigna Pearson
  correlations (#156). Shared coefficients do not establish equivalent
  preprocessing, calibration or clinical validation. Calculations are
  unchanged; these are research outputs, not clinical Prosigna results.
* Temporarily defer the known original Parker-derived ROR-S, ROR-P, ROR-C
  and `ROR-PC (Subtype + Clinic + Prolif)` categorical-output defect at exact
  equality to unscaled cutoffs. Strict comparisons can leave numeric text
  instead of `low`, `med` or `high`. The unscaled cutoff pairs are
  -0.15/0.10, -0.25/0.10, -0.10/0.20 and -0.20/0.20, respectively, before
  `100 * (score + 0.35) / 0.85`. Numeric-text labels are invalid categories,
  not scientifically accepted assignments. Primary evidence does not
  conclusively define equality assignments; a separate scientific review
  is required. Preserve all four grouping implementations, scores,
  coefficients and thresholds. The separate NODE-dependent correction
  above remains in place.

# BreastSubtypeR 1.5.7

## Bug fixes and documentation

- Fixed `BS_AIMS()` for a single mapped sample while retaining the exact
  sample identifier and matrix dimensions. Historical multi-sample calls,
  probabilities and other returned components are unchanged.
- Corrected the AIMS return-value documentation: calls and assigned-call
  probabilities are sample-by-one-column matrices; `all.probs` is a named
  list containing the five-subtype probability matrix.
- Fixed stale results in `iBreastSubtypeR()`. Input replacement, example
  reload and relevant parameter changes invalidate dependent predictions,
  plots and exports. Mapping is retained when only prediction settings change.
- Added server-side download checks for absent, failed or outdated
  predictions. Download controls appear only for supported current results;
  unavailable full metrics are explained rather than replaced with calls.
- Invalid Shiny inputs now produce actionable, recoverable errors. Correcting
  prediction settings permits a rerun without remapping valid inputs or
  reloading the session. Failed predictions cannot expose previous results.
- Clarified the cIHC.itr ER-positive : ER-negative target ratio (default
  `54 / 64`): it must be finite and strictly positive, and target mode permits
  values above one. The Shiny iteration control retains default 100, minimum
  10 and step 10; other whole numbers of at least 10 are accepted.
- Added AIMS and Shiny regression tests. Validation also included real-browser
  checks covering input invalidation, downloads, error recovery, uploaded
  sample identifiers and independent sessions. Synchronized installation,
  workflow and output documentation without changing scientific algorithms
  or the public API.

# BreastSubtypeR 1.5.6

## Bug fixes (behaviour change)

- cIHC.itr (`BS_cIHC.itr()`) now balances ER groups towards the ER+ : ER-
  target `ratio` (default `54 / 64`, the ER distribution of the UNC PAM50
  reference set). In every iteration the ER group that is over-represented
  relative to the target is subsampled without replacement, to
  `ceiling(ratio * n_ERneg)` ER+ or `ceiling(n_ERpos / ratio)` ER- samples,
  and the other group is kept in full; a cohort already at the target is
  used in full. Previously the numerically larger ER group was always
  subsampled, to `ceiling(ratio * n_smaller)` samples, so the balancing
  subset followed the target only when ER+ samples were at least as many as
  ER- samples and was inverted otherwise (for example, 54 ER+ and 64 ER-
  samples were subsampled to 54 : 46 instead of being used in full).
- The `ceiling()` rounding convention is kept, so with the default `ratio`
  results are identical to 1.5.5 whenever ER+ samples are at least as many
  as ER- samples. With the default `ratio = 54 / 64`, the balancing subset
  can change for cohorts or subsets in which ER- samples outnumber ER+
  samples, that is, cohorts with an ER-negative majority. This can affect
  the calls, scores and ROR values of `BS_cIHC.itr()`, and the `cIHC.itr`
  and `entropy` columns of `BS_Multi()`: in manual calls that include
  `"cIHC.itr"`, and in AUTO calls that select cIHC.itr under the existing
  routing and size gates. AUTO selects cIHC.itr in both branches of the
  eligible mixed cohort (at least 15 ER+ and 18 ER- samples, outside the
  TN-only and HER2+ routes), with ten methods when the ER+ fraction is
  strictly between 0.39 and 0.69 and with eight methods otherwise, so
  eligible cohorts with an ER+ fraction of 0.39 or less are included. A
  changed balancing subset need not change the subtype call of every
  cohort or sample. AUTO method selection, its thresholds and all other
  methods are unchanged.
- New argument `er_subsample` in `BS_cIHC.itr()` and `BS_Multi()` (passed
  to cIHC.itr only). `er_subsample = "larger"` restores the previous rule,
  sampling order and outputs exactly, to reproduce results obtained with
  BreastSubtypeR 1.5.5 and earlier; the default is `"target"`.
- `ratio` must be a single finite positive number (write `54 / 64`, not
  `54:64`). `54:64` and `ratio = 0`, which previously gave a base-R error or
  centred on the smaller ER group only, now stop with an informative error.
  With the default `er_subsample = "target"` any such `ratio` is feasible.
- Regression tests compare both modes with outputs retained from
  BreastSubtypeR 1.5.5 and check that `BS_Multi()` changes no other method.

# BreastSubtypeR 1.5.5

## Bug fixes

- Fixed standalone reading of `inst/CITATION` without package metadata.
  The citation file now uses an ASCII-only comment and generates the
  package-manual entry only when package metadata is available. The
  published article citation and the metadata-derived manual entry are
  preserved during normal installed-package use.

# BreastSubtypeR 1.5.4

## Documentation

- The `URL` and `BugReports` fields in DESCRIPTION, the README and the
  vignette now point to the development repository
  https://github.com/JohanHartmanGroupBioteam/BreastSubtypeR and its issue
  tracker instead of the former personal repository, which redirects there.
  The mirroring workflow `.github/workflows/sync-to-org.yml`, obsolete since
  the repository moved to the organization, was removed. No code changes.

# BreastSubtypeR 1.5.3

## Maintainer

- Emmanouil G. Sifakis (emmanouil.sifakis@ki.se) is the package maintainer
  from this version on; Qiao Yang remains an author. The public repository
  moved to https://github.com/JohanHartmanGroupBioteam/BreastSubtypeR (the
  previous address redirects). No code changes.

# BreastSubtypeR 1.5.2

## AUTO routing

- AUTO now decides whether a cohort is HER2+ from the HER2 column alone: the
  HER2+ branch requires at least one evaluable HER2 value and every evaluable
  HER2 value to be `HER2+`. Previously the branch was entered whenever the
  joint ER/HER2 counts held no HER2-negative sample, which had two
  consequences. A cohort with HER2 entirely missing, or coded with
  unrecognized values, was routed as HER2+ and ran `AIMS` and `sspbc` only;
  such cohorts now follow the ER-based rules (for example a 60/40 ER cohort
  gets the balanced mixed panel, with `ssBC.v2` returning `NA`). And a
  HER2-negative sample whose ER value was missing was left out of the joint
  counts, so one such sample did not end the HER2+ classification while an
  otherwise identical sample with a known ER value did; the classification no
  longer depends on whether ER was recorded. Missing values and equivocal
  codes such as `"2+"` are treated alike, as not evaluable, and AUTO now
  reports how many samples had no evaluable HER2 value and therefore took no
  part in the decision. Cohorts with complete HER2 information are routed
  exactly as before. Added tests, including regression tests for the
  size-gated ER/HER2-defined cohorts.
- AUTO sample subsetting for `ssBC` and `ssBC.v2` now uses the same inclusive
  minimums as method selection (`>=`). Previously a subgroup whose size was
  exactly the minimum (ER+ 15, ER- 18, HER2 subgroups 8/9) was accepted by the
  selection step but excluded from the subset, so its samples received `NA`
  from `ssBC`/`ssBC.v2` while the method was reported as running. The packaged
  OSLO2-EMIT0 example is such a case (18 ER- tumors): with the strict
  comparison a fresh `BS_Multi(methods = "AUTO")` run did not reproduce the
  packaged frozen result `OSLO2EMIT0obj$res`, which carries `ssBC` calls for
  those tumors; with `>=` it does. A regression test asserts this. Note that
  the pseudocode in the supplementary methods of the published paper (Yang et
  al. 2025, NAR Genomics and Bioinformatics, doi:10.1093/nargab/lqaf131)
  records the strict comparison in this subsetting step while using `>=` for
  method selection; the package now applies `>=` in both, and the shipped
  example is the behavior the package reproduces.
- README and vignette: the AUTO bullet for ER/HER2-defined cohorts now states
  the size gating (ER group minimum and HER2 subgroup minimum; smaller cohorts
  run AIMS and sspbc only), the HER2+ cohort detection rule, and the treatment
  of cohorts without evaluable HER2 values.

## Bug fixes

- Fixed the `ROR-C Group (Subtype + Clinic)` column in the ROR output of the
  nearest-centroid methods when `hasClinical = TRUE`: it was a copy of the
  `ROR-PC Group (Subtype + Clinic + Prolif)` column instead of the risk group
  derived from the ROR-C score (thresholds -0.1 and 0.2). The ROR-C score and
  all other columns are unchanged. Added a synthetic regression test.
- Fixed the AUTO-mode sample subsets passed to ssBC and ssBC.v2: samples with
  a missing ER (or HER2) value produced NA sample names in `samples_ER.icd` /
  `samples_ERHER2.icd`, which either failed the SummarizedExperiment subsetting
  in `BS_Multi()` ("index out of bounds: NA") or, when the padded vector was as
  long as the cohort, silently skipped the intended subsetting. Missing values
  are now dropped with `which()`, matching `makeCalls.ssBC()`. Cohorts without
  missing ER/HER2 values are unaffected. Added regression tests.
- AUTO now reports when no cohort rule matches the ER/HER2 subgroup sizes and
  it falls back to the single-sample predictors AIMS and sspbc (previously a
  silent fallback). The selected methods are unchanged.
- `Mapping()` now accepts a `SummarizedExperiment` with a single sample.
  `duplicate_genes()`, `prepare_nc_matrix()` and the probe filter dropped the
  matrix dimensions of one-sample input and failed with "dim(X) must have a
  positive length". The collapsed gene-by-sample matrix is now built
  explicitly; results for multi-sample input are identical for all `method`
  values. Present in 1.4.0 and 1.5.1. Added regression tests.
- Phenotype tables with factor `ER`, `HER2` or `TN` columns are now normalized
  exactly like character columns; previously unmatched factor levels were
  replaced by their integer codes (for example "Unknown" -> "3") on the
  `BS_Multi()` / `get_methods()` path, which bypasses the factor conversion in
  `Mapping()`. Added regression tests.
- `Mapping(RawCounts = TRUE, impute = TRUE)`: the FPKM matrix is now checked
  for missing values itself before imputation (the guard tested the already
  imputed log-CPM matrix, so the FPKM matrix was never imputed).
- `BS_Multi()`: the warning issued when PCAPAM50 fails now includes the
  underlying error message (it previously ended after "failed in this
  iteration: ").
- `iBreastSubtypeR()`: the launcher's dependency helper called
  `requireNamespace()` with unsupported arguments, so it failed silently and
  loaded nothing (the app still started because `shiny::runApp()` attaches
  shiny itself). The helper (`.load_app_dependencies()`) now loads the
  requested namespaces and stops with a clear message if a package is missing
  or cannot be loaded. Added a test.
- Shiny app: the cohort preflight reported "ready" when only one of the ER and
  HER2 columns (or only TN) was present, and the AUTO run then failed with
  "requires both 'ER' and 'HER2' columns". The preflight now names the missing
  column(s) and blocks the run. Added tests.
- Documentation: the `@return` sections of `BS_cIHC()`, `BS_cIHC.itr()`,
  `BS_PCAPAM50()` and `BS_ssBC()` now describe the list that is actually
  returned (`BS.all`, `score.ROR`, `mdns`/`mdns.fl`, `outList`, and for
  `BS_cIHC.itr()` the per-iteration call matrices); they previously described
  a character vector, a data.frame or non-existent elements. The
  `BS_cIHC.itr()` `ratio` argument is now documented as applied to the larger
  ER group relative to the smaller one, and the alphabetical tie-break of its
  consensus call is stated. No code changes.
- Documentation: `Mapping(method = "mean")` and `"median"` were described as
  keeping "the probe with the highest mean/median expression"; the
  implementation collapses the duplicate probes of a gene into their per-sample
  mean or median, exactly as `collapseIDs()` in the original PAM50
  bioclassifier code (Parker et al., 2009), from which the `"mean"`,
  `"median"`, `"iqr"` and `"stdev"` options derive. The help text now says so
  and notes that `"max"` (largest row sum) is an addition of this package. No
  code or mapped values changed.
- The per-sample entropy reported by `BS_Multi()` is now `NA` when no executed
  method returned a call for that sample. `table()` drops missing values, so
  such a row previously gave `-sum(numeric(0))`, that is 0 - the same value as
  unanimous agreement. Rows with at least one call are unchanged, and the
  statistic itself (raw, unnormalized Shannon entropy in bits) is untouched.
  The situation cannot arise in AUTO, whose panels always include AIMS and
  sspbc;
  it arises in manual runs of nearest-centroid-only method sets, where a sample
  with unknown ER/HER2 status receives no call. Note that `Vis_Multi()` orders
  by entropy, so such samples now sort last instead of appearing among the
  unanimous ones. The `entropy` column is also documented in `?BS_Multi` for the
  first time. Added tests.
- `BS_PCAPAM50()` no longer treats samples with an unknown ER status as
  ER-negative. The IHC label is derived from the `ER` column, so a missing or
  non-canonical ER value yields `NA`, and the reference test
  `!grepl("^L", IHC)` is `TRUE` for `NA`; such samples were therefore counted
  as ER-negative in the PC1 axis check, in the misclassification-minimizing
  cutoff search and in the ER-balanced gene-centering set. They are now
  excluded from those three steps, a message reports how many were excluded,
  and they are still classified. Cohorts with a complete ER column are
  unaffected: the packaged OSLO2-EMIT0 example reproduces its stored PCAPAM50
  calls exactly. This also aligns `BS_PCAPAM50()` with `BS_cIHC()` and
  `BS_cIHC.itr()`, which already exclude unknown ER from their balancing.
  Added tests.
- `BS_cIHC()`, `BS_cIHC.itr()` and `BS_PCAPAM50()` now stop with a clear
  message when the cohort lacks one of the two ER (or luminal/non-luminal IHC)
  groups needed for ER balancing, instead of failing with "undefined columns
  selected" or "arguments imply differing number of rows". Cohorts with both
  groups are unaffected. Added tests.

## Tests

- `tests/testthat.R` now calls `test_check("BreastSubtypeR")`, so the files
  under `tests/testthat/` run during `R CMD check`. Previously the file held a
  single inline test and `test-tn-cohort-detection.R` was never executed. The
  inline `BS_Multi()` test moved to `tests/testthat/test-bs-multi-manual.R`.
  testthat edition 3 is declared (`Config/testthat/edition`), and the Suggests
  entry requires testthat >= 3.2.0.

# BreastSubtypeR 1.5.1

## Bug fixes

- Fixed TNBC cohort detection in cohort-specific method selection (#133). Cohorts are now classified as TNBC only when all evaluable non-missing `TN` annotations indicate TN. Mixed cohorts containing both TN and nonTN samples are no longer routed to the TNBC-specific branch solely because some samples are TN.

## Tests

- Added synthetic phenotype-table tests for TN-only, mixed, UNC-like mixed, and missing or partly missing TN annotations.

# BreastSubtypeR 1.5.0

## Highlights (from v1.1.3 onward)

- Paper published in *NAR Genomics and Bioinformatics* (2025), **Editor's Choice** (DOI: 10.1093/nargab/lqaf131).
- Support for raw RNA-seq counts (requires gene lengths).
- **iBreastSubtypeR** refresh: cleaner UX, smarter AUTO guidance, consistent exports.

### Enhancements

- **ssBC/ssBC.v2: singleton subgroup robustness:** Subgroups with `n=1` no longer error:
  - Keeps matrix shape (`drop=FALSE`) and hardens dimnames/types.
  - Primary path: original `sspPredict()`.
    Fallback: nearest-centroid (Spearman) when needed.
  - If there are **0 common PAM50 genes**, returns `NA` labels with shaped `distances`/`dist.RORSubtype` to avoid downstream errors.
  - ROR computation guarded for incomplete inputs.
- **SSPBC output now "full":** `BS_sspbc()` and Shiny "sspbc" runs return a full metrics table (not calls-only).
  - Exports map core label columns to the standard names (`Call_5class` / `Call_4class` when applicable).
- **Shiny: "Load example data..." button**
  - One-click load of a small demo dataset from `inst/RshinyTest/` to explore the UI without uploads.
  - Shows a notification on success; users can immediately run **Preprocess & map** and analyses.
- **AUTO preflight UI (Shiny):** Now detects cohort kind (`TN`, `ER/HER2`, `ER`-only, `HER2`-only) and shows compact stats:
  - ER/HER2 subgroups: **ER+/HER2-**, **ER-/HER2-**, **ER+/HER2+**, **ER-/HER2+**
  - TN cohorts: **TN** vs **nonTN**
  - Readiness uses the same minimums used by AUTO (sourced programmatically; no duplicated thresholds).
- **Shorter notifications.**
  - Routine toasts (e.g., "Step 1 complete. Proceed to Step 2.") now auto-dismiss sooner to reduce UI clutter.
- **Phenodata normalization (Mapping):** Accepts flexible ER/HER2/TN encodings and normalizes to canonical forms (`ER+/ER-`, `HER2+/HER2-`, `TN/nonTN`). Ambiguous `HER2="2+"` remains as-is and raises a warning.

### Bug fixes

- **TN cohorts + ssBC**: `BS_Multi()` now respects TN cohorts when methods are specified manually; `ssBC`/`ssBC.v2` switch to `s = "TN"` / `"TN.v2"` when a `TN` column indicates a TN cohort. Falls back to `s = "ER"` / `"ER.v2"` otherwise.
- **AUTO**: Fixed a crash in `BS_Multi(methods = "AUTO")` when ER and/or HER2 contained missing values (NA).
- **AUTO internals**: fixed variable name typo (`samples_ERHER2.icd`).
- **Mapping():** Robust ENTREZID coercion (from `as.character()` to `as.integer()` with suppressed warnings).
- **cIHC.itr**: `outList$distances` now returned as numeric matrix.

### Shiny

- **Surface method warnings as toasts:**
  - Runs are wrapped in a warning handler; package warnings (e.g., ssBC.v2 singleton fallbacks) appear as yellow notifications.
  - Warnings include subgroup, `n`, and example sample IDs for quick triage.
- **Shiny preflight reset**:
  - Fixed a stale cohort summary after switching data sources (manual uploads <-> example). The preflight panel now revalidates once inputs change.

### Developer notes

- Added lightweight internal logger `._msg()` and replaced scattered `message()` calls in AUTO to standardize package output without affecting CRAN/Bioc checks.

### Documentation

- **README/vignette:** brief note on the example-data button and expected file locations.
- **Mapping(): Column metadata clarified.** Added explicit requirements for receptor fields used by AUTO and ER/HER2/TN-dependent methods (`ssBC`, `cIHC`/`cIHC.itr`, `PCAPAM50`) and for ROR covariates (`TSIZE`, `NODE` as numeric 0/1). Documented preferred coding and automatic normalization behavior.

### Compatibility Notes

- SSPBC "full" output keeps previous columns for calls; additional metrics may appear.

## Upgrade Notes

- Raw RNA-seq counts are supported **from v1.1.3 onward** (requires gene lengths).
- If you previously parsed `BS` / `BS.Subtype`, switch to **`Call_5class` / `Call_4class`**.
- Package API unchanged.
