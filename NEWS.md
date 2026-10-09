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
