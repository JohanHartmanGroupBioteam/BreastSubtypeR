# BreastSubtypeR 1.4.2

## Bug fixes

* Preserve singleton gene rows under `Mapping(method = "max")` without
  ranking them (#178). In the nonraw `impute = FALSE` path, remaining `NA`
  and `NaN` measurements are retained instead of causing a dimension error.
  Selection among multiple probes and all other mapping methods are unchanged.

# BreastSubtypeR 1.4.1

Bug-fix release for Bioconductor 3.23. Each fix below was first made on the
development branch (1.5.2) and is cherry-picked here without changes to the
documented behaviour of correctly annotated cohorts; the AUTO routing,
PCAPAM50 and entropy changes of 1.5.2 are not included and ship with
Bioconductor 3.24.

## Bug fixes

- Fixed TNBC cohort detection in cohort-specific method selection (#133).
  Cohorts are now classified as TNBC only when all evaluable non-missing `TN`
  annotations indicate TN. Mixed cohorts containing both TN and nonTN samples
  are no longer routed to the TNBC-specific branch solely because some
  samples are TN. Added synthetic phenotype-table tests for TN-only, mixed,
  UNC-like mixed, and missing or partly missing TN annotations.
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
- `Mapping()` now accepts a `SummarizedExperiment` with a single sample.
  `duplicate_genes()`, `prepare_nc_matrix()` and the probe filter dropped the
  matrix dimensions of one-sample input and failed with "dim(X) must have a
  positive length". The collapsed gene-by-sample matrix is now built
  explicitly; results for multi-sample input are identical for all `method`
  values. Added regression tests.
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
- `BS_cIHC()`, `BS_cIHC.itr()` and `BS_PCAPAM50()` now stop with a clear
  message when the cohort lacks one of the two ER (or luminal/non-luminal IHC)
  groups needed for ER balancing, instead of failing with "undefined columns
  selected" or "arguments imply differing number of rows". Cohorts with both
  groups are unaffected. Added tests.
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

## Tests

- `tests/testthat.R` now calls `test_check("BreastSubtypeR")`, so the files
  under `tests/testthat/` run during `R CMD check`. Previously the file held a
  single inline test and `test-tn-cohort-detection.R` was never executed. The
  inline `BS_Multi()` test moved to `tests/testthat/test-bs-multi-manual.R`.
  testthat edition 3 is declared (`Config/testthat/edition`), and the Suggests
  entry requires testthat >= 3.2.0.

# BreastSubtypeR 1.3.2

## Highlights (from v1.1.3 onward)

-   Paper published in *NAR Genomics and Bioinformatics* (2025), **Editor’s Choice** (DOI: 10.1093/nargab/lqaf131).
-   Support for raw RNA-seq counts (requires gene lengths).
-   **iBreastSubtypeR** refresh: cleaner UX, smarter AUTO guidance, consistent exports.

### Enhancements

-   **ssBC/ssBC.v2: singleton subgroup robustness:** Subgroups with `n=1` no longer error:
    -   Keeps matrix shape (`drop=FALSE`) and hardens dimnames/types.
    -   Primary path: original `sspPredict()`.\
        Fallback: nearest-centroid (Spearman) when needed.
    -   If there are **0 common PAM50 genes**, returns `NA` labels with shaped `distances`/`dist.RORSubtype` to avoid downstream errors.
    -   ROR computation guarded for incomplete inputs.
-   **SSPBC output now "full":** `BS_sspbc()` and Shiny "sspbc" runs return a full metrics table (not calls-only).
    -   Exports map core label columns to the standard names (`Call_5class` / `Call_4class` when applicable).
-   **Shiny: "Load example data…" button**
    -   One-click load of a small demo dataset from `inst/RshinyTest/` to explore the UI without uploads.\
    -   Shows a notification on success; users can immediately run **Preprocess & map** and analyses.
-   **AUTO preflight UI (Shiny):** Now detects cohort kind (`TN`, `ER/HER2`, `ER`-only, `HER2`-only) and shows compact stats:
    -   ER/HER2 subgroups: **ER+/HER2-**, **ER-/HER2-**, **ER+/HER2+**, **ER-/HER2+**\
    -   TN cohorts: **TN** vs **nonTN**\
    -   Readiness uses the same minimums used by AUTO (sourced programmatically; no duplicated thresholds).
-   **Shorter notifications.**\
    -Routine toasts (e.g., "Step 1 complete. Proceed to Step 2.") now auto-dismiss sooner to reduce UI clutter.
-   **Phenodata normalization (Mapping):** Accepts flexible ER/HER2/TN encodings and normalizes to canonical forms (`ER+/ER-`, `HER2+/HER2-`, `TN/nonTN`). Ambiguous `HER2="2+"` remains as-is and raises a warning.

### Bug fixes

-   **TN cohorts + ssBC**: `BS_Multi()` now respects TN cohorts when methods are specified manually;
`ssBC`/`ssBC.v2` switch to `s = "TN"` / `"TN.v2"` when a `TN` column indicates a TN cohort.
Falls back to `s = "ER"` / `"ER.v2"` otherwise.
-   **AUTO**: Fixed a crash in BS_Multi(methods="AUTO") when ER and/or HER2 contained missing values (NA).
-   **AUTO internals**: fixed variable name typo (`samples_ERHER2.icd`).
-   **Mapping():** Robust ENTREZID coercion (from `as.character()` to `as.integer()` with suppressed warnings).
-   **cIHC.itr**: `outList$distances` now returned as numeric matrix.

### Shiny

-   **Surface method warnings as toasts:**
    -   Runs are wrapped in a warning handler; package warnings (e.g., ssBC.v2 singleton fallbacks) appear as yellow notifications.\
    -   Warnings include subgroup, `n`, and example sample IDs for quick triage.
-   **Shiny preflight reset**:
    -   Fixed a stale cohort summary after switching data sources (manual uploads ↔ example). The preflight panel now revalidates once inputs change.

### Developer notes

-   Added lightweight internal logger `._msg()` and replaced scattered `message()` calls in AUTO to standardize package output without affecting CRAN/Bioc checks.

### Documentation

-   **README/vignette:** brief note on the example-data button and expected file locations.
-   **Mapping(): Column metadata clarified.** Added explicit requirements for receptor fields used by AUTO and ER/HER2/TN-dependent methods (`ssBC`, `cIHC`/`cIHC.itr`, `PCAPAM50`) and for ROR covariates (`TSIZE`, `NODE` as numeric 0/1). Documented preferred coding and automatic normalization behavior.

### Compatibility Notes

-   SSPBC "full" output keeps previous columns for calls; additional metrics may appear.

## Upgrade Notes

-   Raw RNA-seq counts are supported **from v1.1.3 onward** (requires gene lengths).
-   If you previously parsed `BS` / `BS.Subtype`, switch to **`Call_5class` / `Call_4class`**.
-   Package API unchanged.
