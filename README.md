# BreastSubtypeR <a href='https://github.com/JohanHartmanGroupBioteam/BreastSubtypeR'><img src="inst/ShinyBreastSubtypeR/logo.svg" alt="BreastSubtypeR logo" align="right" height="110"/></a>

<!-- badges: start -->
[![Bioconductor Release](https://bioconductor.org/shields/years-in-bioc/BreastSubtypeR.svg)](https://bioconductor.org/packages/BreastSubtypeR)
[![Bioconductor Devel](https://bioconductor.org/shields/build/devel/bioc/BreastSubtypeR.svg)](https://bioconductor.org/packages/devel/bioc/html/BreastSubtypeR.html)
[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://opensource.org/licenses/GPL-3.0)
[![Paper: NAR Genomics & Bioinformatics](https://img.shields.io/badge/Paper-NAR%20Genomics%20%26%20Bioinformatics-0a7)](https://doi.org/10.1093/nargab/lqaf131)
<!-- badges: end -->

**BreastSubtypeR** is an assumption-aware, multi-method R/Bioconductor package with a local Shiny app.
It consolidates published intrinsic subtyping methods under one API and lets you run multiple classifiers at once.
**AUTO** inspects cohort diagnostics to select compatible methods and reduce misclassification.

> *Research use only; in clinical practice, intrinsic molecular subtyping is standardized via approved diagnostics (e.g., Prosigna®).*

📄 **Publication:** *NAR Genomics and Bioinformatics* (2025), **Editor’s Choice** → [doi:10.1093/nargab/lqaf131](https://doi.org/10.1093/nargab/lqaf131)

<details>
<summary><strong>How to cite</strong> (plain text &amp; BibTeX)</summary>

**Plain text**

Yang Q., Hartman J., Sifakis E.G. (2025). BreastSubtypeR: a unified R/Bioconductor package for intrinsic molecular subtyping in breast cancer research. *NAR Genomics and Bioinformatics*. <https://doi.org/10.1093/nargab/lqaf131>

**BibTeX**
```bibtex
@article{Yang2025BreastSubtypeR,
  author  = {Yang, Qiao and Hartman, Johan and Sifakis, Emmanouil G.},
  title   = {{BreastSubtypeR}: a unified R/Bioconductor package for intrinsic molecular subtyping in breast cancer research},
  journal = {NAR Genomics and Bioinformatics},
  year    = {2025},
  volume  = {7},
  number  = {4},
  pages   = {lqaf131},
  doi     = {10.1093/nargab/lqaf131},
  url     = {https://doi.org/10.1093/nargab/lqaf131}
}

```
</details>

------------------------------------------------------------------------

## Features

- **Unified interface for published methods:** consolidates PAM50 variants, AIMS, ssBC/sspbc, and others under one consistent API.
- **Run multiple methods at once (`BS_Multi`):** execute several classifiers in a single call and compare results side by side.
- **AUTO (cohort-aware selection):** checks ER/HER2 distribution, subtype purity, and subgroup sizes; disables incompatible classifiers.
- **Method-specific pre-processing:** routes raw RNA-seq counts, or log2-normalized RNA-seq/microarray/nCounter expression with `RawCounts = FALSE`.
- **Robust mapping:** Entrez ID–based gene mapping with conflict resolution.
- **Local Shiny app (`iBreastSubtypeR`):** point-and-click analysis; data stay on your machine.
- **Reproducibility:** Bioconductor distribution, unit tests, vignettes, and `SummarizedExperiment` compatibility.

------------------------------------------------------------------------

### Methods included (single-method implementations)

| **Method id** | **Short description** | **Group** | **Reference** |
|------------------|-------------------|------------------|------------------|
| `parker.original` | Original PAM50 by Parker et al., 2009 | NC-based | [Parker et al., 2009](https://doi.org/10.1200/JCO.2008.18.1370) |
| `genefu.scale` | PAM50 implementation as in the genefu R package (scaled version) | NC-based | [Gendoo et al., 2016](https://doi.org/10.1093/bioinformatics/btv693) |
| `genefu.robust` | PAM50 implementation as in the genefu R package (robust version) | NC-based | [Gendoo et al., 2016](https://doi.org/10.1093/bioinformatics/btv693) |
| `cIHC` | Conventional ER-balancing using immunohistochemistry (IHC) | NC-based | [Ciriello et al., 2015](https://doi.org/10.1016/j.cell.2015.09.033) |
| `cIHC.itr` | Iterative version of cIHC | NC-based | [Curtis et al., 2012](https://doi.org/10.1038/nature10983) |
| `PCAPAM50` | Selects IHC-defined ER subsets, then uses Principal Component Analysis (PCA) to create ESR1 expression-based ER-balancing | NC-based | [Raj-Kumar et al., 2019](https://doi.org/10.1038/s41598-019-44339-4) |
| `ssBC` | Subgroup-specific gene-centering PAM50 | NC-based | [Zhao et al., 2015](https://doi.org/10.1186/s13058-015-0520-4) |
| `ssBC.v2` | Updated subgroup-specific gene-centering PAM50 | NC-based | [Fernandez-Martinez et al., 2020](https://doi.org/10.1200/JCO.20.01276), [2024](https://doi.org/10.1001/jamaoncol.2023.7304) |
| `AIMS` | Absolute Intrinsic Molecular Subtyping (AIMS) method | SSP-based | [Paquet & Hallett, 2015](https://doi.org/10.1093/jnci/dju357) |
| `sspbc` | Single-Sample Predictors for Breast Cancer (AIMS adaptation) | SSP-based | [Staaf et al., 2022](https://doi.org/10.1038/s41523-022-00465-3) |

(See the vignette for implementation details.)

Zhao et al. (2015) introduced subgroup-specific gene centering.
Fernandez-Martinez et al. (2020) described its HER2/ER-specific application.
The packaged `ssBC.v2` quantiles match the public `SIGMA.txt` reference in the
[Perou Lab normalization repository](https://github.com/perou-lab/PAM50_ER_HER2_normalization),
associated with Fernandez-Martinez et al. (2024). Across the 50 PAM50 genes,
`ER.v2` uses `ERneg_HER2neg`, `ERpos_HER2neg`, `HER2pos_ERneg` and
`HER2pos_ERpos`; `TN.v2` uses `TNBC`. The original `ER` and `TN` quantiles
remain separate. See `?BreastSubtypeRobj` for the pinned numerical reference.

------------------------------------------------------------------------

## Installation

Install the package from the Bioconductor release compatible with your R installation:

``` r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("BreastSubtypeR")
```

Bioconductor release, Bioconductor devel, and GitHub are separate distribution
channels. GitHub may contain newer development changes that are not yet
available through Bioconductor. Consult the official
[release package page](https://bioconductor.org/packages/release/bioc/html/BreastSubtypeR.html)
and [devel package page](https://bioconductor.org/packages/devel/bioc/html/BreastSubtypeR.html)
for available versions and their R requirements. Use the R version compatible
with your selected Bioconductor release or devel environment; the package's
minimum R requirement does not replace that channel's requirements.
After installing from any channel, check your installed version with
`packageVersion("BreastSubtypeR")`.

For the development Bioconductor environment, use a separate R library and the
R version required by [Bioconductor devel](https://bioconductor.org/packages/devel/bioc/html/BreastSubtypeR.html):

``` r
BiocManager::install(version = "devel")
BiocManager::install("BreastSubtypeR")
```

Or install the development source from GitHub:

``` r
if (!require("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("JohanHartmanGroupBioteam/BreastSubtypeR")
packageVersion("BreastSubtypeR")
```

## Quick start

> *These examples use datasets shipped with the package. For your own data, provide a `SummarizedExperiment` **with clinical metadata in `colData`** (e.g., `PatientID`, `ER`, `HER2`; for ROR: `TSIZE`, `NODE`).*

``` r
library(BreastSubtypeR)

# Example data
data("BreastSubtypeRobj")
data("OSLO2EMIT0obj")
```

**1) Map & prepare (method-specific pre-processing + mapping)**

``` r
data_input <- Mapping(OSLO2EMIT0obj$se_obj, RawCounts = FALSE, method = "max", impute = TRUE)
```

**2) Multi-method run (user-defined)**

``` r
res <- BS_Multi(data_input = data_input, methods = c("parker.original","PCAPAM50","sspbc"))
head(res$res_subtypes, 5)
```

**3) AUTO mode (cohort-aware selection) + visualize**

``` r
res_auto <- BS_Multi(data_input = data_input, methods = "AUTO")
Vis_Multi(res_auto$res_subtypes)
```

***AUTO logic (clarifications)***

- **TNBC cohorts** take precedence when all evaluable TN annotations indicate
TN. At least 18 TN samples enable **ssBC and ssBC.v2** + SSP (AIMS, sspbc);
smaller TN cohorts run AIMS and sspbc only.
- A cohort is treated as HER2+ when it has at least one evaluable HER2 value
and every evaluable HER2 value is HER2+. The test reads the HER2 column
alone: one HER2- sample ends the HER2+ classification whether or not its ER
value is known, and a sample whose HER2 value is missing or equivocal ("2+")
neither creates nor removes a HER2+ cohort. AUTO reports how many samples
had no evaluable HER2 value and therefore took no part in the decision.
- In this **HER2+ cohort** route, **ssBC.v2** + SSP run when there are at least
8 ER+/HER2+ samples or 9 ER-/HER2+ samples; otherwise only SSP run. This route
does not additionally require the total-ER minima of 15 and 18.
- Cohorts without any evaluable HER2 value are not treated as HER2+; they
follow the ER-based rules (ssBC.v2 then returns NA for every sample).
- For the remaining cohorts, both ER totals below their minima run SSP only.
If only one ER group reaches its minimum, its HER2 subgroups determine whether
**ssBC and/or ssBC.v2** run alongside SSP.
- **Mixed cohorts** with at least 15 ER+ and 18 ER- samples run ten methods when
`0.39 < n_ERpos / (n_ERpos + n_ERneg) < 0.69`. The inequalities are strict;
at or outside either boundary, eight methods run, excluding `parker.original`
and `genefu.scale`. Both branches include `cIHC.itr`.
- Minimum sizes (defaults): ER+ = 15, ER− = 18, TN = 18. 
Subgroups use half of ER totals (rounded):
    -   ER+: `n_ERposHER2pos_threshold = n_ERposHER2neg_threshold = round(15 / 2)`
    -   ER−: `n_ERnegHER2pos_threshold = n_ERnegHER2neg_threshold = round(18 / 2)`

*Provenance*: ER+/ER− minimums are simulation-based defaults. 
TN currently uses the ER− minimum (18). 
Subgroup gates and TN minimum may be refined in future releases 
as additional simulations become available.

**4) Classify a single sample with AIMS**

``` r
one_sample <- Mapping(OSLO2EMIT0obj$se_obj[, 1], RawCounts = FALSE)
aims_one <- BS_AIMS(one_sample$se_SSP)
aims_one$cl
stopifnot(
    identical(dim(aims_one$cl), c(1L, 1L)),
    identical(rownames(aims_one$cl), colnames(one_sample$se_SSP)),
    identical(dim(aims_one$prob), c(1L, 1L)),
    identical(dim(aims_one$all.probs[["20"]]), c(1L, 5L))
)
```

AIMS assigns one of five subtypes independently to each sample. Its `cl` and
`prob` components retain sample-by-one matrix dimensions and sample IDs for
one or many samples. In `BS_Multi(Subtype = TRUE)`, the five-class summary is
retained and a separate four-class summary is added; AIMS Normal calls become
`NA` in that four-class summary because AIMS has no four-class model.

**5) Configure iterative ER balancing**

`BS_cIHC.itr()` defaults to `iteration = 100`, `ratio = 54 / 64`, and
`er_subsample = "target"`. The ratio is the target number of ER-positive
samples per ER-negative sample, not the ER-positive fraction of a cohort.
It must be one finite, strictly positive number; target mode accepts values
above 1. Write `54 / 64`, not `54:64`.

Target mode subsamples the group over-represented relative to the target,
using ceiling rounding, and retains the other group. For reproduction of
BreastSubtypeR 1.5.5 and earlier, `er_subsample = "larger"` instead subsamples
the numerically larger group; its ratio cannot exceed `n_larger / n_smaller`.
See `?BS_cIHC.itr` for the exact rules and output structures.

**6) Launch the local Shiny app**

``` r
BreastSubtypeR::iBreastSubtypeR() # interactive GUI (local)
```

***Notes:***

- The app runs locally; no data leave your machine.

- If you see a missing UI dependency:

``` r
install.packages(c("shiny","bslib"))
```

***Example data (for Shiny & scripts)***

You can explore the app without preparing files:

- In **iBreastSubtypeR**, click **“Load example data…”** to pre-fill
  the UI with a small demo cohort. Then click **Preprocess & map**, choose a
  method and settings in Step 2, and run subtyping.

- For uploaded data, supply expression, clinical and annotation tables with
  matching sample/feature identifiers. Select normalized **log2** expression
  or raw RNA-seq counts; raw counts also require gene lengths in the annotation.

- Reloading the example, replacing any input file, or changing the input data
  type clears the mapped data, previous predictions, plots and downloads.
  Preprocess and run again. Changing relevant prediction settings clears predictions,
  plots and downloads while retaining valid mapped data when remapping is
  unnecessary.

- Downloads appear only after a successful current prediction. Calls-only
  exports are available for every method; Full metrics is available for
  supported single NC methods, not AUTO, AIMS or sspbc. Old download URLs are
  also rejected by the server after invalidation.

- NC Full metrics include ROR-S and ROR-P. Enable **Use clinical variables
  (ROR)** with valid numeric `TSIZE` and `NODE` to add clinical ROR-C/ROR-PC
  outputs. `TSIZE` is binary (0 for at most 2 cm, 1 for more than 2 cm).
  For `NODE`, supply the actual nonnegative integer number of positive nodes
  when available. Binary 0/1 is accepted but cannot distinguish 1--3 from
  4 or more nodes for `ROR-PC Group (Subtype + Clinic + Prolif.Subtype)`.
  NODE affects only that additional risk group, not numerical ROR scores.
  The R API does not enforce nonnegative integer counts; see `?Mapping`
  for missing-value behavior. Exact-cutoff behavior is tracked separately
  in [issue #154](https://github.com/JohanHartmanGroupBioteam/BreastSubtypeR/issues/154).
  AUTO provides calls and entropy rather than a Full metrics export.

- Shiny cIHC.itr uses target-mode balancing and the ratio semantics above.
  Iterations default to 100, with a minimum of 10 and a UI step of 10; other
  whole numbers such as 11 are allowed. Invalid ratios or iterations show an
  actionable error and leave the session connected. Correct the setting and
  rerun using the valid mapped input; failed runs cannot export old results.

- Programmatically, the same files are shipped inside the package:
```r
  exdir   <- system.file("RshinyTest", package = "BreastSubtypeR")
  gex     <- file.path(exdir, "OSLO2EMIT0_GEX_log2.FPKM.txt")
  clin    <- file.path(exdir, "OSLO2EMIT0_clinical.txt")
  anno    <- file.path(exdir, "OSLO2EMIT0_anno.txt")
  stopifnot(all(file.exists(gex, clin, anno)))
```

## Vignette & documentation

A comprehensive usage guide (input types, AUTO details, full method descriptions) is included as a **vignette**.

See function help pages for specifics (e.g., `?BS_Multi`, `?Mapping`, `?iBreastSubtypeR`).

## Contributing & issues

**Canonical source**  
Bioconductor package page: <https://bioconductor.org/packages/BreastSubtypeR>  
Bioconductor DOI: <https://doi.org/10.18129/B9.bioc.BreastSubtypeR>  
Development repository: <https://github.com/JohanHartmanGroupBioteam/BreastSubtypeR>

**Support & bugs**  
Bugs/PRs: <https://github.com/JohanHartmanGroupBioteam/BreastSubtypeR/issues>

## License

GPL-3
