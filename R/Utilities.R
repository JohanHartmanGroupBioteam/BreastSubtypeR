#'
#' Functions for BreastSubtypeR package
#' @keywords internal
#' @import ComplexHeatmap
#' @import RColorBrewer
#' @import circlize
#' @import ggrepel
#' @import magrittr
#' @import impute
#' @importFrom dplyr mutate across everything
#' @importFrom ggplot2 ggplot aes geom_point geom_bar geom_text geom_line geom_vline geom_hline
#' @importFrom ggplot2 geom_boxplot scale_x_continuous scale_color_manual scale_fill_manual
#' @importFrom ggplot2 coord_polar labs theme theme_classic theme_minimal theme_void element_text
#' @importFrom SummarizedExperiment SummarizedExperiment colData rowData assay
#'
#' @noRd
NULL
NULL

#' Internal message helper
#' @noRd
.msg <- function(fmt, ..., origin = c("AUTO", "MANUAL")) {
    origin <- match.arg(origin)
    if (isTRUE(getOption("BreastSubtypeR.verbose", TRUE))) {
        message(sprintf("[%s] %s", origin, sprintf(fmt, ...)))
    }
}

#' Impute missing values using k-nearest neighbors
#'
#' @param x Numeric matrix with missing values (NA).
#' @param verbose Logical; if TRUE, prints missing probes/samples.
#' @return Imputed matrix.
#' @noRd
impute_missing <- function(x, verbose = FALSE) {
    if (anyNA(x)) {
        if (verbose) {
            probeid_NA <- rownames(x)[rowSums(is.na(x)) > 0]
            sample_NA <- colnames(x)[colSums(is.na(x)) > 0]
            message("Missing probes: ", paste(probeid_NA, collapse = ", "))
            message("Missing samples: ", paste(sample_NA, collapse = ", "))
        }
        x <- impute::impute.knn(x)$data
    }
    x
}


#' Collapse duplicate genes by a summary statistic
#'
#' @param x Numeric matrix (rows = probes/genes, cols = samples).
#' @param y Dataframe with `probe` and `ENTREZID` columns.
#' @param method Method for resolving duplicates ("max", "mean", etc.).
#' @return Deduplicated matrix (rows = unique Entrez IDs).
#' @noRd
duplicate_genes <- function(x, y, method) {
    probeid <- rownames(x)
    entrezid <- y$ENTREZID
    names(entrezid) <- y$probe

    ## process probeid in input data
    entrezid <- entrezid[probeid]
    ## remove NA
    entrezid <- entrezid[!(is.na(entrezid))]
    x <- x[names(entrezid), , drop = FALSE]
    entrezid <- factor(entrezid, levels = unique(entrezid))
    ## names are unique probeid and content are redundant entrezid

    ## This is for probeID or transcriptID
    ## split expression matrix
    split_mat <- split.data.frame(as.data.frame(x), entrezid, drop = FALSE)

    # function to calculate the desired statistic
    calculate_stat <- function(mat, method) {
        switch(method,
            "mean" = apply(mat, 2, mean, na.rm = TRUE),
            "median" = apply(mat, 2, median, na.rm = TRUE),
            "max" = mat[which.max(rowSums(mat)), ],
            "stdev" = {
                temp <- mat
                stdevs <- apply(mat, 1, sd, na.rm = TRUE)
                vals <- temp[match(max(stdevs), stdevs), ]
                return(vals)
            },
            "iqr" = {
                temp <- mat
                iqrs <- apply(mat, 1, function(x) {
                    quantile(x, 0.75, na.rm = TRUE) -
                        quantile(x, 0.25, na.rm = TRUE)
                })
                vals <- temp[match(max(iqrs), iqrs), ]
                return(vals)
            }
        )
    }

    ## one row per Entrez ID, one column per sample; built explicitly so that
    ## single-sample matrices keep their dimensions
    collapsed <- lapply(split_mat, function(mat) {
        as.numeric(unlist(calculate_stat(mat, method), use.names = FALSE))
    })
    out <- do.call(rbind, collapsed)
    dimnames(out) <- list(names(split_mat), colnames(x))

    return(out)
}


#' Collapse duplicate genes by a summary statistic
#'
#' @param x Numeric matrix (rows = probes/genes, cols = samples).
#' @param y Dataframe with `probe` and `ENTREZID` columns.
#' @return Deduplicated matrix (rows = unique Entrez IDs).
#' @noRd
prepare_nc_matrix <- function(x, genes.sig50, samplenames, verbose) {
    ## print necessary information
    ## Parker
    missing_ID_parker <- setdiff(genes.sig50$EntrezGene.ID, rownames(x))

    if (verbose && length(missing_ID_parker) > 0L) {
        message("Absent PAM50 reference Entrez IDs: ",
            paste(missing_ID_parker, collapse = ", "))
    }

    ## get matrix for NC (symbol as rows, sample as col)
    genes_nc <- genes.sig50$EntrezGene.ID
    x_NC <- x[na.omit(match(genes_nc, rownames(x))), , drop = FALSE]
    rownames(x_NC) <- genes.sig50$Symbol[match(rownames(x_NC), genes_nc)]
    x_NC <- data.frame(x_NC)
    colnames(x_NC) <- samplenames

    return(x_NC)
}


#' Prepare matrix for SSP-based methods
#'
#' @param x Numeric matrix (rows = Entrez IDs, cols = samples).
#' @param genes.s Dataframe with gene signatures.
#' @param RawCounts Logical; if TRUE, input is raw counts.
#' @param verbose Logical; print missing genes.
#' @return Processed matrix for SSP.
#' @noRd
prepare_ssp_matrix <- function(x, genes.s, RawCounts, samplenames, verbose) {
    ## AIMS
    genes_ssp <- genes.s$EntrezGene.ID[genes.s$SSP_based == "Yes"]
    missing <- setdiff(genes_ssp, rownames(x))

    if (verbose && length(missing) > 0) {
        message("Absent Entrez IDs in the shared SSP input: ",
            paste(missing, collapse = ", "))
    }

    ## get matrix for AIMS (entrezID as colnames)
    x_SSP <- x[intersect(genes_ssp, rownames(x)), , drop = FALSE]
    colnames(x_SSP) <- samplenames

    ## exponential transformation for microarray/RNAseq
    if (!RawCounts) {
        x_SSP <- 2^x_SSP
    }

    return(x_SSP)
}

#' Helper to normalize ER/HER2/TN columns in Column metadata 
#'
#' @keywords internal
#' @noRd
.normalize_er_her2_tn <- function(df, quiet = TRUE) {
  canon <- function(x) trimws(toupper(as.character(x)))
  map_bin <- function(x, pos_label, neg_label) {
    x0  <- canon(x)
    pos <- c("POSITIVE","POS","1","TRUE","T","YES","Y",
             gsub("\\+","POSITIVE", pos_label), pos_label)
    neg <- c("NEGATIVE","NEG","0","FALSE","F","NO","N",
             gsub("-","NEGATIVE", neg_label), neg_label)
    ## as.character(): with a factor column the ifelse() fallback would
    ## otherwise return the integer codes of unmatched levels
    out <- ifelse(x0 %in% pos, pos_label,
                  ifelse(x0 %in% neg, neg_label, as.character(x)))
    out
  }
  
  if ("ER" %in% names(df)) {
    old <- as.character(df$ER)
    df$ER <- map_bin(old, "ER+", "ER-")
    if (!identical(old, df$ER))
      warning("Phenodata: coerced ER values to {ER+, ER-}.", call. = FALSE)
  }
  if ("HER2" %in% names(df)) {
    old <- as.character(df$HER2)
    has2p <- grepl("\\b2\\+\\b", old, ignore.case = TRUE)
    df$HER2 <- ifelse(has2p, old, map_bin(old, "HER2+", "HER2-"))
    if (!identical(old, df$HER2))
      warning("Phenodata: coerced HER2 values to {HER2+, HER2-} (skipped '2+').", call. = FALSE)
  }
  if ("TN" %in% names(df)) {
    old <- as.character(df$TN)
    x0  <- canon(old)
    pos <- c("TRUE","T","YES","Y","1","TN","TNBC")
    neg <- c("FALSE","F","NO","N","0","NON-TN","NONTN","NON_TN")
    df$TN <- ifelse(x0 %in% pos, "TN",
                    ifelse(x0 %in% neg, "nonTN", old))
    if (!identical(old, df$TN))
      warning("Phenodata: coerced TN values to {TN, nonTN}.", call. = FALSE)
  }
  note <- function(msg) if (!quiet) .msg(msg)
  df
}

#' Summarize TN annotations for cohort-level TN-only detection
#'
#' @keywords internal
#' @noRd
.tn_only_summary <- function(tn_flag) {
    tn_chr <- toupper(trimws(as.character(tn_flag)))
    tn_evaluable <- !is.na(tn_chr) & nzchar(tn_chr)
    tn_truthy <- tn_chr %in% c("TN", "TRUE", "T", "YES", "Y", "1")

    n_tn <- sum(tn_truthy & tn_evaluable)
    n_evaluable_tn <- sum(tn_evaluable)

    list(
        n_tn = n_tn,
        n_evaluable_tn = n_evaluable_tn,
        is_tn_only = n_evaluable_tn > 0L && n_tn == n_evaluable_tn
    )
}


#' Map Gene IDs and Handle missing data
#'
#' @param method A string specifying the method for resolving duplicate probes
#' in microarray or RNA-seq data (see `duplicate_genes()`; "mean", "median",
#' "iqr" and "stdev" follow collapseIDs() of the original PAM50 code):
#'   - `"mean"`: per-sample mean across the duplicate probes.
#'   - `"median"`: per-sample median across the duplicate probes.
#'   - `"iqr"`: keeps the probe with the highest interquartile range (IQR),
#'   typically used for short-oligo arrays (e.g., Affymetrix).
#'   - `"stdev"`: keeps the probe with the highest standard deviation.
#'   - `"max"`: keeps the probe with the largest row sum across samples,
#'   often used for RNA-seq data.
#' @noRd

domapping <- function(
        se_obj,
        RawCounts = FALSE,
        method = "max",
        impute = TRUE,
        verbose = TRUE) {
    ## 1. Input raw counts
    if (RawCounts && !"Length" %in% colnames(rowData(se_obj))) {
        stop(
            "Missing gene length information: for RawCounts=TRUE, ",
            "rowData(se_obj)$Length must be provided (in base pairs) ",
            "to enable FPKM calculation."
        )
    }

    # 2. Load gene signatures
    data_env <- new.env(parent = emptyenv())
    data("BreastSubtypeRobj", envir = data_env, package = "BreastSubtypeR")
    BreastSubtypeRobj <- data_env[["BreastSubtypeRobj"]]
    genes.s <- BreastSubtypeRobj$genes.signature

    # 3. Normalize raw counts (if applicable)
    if (RawCounts) {
        # Pull counts from SummarizedExperiment and build a DGEList
        counts <- SummarizedExperiment::assay(se_obj)
        y_all <- edgeR::DGEList(counts = round(as.matrix(counts)))

        # Light positive-only subset (very gentle): CPM > 0 in ≥30% samples
        prop_pos <- rowMeans(edgeR::cpm(y_all) > 0.1)
        keep_sub <- prop_pos >= 0.30
        if (!any(keep_sub)) stop("Cannot form a positive subset for UQ normalization.")
        y_sub <- y_all[keep_sub, , keep.lib.sizes = TRUE]

        # Estimate UQ factors on the subset (p = 0.75)
        y_sub <- edgeR::calcNormFactors(y_sub, method = "upperquartile", p = 0.75)

        # Apply those factors to the full data
        y_all$samples$norm.factors <- y_sub$samples$norm.factors

        # NC-based input for downstream: log-CPM on the FULL matrix (uses effective lib sizes)
        x <- edgeR::cpm(y_all, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)

        # FPKM from raw counts + gene lengths (in base pairs)
        gl <- as.numeric(SummarizedExperiment::rowData(se_obj)$Length)
        names(gl) <- rownames(se_obj)
        gl <- gl[rownames(y_all)]
        counts.fpkm <- edgeR::rpkm(y_all,
            gene.length = gl,
            normalized.lib.sizes = FALSE, log = FALSE
        )
    } else {
        x <- SummarizedExperiment::assay(se_obj)
    }

    # 4. Check feature data
    ## Extract data from SummarizedExperiment
    y <- rowData(se_obj) %>% data.frame()
    ## check feature data
    if (length(y) == 0) {
        stop("Please provide feature annotation to do probeID mapping ")
    }
    ## change data type
    y$ENTREZID <- suppressWarnings(as.integer(as.character(y$ENTREZID)))
    n_bad <- sum(is.na(y$ENTREZID))
    if (verbose && n_bad > 0L) {
        message(sprintf("Feature rows with non-integer ENTREZID dropped before mapping: %d", n_bad))
    }
    samplenames <- colnames(x)

    # 5. Filter by signature genes and impute
    ## filter by ENTREZID
    y <- y[y$ENTREZID %in% genes.s$EntrezGene.ID, ]
    x <- x[y$probe, , drop = FALSE]
    if (impute && anyNA(x)) x <- impute_missing(x, verbose)
    if (RawCounts && impute && anyNA(counts.fpkm)) {
        counts.fpkm <- impute_missing(counts.fpkm, verbose)
    }

    # 6. Deduplicate by EntrezID
    x <- duplicate_genes(x, y, method)
    if (RawCounts) {
        counts.fpkm <- duplicate_genes(counts.fpkm, y, method)
    }

    # 7. Prepare output (NC and SSP matrices)
    x_NC <- prepare_nc_matrix(x, BreastSubtypeRobj$genes.sig50, samplenames, verbose)
    x_SSP <- prepare_ssp_matrix(
        if (RawCounts) counts.fpkm else x,
        genes.s,
        RawCounts,
        samplenames,
        verbose
    )

    result <- list(x_NC = x_NC, x_SSP = x_SSP)

    return(result)
}

#' Function implementing AUTO decision logic
#' @noRd
get_methods <- function(pheno) {
    #### AUTO mode
    cohort.select <- "ERpos"
    samples_ER.icd <- NULL
    samples_ERHER2.icd <- NULL
    methods <- NULL

    if (ncol(pheno) == 0) {
        .msg("The pheno table has not been detected.", origin = "AUTO")
        .msg("Running methods: AIMS, sspbc", origin = "AUTO")
        methods <- c("AIMS", "sspbc")
    } else if (!all(c("ER", "HER2") %in% colnames(pheno))) {
        stop("The 'AUTO' mode requires both 'ER' and 'HER2' columns in the 'pheno' dataframe.")
    } else {
      
        pheno <- .normalize_er_her2_tn(pheno)
        
        # Detect missing receptor annotations (optional)
        n_na_ER   <- sum(is.na(pheno$ER))
        n_na_HER2 <- sum(is.na(pheno$HER2))
        if (n_na_ER > 0L || n_na_HER2 > 0L) {
          .msg("AUTO: missing ER/HER2 values detected (ER NA=%d, HER2 NA=%d). These samples are excluded from subgroup counts.",
               n_na_ER, n_na_HER2, origin = "AUTO")
        }
      
        # Calculate sample sizes
        sample_counts <- with(pheno, table(ER, HER2)) # (kept for future debugging)
        n_ERpos <- sum(pheno$ER == "ER+", na.rm = TRUE)
        n_ERneg <- sum(pheno$ER == "ER-", na.rm = TRUE)
        
        n_ERnegHER2pos <- sum(pheno$ER == "ER-" & pheno$HER2 == "HER2+", na.rm = TRUE)
        n_ERnegHER2neg <- sum(pheno$ER == "ER-" & pheno$HER2 == "HER2-", na.rm = TRUE)
        n_ERposHER2pos <- sum(pheno$ER == "ER+" & pheno$HER2 == "HER2+", na.rm = TRUE)
        n_ERposHER2neg <- sum(pheno$ER == "ER+" & pheno$HER2 == "HER2-", na.rm = TRUE)

        # Evaluable HER2 values (AUTO HER2+ cohort detection). Values
        # outside {HER2+, HER2-} are not evaluable: missing values and
        # equivocal codes such as "2+", which .normalize_er_her2_tn() leaves
        # unchanged, both fall in this group.
        n_HER2_known <- sum(pheno$HER2 %in% c("HER2+", "HER2-"))
        n_HER2pos_known <- sum(pheno$HER2 %in% "HER2+")
        n_HER2_unevaluable <- nrow(pheno) - n_HER2_known
        
        # Set thresholds
        n_ERpos_threshold <- 15 # simulation-based cut-off
        n_ERneg_threshold <- 18 # simulation-based cut-off
        n_TN_threshold <- 18 # same as ER- cut-off (initial heuristic; easy to update)
        n_ERposHER2pos_threshold <- round(n_ERpos_threshold / 2) # initial heuristics
        n_ERposHER2neg_threshold <- round(n_ERpos_threshold / 2) # initial heuristics
        n_ERnegHER2pos_threshold <- round(n_ERneg_threshold / 2) # initial heuristics
        n_ERnegHER2neg_threshold <- round(n_ERneg_threshold / 2) # initial heuristics
        upper_ratio <- 0.69 # simulation-based cut-off
        lower_ratio <- 0.39 # simulation-based cut-off

        ## ---- TNBC handling ------------------------------------------------
        if ("TN" %in% colnames(pheno)) {
            # Use TNBC cohort handling only when all evaluable TN annotations are TN.
            tn_summary <- .tn_only_summary(pheno$TN)
            n_TN <- tn_summary$n_tn

            if (tn_summary$is_tn_only) {
                cohort.select <- "TNBC"

                if (n_TN >= n_TN_threshold) {
                    .msg("A TNBC cohort has been detected.", origin = "AUTO")
                    .msg("Running methods: ssBC (TN), ssBC.v2 (TN), AIMS, sspbc", origin = "AUTO")
                    methods <- c("ssBC", "ssBC.v2", "AIMS", "sspbc")
                } else {
                    .msg("A small number of TNBC samples has been detected.", origin = "AUTO")
                    .msg("Running methods: AIMS, sspbc", origin = "AUTO")
                    methods <- c("AIMS", "sspbc")
                }
            } else if (tn_summary$n_evaluable_tn > 0L) {
                .msg("TN column detected but not all evaluable samples are TN;
              proceeding with non-TNBC logic.", origin = "AUTO")
            } else {
                .msg("TN column detected but no evaluable samples flagged as TN;
              proceeding with non-TNBC logic.", origin = "AUTO")
            }
        }

        ## ---- main panel (non-TNBC) ---------------------------------------
        if (is.null(methods)) { # only if TNBC branch did not set methods

            if (n_HER2_known > 0L && n_HER2pos_known == n_HER2_known) {
                ## HER2+ cohort: at least one evaluable HER2 value, and every
                ## evaluable HER2 value is HER2+. The test reads the HER2
                ## column alone. A HER2- sample therefore ends the HER2+
                ## classification whether or not its ER value is known, and a
                ## sample without an evaluable HER2 value neither creates nor
                ## removes a HER2+ cohort. Cohorts without any evaluable HER2
                ## value fall through to the ER-based rules below.
                .msg("A HER2+ cohort has been detected.", origin = "AUTO")
                if (n_HER2_unevaluable > 0L) {
                    .msg(
                        "%d of %d samples have no evaluable HER2 value and took no part in this decision.",
                        n_HER2_unevaluable, nrow(pheno),
                        origin = "AUTO"
                    )
                }
                cohort.select <- "HER2pos"

                if (n_ERposHER2pos < n_ERposHER2pos_threshold &&
                    n_ERnegHER2pos < n_ERnegHER2pos_threshold) {
                    .msg("A small HER2+ cohort has been detected.", origin = "AUTO")
                    .msg("Running methods: AIMS, sspbc", origin = "AUTO")
                    methods <- c("AIMS", "sspbc")
                } else if (n_ERposHER2pos >= n_ERposHER2pos_threshold &&
                    n_ERnegHER2pos < n_ERnegHER2pos_threshold) {
                    .msg("A ER+/HER2+ cohort has been detected.", origin = "AUTO")
                    .msg("Running methods for HER2+ samples: ssBC.v2, AIMS, sspbc", origin = "AUTO")
                    methods <- c("ssBC.v2", "AIMS", "sspbc")
                } else if (n_ERposHER2pos < n_ERposHER2pos_threshold &&
                    n_ERnegHER2pos >= n_ERnegHER2pos_threshold) {
                    .msg("A ER-/HER2+ cohort has been detected.", origin = "AUTO")
                    .msg("Running methods for HER2+ samples: ssBC.v2, AIMS, sspbc", origin = "AUTO")
                    methods <- c("ssBC.v2", "AIMS", "sspbc")
                } else {
                    .msg("Running methods for HER2+ samples: ssBC.v2, AIMS, sspbc", origin = "AUTO")
                    methods <- c("ssBC.v2", "AIMS", "sspbc")
                }
            } else if (n_ERpos < n_ERpos_threshold && n_ERneg < n_ERneg_threshold) {
                .msg("A small number of ER-/ER+ samples has been detected.", origin = "AUTO")
                .msg("Running methods: AIMS, sspbc", origin = "AUTO")
                methods <- c("AIMS", "sspbc")
            } else if (n_ERpos >= n_ERpos_threshold && n_ERneg < n_ERneg_threshold) {
                cohort.select <- "ERpos"

                if (n_ERposHER2pos >= n_ERposHER2pos_threshold &&
                    n_ERposHER2neg >= n_ERposHER2neg_threshold) {
                    .msg("Running methods for ER+ samples:\n
                  ssBC, ssBC.v2, AIMS, sspbc", origin = "AUTO")
                    methods <- c("ssBC", "ssBC.v2", "AIMS", "sspbc")
                } else if (n_ERposHER2pos >= n_ERposHER2pos_threshold &&
                    n_ERposHER2neg < n_ERposHER2neg_threshold) {
                    .msg("Running methods for ER+/HER2+ samples:\n
                  ssBC.v2, AIMS, sspbc", origin = "AUTO")
                    methods <- c("ssBC.v2", "AIMS", "sspbc")
                } else if (n_ERposHER2pos < n_ERposHER2pos_threshold &&
                    n_ERposHER2neg >= n_ERposHER2neg_threshold) {
                    .msg("Running methods for ER+/HER2- samples:\n
                  ssBC.v2, AIMS, sspbc", origin = "AUTO")
                    methods <- c("ssBC.v2", "AIMS", "sspbc")
                }
            } else if (n_ERpos < n_ERpos_threshold && n_ERneg >= n_ERneg_threshold) {
                cohort.select <- "ERneg"

                if (n_ERnegHER2pos >= n_ERnegHER2pos_threshold &&
                    n_ERnegHER2neg >= n_ERnegHER2neg_threshold) {
                    .msg("Running methods for ER- samples:\n
                  ssBC, ssBC.v2, AIMS, sspbc", origin = "AUTO")
                    methods <- c("ssBC", "ssBC.v2", "AIMS", "sspbc")
                } else if (n_ERnegHER2pos >= n_ERnegHER2pos_threshold &&
                    n_ERnegHER2neg < n_ERnegHER2neg_threshold) {
                    .msg("Running methods for ER-/HER2+ samples:\n
                  ssBC.v2, AIMS, sspbc", origin = "AUTO")
                    methods <- c("ssBC.v2", "AIMS", "sspbc")
                } else if (n_ERnegHER2pos < n_ERnegHER2pos_threshold &&
                    n_ERnegHER2neg >= n_ERnegHER2neg_threshold) {
                    .msg("Running methods for ER-/HER2- samples:\n
                  ssBC.v2, AIMS, sspbc", origin = "AUTO")
                    methods <- c("ssBC.v2", "AIMS", "sspbc")
                }
            } else if (n_ERpos >= n_ERpos_threshold && n_ERneg >= n_ERneg_threshold) {
                cohort.select <- "mixed"
                ratio_ER <- n_ERpos / (n_ERneg + n_ERpos)

                if (ratio_ER > lower_ratio && ratio_ER < upper_ratio) {
                    .msg("Running methods:\n
                  parker.original, genefu.scale, genefu.robust, ssBC, ssBC.v2,
                  cIHC, cIHC.itr, PCAPAM50, AIMS, sspbc", origin = "AUTO")
                    methods <- c(
                        "parker.original", "genefu.scale", "genefu.robust",
                        "ssBC", "ssBC.v2", "cIHC", "cIHC.itr", "PCAPAM50", "AIMS", "sspbc"
                    )
                } else {
                    .msg("The ER+ ratio in the current dataset differs from that observed in the UNC232 training cohort.", origin = "AUTO")
                    .msg("Running methods:\n
                  genefu.robust, ssBC, ssBC.v2, cIHC, cIHC.itr, PCAPAM50, AIMS, sspbc", origin = "AUTO")
                    methods <- c(
                        "genefu.robust", "ssBC", "ssBC.v2", "cIHC", "cIHC.itr",
                        "PCAPAM50", "AIMS", "sspbc"
                    )
                }
            }
        }

        ## ---- sample subsetting for ssBC variants (skip TNBC) --------------
        if (cohort.select != "TNBC") {
            # counts (unchanged)
            ERHER2_counts <- c(
                ERpos        = n_ERpos,
                ERneg        = n_ERneg,
                ERposHER2pos = n_ERposHER2pos,
                ERposHER2neg = n_ERposHER2neg,
                ERnegHER2pos = n_ERnegHER2pos,
                ERnegHER2neg = n_ERnegHER2neg
            )

            ## >= : the same minimums that select the methods above also
            ## admit a subgroup to the ssBC / ssBC.v2 sample subsets
            er_idx <- c(
                ERpos = ERHER2_counts["ERpos"] >= n_ERpos_threshold,
                ERneg = ERHER2_counts["ERneg"] >= n_ERneg_threshold
            )
            samples_ER <- names(ERHER2_counts)[seq(1, 2)][er_idx]

            erher2_idx <- c(
                ERposHER2pos = ERHER2_counts["ERposHER2pos"] >= n_ERposHER2pos_threshold,
                ERposHER2neg = ERHER2_counts["ERposHER2neg"] >= n_ERposHER2neg_threshold,
                ERnegHER2pos = ERHER2_counts["ERnegHER2pos"] >= n_ERnegHER2pos_threshold,
                ERnegHER2neg = ERHER2_counts["ERnegHER2neg"] >= n_ERnegHER2neg_threshold
            )
            samples_ERHER2 <- names(ERHER2_counts)[seq(3, 6)][erher2_idx]

            ## ssBC indices/messages only if ssBC was actually selected,
            ## and never in pure HER2+ cohorts
            if ("ssBC" %in% methods && cohort.select != "HER2pos") {
                if (length(samples_ER) > 0) {
                    .msg("ssBC for samples: %s", paste(samples_ER, collapse = ", "), origin = "AUTO")
                    samples_ER.icd <- unlist(lapply(samples_ER, function(subtype) {
                        subtype <- stringr::str_replace_all(subtype, "pos", "+")
                        subtype <- stringr::str_replace_all(subtype, "neg", "-")
                        ## which() drops samples with missing ER, which would
                        ## otherwise contribute NA sample names
                        rownames(pheno)[which(pheno$ER == subtype)]
                    }))
                }
            }

            ## ssBC.v2 indices/messages only if ssBC.v2 was selected
            if ("ssBC.v2" %in% methods) {
                if (length(samples_ERHER2) > 0) {
                    .msg("ssBC.v2 for samples: %s", paste(samples_ERHER2, collapse = ", "), origin = "AUTO")
                    samples_ERHER2.icd <- unlist(lapply(samples_ERHER2, function(subtype) {
                        subtype <- subtype |>
                            stringr::str_replace_all("pos", "+") |>
                            stringr::str_replace_all("neg", "-")
                        ER_sts <- substr(subtype, 1, 3)
                        HER2_sts <- substr(subtype, 4, 8)
                        rownames(pheno)[which(pheno$ER == ER_sts & pheno$HER2 == HER2_sts)]
                    }))
                }
            }
        }
    }

    if (length(samples_ER.icd) == 0) samples_ER.icd <- NULL
    if (length(samples_ERHER2.icd) == 0) samples_ERHER2.icd <- NULL
    if (is.null(methods)) {
        ## safety net: reached by ER+-only or ER--only cohorts whose HER2
        ## subgroups are all below their thresholds (for example when many
        ## HER2 values are missing); say so instead of falling back silently
        .msg(
            "No cohort rule matched the ER/HER2 subgroup sizes; running the single-sample predictors AIMS and sspbc only.",
            origin = "AUTO"
        )
        methods <- c("AIMS", "sspbc")
    }

    list(
        samples_ER.icd = samples_ER.icd,
        samples_ERHER2.icd = samples_ERHER2.icd,
        methods = methods,
        cohort.select = cohort.select
    )
}


#' Function for consensus subtype
#' @noRd
get_consensus_subtype <- function(patient_row) {
    patient_row <- unlist(patient_row, use.names = FALSE)
    counts <- table(patient_row)
    max_subtype <- names(counts)[which.max(counts)]
    return(max_subtype)
}

#' Function for entropy calculation
#'
#' Raw Shannon entropy in bits over the calls that are actually present in the
#' row. `table()` drops missing values, so a row in which every executed method
#' returned `NA` used to give `-sum(numeric(0))`, that is 0 - the same value a
#' row gets when every method agrees. Such a row has no call distribution and
#' therefore no entropy, and is reported as `NA` instead. Rows with at least
#' one call keep exactly the value computed before.
#'
#' @param patient_row Character vector of per-method calls for one sample.
#' @return Raw Shannon entropy in bits, or `NA_real_` when no call is present.
#' @noRd
get_entropy <- function(patient_row) {
    calls <- patient_row[!is.na(patient_row)]
    if (length(calls) == 0L) {
        return(NA_real_)
    }
    freq <- table(calls)
    prob <- freq / sum(freq)
    entropy <- -sum(prob * log2(prob))
    return(entropy)
}


#' Function to get the average correlation and ROR
#' @noRd
get_average_subtype <- function(res_ihc_iterative, consensus_subtypes) {
    ## correlation and ROR to be averaged
    sum_colnames <- c("Basal", "Her2", "LumA", "LumB", "Normal")

    all_patients <- names(res_ihc_iterative[[1]]$predictions)

    sum_cols_list <- mapply(
        function(res_ihc) {
            # res_ihc = as.data.frame( res_ihc_iterative[[1]])

            res_ihc$distances <- as.data.frame(res_ihc$distances)

            ## if FALSE, make the cell as NULL
            keep <- res_ihc$predictions == consensus_subtypes
            res_ihc$distances[!keep, ] <- as.list(rep(NA, ncol(res_ihc$distances)))

            ## convert every column to numeric (preserving NA), same as before
            res <- dplyr::mutate(
                res_ihc$distances,
                dplyr::across(
                    dplyr::everything(),
                    ~ ifelse(!is.na(.x), as.numeric(as.character(.x)), NA_real_)
                )
            )

            return(res)
        },
        res_ihc_iterative,
        SIMPLIFY = FALSE,
        USE.NAMES = FALSE
    )


    ## count_na
    count_weight_save <- Reduce(`+`, lapply(sum_cols_list, function(x) {
        x[!is.na(x)] <- 1
        x[is.na(x)] <- 0
        return(x)
    }))

    ## sum all for each cell
    sum_cols_save <- Reduce(`+`, lapply(sum_cols_list, function(x) {
        ## change NA cell to 0 cell
        x[is.na(x)] <- 0

        return(x)
    }))

    ## get the mean for each cell
    ## only when subtype is supported by consensus_subtypes
    ## for each iteration and each patient
    mean_cols_save <- sum_cols_save / count_weight_save

    ## get the mean testdata
    sum_cols_list.testdata <- mapply(
        function(res_ihc) {
            res_ihc$testData
        },
        res_ihc_iterative,
        SIMPLIFY = FALSE,
        USE.NAMES = FALSE
    )

    ## sum all for each cell
    mean_cols_save.testdata <- Reduce(
        `+`,
        sum_cols_list.testdata
    ) / length(sum_cols_list.testdata)


    ## distances.Subtype for ROR
    sum_cols_list.Subtype <- mapply(
        function(res_ihc) {
            # res_ihc = res_ihc_iterative[[1]]

            res_ihc$distances.Subtype <- data.frame(res_ihc$dist.RORSubtype)

            ## if FALSE, make the cell as NULL
            keep <- res_ihc$predictions == consensus_subtypes
            res_ihc$distances.Subtype[!keep, ] <-
                as.list(rep(NA, ncol(res_ihc$distances.Subtype)))

            res <- dplyr::mutate(
                res_ihc$distances.Subtype,
                dplyr::across(
                    dplyr::everything(),
                    ~ ifelse(!is.na(.x), as.numeric(as.character(.x)), NA_real_)
                )
            )

            return(res)
        },
        res_ihc_iterative,
        SIMPLIFY = FALSE,
        USE.NAMES = FALSE
    )

    ## count_na
    count_weight_save.Subtype <- Reduce(
        `+`,
        lapply(sum_cols_list.Subtype, function(x) {
            x[!is.na(x)] <- 1
            x[is.na(x)] <- 0
            return(x)
        })
    )


    ## sum all for each cell
    sum_cols_save.Subtype <- Reduce(
        `+`,
        lapply(sum_cols_list.Subtype, function(x) {
            ## change NA cell to 0 cell
            x[is.na(x)] <- 0
            return(x)
        })
    )

    ## get the mean for each cell
    ## only when subtype is supported by consensus_subtypes
    ## for each iteration and each patient
    sum_cols_save.Subtype <- sum_cols_save.Subtype / count_weight_save.Subtype

    res <- list(
        mean_distance = mean_cols_save,
        mean_distance.Subtype = sum_cols_save.Subtype,
        testdata = mean_cols_save.testdata
    )


    return(res)
}

#### Visualization  ####

#' Boxplot of Correlation per Subtype
#'
#' @name Vis_boxplot
#' @description This function generates a boxplot to visualize the correlation
#' distribution between different subtypes of breast cancer, based on the
#' provided correlation table and subtype information.
#'
#' @param out A data frame containing the columns `"PatientID"` and `"Subtype"`.
#'   The `"PatientID"` column should have unique identifiers for each patient,
#'   and the `"Subtype"` column should specify the assigned subtype for each
#'   patient.
#' @param correlations A data frame or matrix containing the correlation values
#'   computed from NC-based methods.
#'
#' @return A `ggplot` object representing the boxplot visualization of the
#'   correlation distributions across the different subtypes.
#'
#' @details
#' Each plotted value is the maximum across one row of `correlations`.
#' Rows are paired positionally with `out$Subtype`; keep the two inputs in
#' the same sample order. Neutral boxes retain the usual quartiles, median,
#' whiskers and outliers. The maximum correlation is not a probability.
#'
#' @examples
#' data("OSLO2EMIT0obj")
#' res <- OSLO2EMIT0obj$res
#'
#' # Prepare data: Subtype information and correlation matrix
#' out <- data.frame(
#'     PatientID = res$results$genefu.robust$BS.all$PatientID,
#'     Subtype = res$results$genefu.robust$BS.all$BS
#' )
#'
#' correlations <- res$results$genefu.robust$outList$distances
#'
#' # Generate the boxplot
#' p <- Vis_boxplot(out, correlations)
#' plot(p)
#'
#' @export

Vis_boxplot <- function(out, correlations) {
    df <- data.frame(
        predictions = out$Subtype,
        cor = apply(correlations, 1, max)
    )

    plot <- ggplot(df, aes(x = .data$predictions, y = .data$cor)) +
        geom_boxplot() +
        labs(x = "", y = "Maximum centroid correlation") +
        .bsr_plot_theme()

    return(plot)
}


#' Heatmap Visualization of Gene Expression by Subtype
#'
#' @name Vis_heatmap
#' @description This function generates a heatmap to visualize gene expression
#' patterns across breast cancer subtypes, based on the provided gene expression
#' matrix and subtype information.
#'
#' @param x A gene expression matrix, where genes are rows and samples are
#'   columns. The data should be log2 transformed.
#' @param out A data frame containing two columns: `"PatientID"` and
#'   `"Subtype"`. The `"PatientID"` column should contain unique patient
#'   identifiers, and the `"Subtype"` column should specify the assigned subtype
#'   for each patient.
#'
#' @return A `ComplexHeatmap::Heatmap` object representing gene expression
#'   across the subtype-specific sample slices.
#'
#' @details
#' Genes are scaled with `t(scale(t(x[, out$PatientID])))`. The blue, white
#' and red expression scale is anchored at the observed minimum, zero and
#' maximum gene z-score, respectively; the scale is not clipped or made
#' symmetric. Subtype annotations retain the standard subtype colors.
#' Missing subtype calls are grey and keep their expression measurements.
#'
#' Gene and within-slice sample clustering use Spearman distance and average
#' linkage. The existing dendrogram reordering and ComplexHeatmap clustering
#' of column-slice means are retained. Both gene and sample dendrograms are
#' displayed. Horizontal legends occupy a reserved footer, without changing
#' the returned class. For 50 genes, a device around 10 by 7.4 inches is
#' recommended; use RStudio's Zoom view for a small plot pane, PDF for
#' vector output or PNG at 300 dpi.
#'
#' @examples
#' library(SummarizedExperiment)
#' data("OSLO2EMIT0obj")
#' res <- OSLO2EMIT0obj$res
#'
#' # Prepare data: Gene expression matrix and subtype information
#' x <- assay(OSLO2EMIT0obj$data_input$se_NC)
#' out <- data.frame(
#'     PatientID = res$results$genefu.robust$BS.all$PatientID,
#'     Subtype = res$results$genefu.robust$BS.all$BS
#' )
#'
#' # Generate the heatmap
#' p <- Vis_heatmap(x, out)
#' plot(p)
#'
#' @export

Vis_heatmap <- function(x, out) {
    scaled_mat <- t(scale(t(x[, out$PatientID])))

    ## color
    scaled_dat <- c(min(scaled_mat), 0, max(scaled_mat))
    col_fun <- circlize::colorRamp2(
        scaled_dat, c("#2166AC", "#FFFFFF", "#B2182B"))

    ## column annotation
    col_anno <- data.frame(
        row.names = out$PatientID,
        Subtype = out$Subtype
    )
    anno_col <- ComplexHeatmap::HeatmapAnnotation(
        df = col_anno,
        show_legend = FALSE,
        annotation_name_gp = grid::gpar(fontsize = 10),
        col = list(
            Subtype = c(
                "Basal" = "red",
                "Her2" = "hotpink",
                "LumA" = "darkblue",
                "LumB" = "skyblue",
                "Normal" = "green"
            )
        )
    )

    subtype_colors <- c(Basal = "red", Her2 = "hotpink", LumA = "darkblue",
                        LumB = "skyblue", Normal = "green")
    present <- names(subtype_colors) %in% out$Subtype
    legends <- if (any(present)) list(.bsr_legend_style(
        "Intrinsic subtype", labels = names(subtype_colors)[present], nrow = 1,
        legend_gp = grid::gpar(fill = subtype_colors[present]))) else list()
    legends[[length(legends) + 1L]] <- .bsr_legend_style(
        "Expression: gene z-score", col_fun = col_fun, at = scaled_dat,
        labels = sprintf("%.2f", scaled_dat), direction = "horizontal",
        legend_width = grid::unit(55, "mm"))
    if (anyNA(out$Subtype)) legends[[length(legends) + 1L]] <- .bsr_missing_key()
    footer <- .bsr_legend_footer("expression_keys", legends)

    heatmap <- ComplexHeatmap::Heatmap(
        scaled_mat,
        name = "Expr",
        col = col_fun,
        ## annotation
        top_annotation = anno_col,
        bottom_annotation = footer$annotation,
        post_fun = footer$draw,

        ## clustering
        ## as original heatmap plot
        clustering_distance_rows = "spearman",
        clustering_method_rows = "average",
        show_row_dend = TRUE,
        row_dend_width = grid::unit(15, "mm"),
        row_dend_gp = grid::gpar(col = "#4B5563", lwd = 0.65),
        column_dend_gp = grid::gpar(col = "#4B5563", lwd = 0.7),
        cluster_column_slices = TRUE,
        column_split = col_anno$Subtype,
        clustering_distance_columns = "spearman",
        clustering_method_columns = "average",

        ## general
        show_column_names = FALSE,
        show_heatmap_legend = FALSE,
        column_title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
        row_names_gp = grid::gpar(fontsize = 8, fontfamily = "sans")
    )

    return(heatmap)
}




#' PCA Plot Visualization of Gene Expression by Subtype
#'
#' @name Vis_PCA
#' @description This function generates a PCA plot to visualize the principal
#' components of gene expression data, colored by the assigned subtypes.
#' Optionally, it can display a scree plot of eigenvalues to evaluate the
#' explained variance.
#'
#' @param x A gene expression matrix, where genes are rows and samples are
#'   columns. The data should be log2 transformed.
#' @param out A data frame containing two columns: `"PatientID"` and
#'   `"Subtype"`. The `"PatientID"` column should contain unique patient
#'   identifiers, and the `"Subtype"` column should specify the assigned subtype
#'   for each patient.
#' @param Eigen Logical. If `TRUE`, the function will display a scree plot
#'   showing the eigenvalues of the principal components.
#'
#' @return A `ggplot` object representing the PCA plot, colored by subtype.
#'   If `Eigen = TRUE`, returns only the scree plot.
#'
#' @details
#' The scatter plot uses a shared sans-serif theme and a bottom subtype
#' legend. PCA centering, scaling, scores, joins, explained variance and
#' reference lines are unchanged; no confidence ellipses or tests are added.
#'
#' @examples
#' library(SummarizedExperiment)
#' data("OSLO2EMIT0obj")
#' res <- OSLO2EMIT0obj$res
#'
#' # Prepare data: Gene expression matrix and subtype information
#' x <- assay(OSLO2EMIT0obj$data_input$se_NC)
#' out <- data.frame(
#'     PatientID = res$results$genefu.robust$BS.all$PatientID,
#'     Subtype = res$results$genefu.robust$BS.all$BS
#' )
#'
#' # Generate the PCA plot
#' p <- Vis_PCA(x = x, out = out)
#' plot(p)
#'
#' # Generate PCA plot with scree plot of eigenvalues
#' p_with_eigen <- Vis_PCA(x = x, out = out, Eigen = TRUE)
#' plot(p_with_eigen)
#'
#' @export

Vis_PCA <- function(x, out, Eigen = FALSE) {
    Subtype.color <- c(
        "Basal" = "red",
        "Her2" = "hotpink",
        "LumA" = "darkblue",
        "LumB" = "skyblue",
        "Normal" = "green"
    )
    pca <- prcomp(t(x), center = TRUE, scale. = TRUE)

    if (Eigen) {
        # Scree plot
        variance <- pca$sdev^2 / sum(pca$sdev^2) * 100
        scree_data <- data.frame(PC = seq_along(variance), Variance = variance)

        screeplot <- ggplot(
            scree_data[seq(1, 10), ],
            aes(x = .data$PC, y = .data$Variance)
        ) +
            geom_bar(
                stat = "identity",
                width = 0.5,
                fill = "steelblue"
            ) +
            geom_line() +
            geom_point(size = 2) +
            scale_x_continuous(breaks = seq(1, 10, 1)) +
            geom_text(
                aes(
                    x = .data$PC,
                    label = paste0(round(.data$Variance, 2), "%")
                ),
                nudge_y = 1
            ) +
            labs(
                x = "Principal Component",
                y = "Percentage of variance Explained"
            ) +
            theme_classic() +
            theme(
                axis.text = element_text(size = 12),
                axis.title = element_text(size = 14)
            )

        return(screeplot)
    } else {
        ## PCA plot
        scores <- as.data.frame(pca$x)
        scores$PatientID <- rownames(scores)
        scores <- dplyr::left_join(scores, out, by = "PatientID")
        rownames(scores) <- scores$PatientID

        pc1_var <- round(100 * summary(pca)$importance[2, 1], 2)
        pc2_var <- round(100 * summary(pca)$importance[2, 2], 2)

        pcaplot <- ggplot(
            data = scores,
            aes(x = .data$PC1, y = .data$PC2, color = .data$Subtype)
        ) +
            geom_point(size = 2) +
            geom_vline(xintercept = 0, linetype = "dashed") +
            geom_hline(yintercept = 0, linetype = "dashed") +
            scale_color_manual(name = "Subtype", values = Subtype.color) +
            labs(
                x = paste0("PC1 (", pc1_var, "% variance)"),
                y = paste0("PC2 (", pc2_var, "% variance)")
            ) +
            .bsr_plot_theme() +
            theme(
                legend.position = "bottom",
                panel.grid.major = ggplot2::element_line(
                    colour = "#ECEFF2", linewidth = 0.25)
            )

        return(pcaplot)
    }
}

#' Pie Chart Visualization of Subtype Distribution
#'
#' @name Vis_pie
#' @description This function generates a pie chart to visualize the
#' distribution of breast cancer subtypes in a cohort, based on the provided
#' `Subtype` data.
#'
#' @param out A data frame containing two columns: `"PatientID"` and
#'   `"Subtype"`. The `"PatientID"` column should contain unique patient
#'   identifiers, and the `"Subtype"` column should specify the assigned subtype
#'   for each patient.
#'
#' @return A `ggplot` object representing a pie chart showing the proportion of
#'   each subtype in the dataset.
#'
#' @details
#' The legend reports subtype counts and percentages. As with `table()`,
#' missing subtype calls are excluded from the assigned-call denominator;
#' the caption reports both counts. Percentages are rounded to two decimals
#' before drawing the slices, preserving the existing slice geometry and
#' subtype order. Missing calls are not added as an extra slice.
#'
#' @examples
#' data("OSLO2EMIT0obj")
#' res <- OSLO2EMIT0obj$res
#'
#' # Prepare data: Subtype information
#' out <- data.frame(
#'     PatientID = res$results$genefu.robust$BS.all$PatientID,
#'     Subtype = res$results$genefu.robust$BS.all$BS
#' )
#'
#' # Generate the pie chart
#' p <- Vis_pie(out = out)
#' plot(p)
#'
#' @export

Vis_pie <- function(out) {
    data <- data.frame(table(out$Subtype))
    data <- data %>%
        dplyr::mutate(
            perc = round(.data$`Freq` / sum(.data$`Freq`) * 100, 2),
            csum = rev(cumsum(rev(.data$perc))),
            pos = .data$perc / 2 + dplyr::lead(.data$csum, 1),
            pos = dplyr::if_else(is.na(.data$pos), .data$perc / 2, .data$pos)
        )


    Subtype.color <- c(
        "Basal" = "red",
        "Her2" = "hotpink",
        "LumA" = "darkblue",
        "LumB" = "skyblue",
        "Normal" = "green"
    )

    plot <- ggplot(data, aes(x = "", y = .data$perc, fill = .data$Var1)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(
            name = "Subtype   n (%)", values = Subtype.color,
            breaks = as.character(data$Var1),
            labels = sprintf("%s   %d (%.2f%%)", data$Var1, data$Freq, data$perc)
        ) +
        coord_polar("y", direction = -1) +
        labs(
            title = "Subtype distribution",
            subtitle = sprintf("%d samples", nrow(out)),
            caption = sprintf(paste0(
                "Denominator: %d assigned calls; %d missing calls excluded.\n",
                "Percentages retain the baseline two-decimal rounding."),
                sum(data$Freq), sum(is.na(out$Subtype)))
        ) +
        theme_void(base_size = 11, base_family = "sans") +
        .bsr_plot_theme() +
        theme(
            axis.line = ggplot2::element_blank(),
            axis.text = ggplot2::element_blank(),
            axis.title = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),
            legend.position = "right"
        )


    return(plot)
}


#' Multi-Method Subtype Heatmap Visualization
#'
#' @name Vis_Multi
#' @description This function generates a heatmap to visualize breast cancer subtypes
#' classified by multiple subtyping methods. It helps users compare how different methods
#' assign subtypes to the same set of samples.
#'
#' @param data Output of the \code{\link{BS_Multi}} function.
#'
#' @return Returns a heatmap visualizing the subtype classifications across multiple methods.
#'
#' @details
#' Returns a `ComplexHeatmap::Heatmap` with separate subtype, method-family
#' and (when needed) missing-call keys in a horizontal footer. The method
#' rows, sample sorting rules, subtype colors and grey missing cells are
#' preserved. The entropy annotation displays the supplied raw Shannon
#' entropy in bits: disagreement among available method calls, not calibrated
#' confidence. No entropy values are recomputed or normalized by this plot.
#' A device around 10 by 4.6 inches leaves room for labels and legends.
#'
#' @examples
#' data("OSLO2EMIT0obj")
#'
#' # Assuming `OSLO2EMIT0obj$res$res_subtypes` contains multi-method subtype results
#' p <- Vis_Multi(OSLO2EMIT0obj$res$res_subtypes)
#' plot(p)
#'
#' @export

Vis_Multi <- function(data) {
    ## order by subtype
    data <- data[order(data[, ncol(data) - 1]), ]
    ## order by entropy
    data <- data[order(data[, ncol(data)], decreasing = FALSE), ]
    mat <- data[, -ncol(data)]

    Labels <- unique(as.vector(as.matrix(mat)))

    ## preset
    ctg <- data.frame(
        Cat = rep(c("NC-based", "SSP-based"), c(8, 2)),
        row.names = c(
            "parker.original",
            "genefu.scale",
            "genefu.robust",
            "cIHC",
            "cIHC.itr",
            "PCAPAM50",
            "ssBC",
            "ssBC.v2",
            "AIMS",
            "sspbc"
        )
    )
    ctg$Cat <- factor(ctg$Cat, levels = c("NC-based", "SSP-based"))

    ## color
    Subtype.color <- c(
        "Basal" = "red",
        "Her2" = "hotpink",
        "LumA" = "darkblue",
        "LumB" = "skyblue",
        "Normal" = "green"
    )
    col_cat <- setNames(c("#984ea3", "#ff7f00"), c("NC-based", "SSP-based"))

    ## make row annotation
    row_anno <- data.frame(
        Category = ctg[colnames(mat), ],
        row.names = colnames(mat)
    )
    row_anno <- HeatmapAnnotation(
        df = row_anno, which = "row", col = list(Category = col_cat),
        show_legend = FALSE, show_annotation_name = FALSE,
        simple_anno_size = grid::unit(4, "mm")
    )
    col_anno <- HeatmapAnnotation(
        which = "column",
        Entropy = anno_barplot(
            data[rownames(data), ncol(data)], bar_width = 1,
            gp = grid::gpar(fill = "#4B5563", col = NA),
            height = grid::unit(20, "mm"),
            axis_param = list(gp = grid::gpar(fontsize = 9))),
        annotation_label = "Entropy\n(bits)",
        annotation_name_gp = grid::gpar(fontsize = 10)
    )

    Subtype.color <- Subtype.color[names(Subtype.color) %in% Labels]
    legends <- if (length(Subtype.color)) list(.bsr_legend_style(
        "Intrinsic subtype", labels = names(Subtype.color), nrow = 1,
        legend_gp = grid::gpar(fill = Subtype.color))) else list()
    category_colors <- col_cat[names(col_cat) %in% ctg[colnames(mat), "Cat"]]
    if (length(category_colors)) legends[[length(legends) + 1L]] <-
        .bsr_legend_style("Method family", labels = names(category_colors),
                          nrow = 1, legend_gp = grid::gpar(fill = category_colors))
    if (anyNA(mat)) legends[[length(legends) + 1L]] <- .bsr_missing_key()
    footer <- .bsr_legend_footer("method_keys", legends)
    p <- ComplexHeatmap::Heatmap(
        t(as.matrix(mat)),
        name = "Subtypes", col = Subtype.color,
        row_names_gp = grid::gpar(fontsize = 10, fontfamily = "sans"),
        right_annotation = row_anno, top_annotation = col_anno,
        bottom_annotation = footer$annotation, post_fun = footer$draw,
        show_column_names = FALSE, show_heatmap_legend = FALSE
    )

    return(p)
}
