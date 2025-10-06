#'
#' Functions adapted from NC-based subtyping methods
#' @import magrittr
#' @import impute
#' @importFrom withr with_seed
#' @importFrom dplyr select
#' @importFrom tidyselect everything
#' @noRd
NULL

#'
#' Function for central median
#' @param x gene expression matrix
#' @noRd
medianCtr <- function(x) {
    annAll <- dimnames(x)
    medians <- apply(x, 1, median, na.rm = TRUE)
    x <- t(scale(t(x), center = medians, scale = FALSE))
    dimnames(x) <- annAll
    return(x)
}

#' Function to rescale values based on quantiles
#' This is adapted from genefu package
#' @param x Gene expression matrix or vector
#' @noRd
rescale <- function(x, na.rm = FALSE, q = 0) {
    if (q == 0) {
        ma <- max(x, na.rm = na.rm)
        mi <- min(x, na.rm = na.rm)
    } else {
        ma <- quantile(x, probs = 1 - (q / 2), na.rm = na.rm, names = FALSE, type = 7)
        mi <- quantile(x, probs = q / 2, na.rm = na.rm, names = FALSE, type = 7)
    }
    rng <- ma - mi
    if (!is.finite(rng) || rng == 0) {
        return(x * 0)
    } # avoid NaN/Inf
    xx <- (x - mi) / rng

    return(xx)
}

#' Function for ordering genes in expression matrix as PAM50 genes
#' @param x centroid matrix
#' @param y Gene expression matrix
#' @noRd
overlapSets <- function(x, y) {
    # subset the two lists to have a commonly ordered gene list
    x <- x[dimnames(x)[[1]] %in% dimnames(y)[[1]], ]
    y <- y[dimnames(y)[[1]] %in% dimnames(x)[[1]], ]

    # and sort such that thing are in the correct order
    x <- x[sort.list(row.names(x)), ]
    y <- y[sort.list(row.names(y)), ]

    return(list(x = x, y = y))
}

# Align multiple expression matrices by common gene rows, then cbind (columns)
.align_and_cbind_by_genes <- function(...) {
    mats <- list(...)
    # keep only non-empty matrices
    mats <- Filter(function(m) !is.null(m) && is.matrix(m) && nrow(m) > 0 && ncol(m) > 0, mats)
    if (!length(mats)) {
        # return an empty, well-formed matrix
        return(matrix(numeric(0),
            nrow = 0, ncol = 0,
            dimnames = list(NULL, NULL)
        ))
    }

    # require rownames on all
    if (any(vapply(mats, function(m) is.null(rownames(m)), logical(1)))) {
        stop("All subgroup testData matrices must have rownames (gene IDs).")
    }

    common <- Reduce(intersect, lapply(mats, rownames))
    if (!length(common)) {
        # no common genes: return a 0-row matrix with all sample columns
        all_cols <- unlist(lapply(mats, colnames), use.names = FALSE)
        return(matrix(numeric(0),
            nrow = 0, ncol = length(all_cols),
            dimnames = list(NULL, all_cols)
        ))
    }

    mats_aligned <- lapply(mats, function(m) m[common, , drop = FALSE])
    # keep a consistent gene order across all mats
    mats_aligned <- lapply(mats_aligned, function(m) m[match(common, rownames(m)), , drop = FALSE])
    do.call(cbind, mats_aligned)
}

#' Function for calibration methods
#' @param y gene x sample matrix
#' @param medians.all matrix/data.frame of reference medians; rownames=genes, cols=platforms (e.g., nCounter, RNAseq.V2, Given.mdns, IHC.mdns, etc.)
#' @param calibration "None", "Internal", or "External"
#' @param internal    For Internal: one of "medianCtr" (default), "meanCtr", "qCtr", or the NAME of a column in medians.all (e.g., "IHC.mdns")
#' @param external    For External: NAME of a column in medians.all (e.g., "RNAseq.V2", "Given.mdns")
#' @noRd
docalibration <- function(
        y,
        medians.all,
        calibration = c("None", "Internal", "External"),
        internal = NA,
        external = NA) {
    calibration <- match.arg(calibration)
    y <- as.matrix(y)
    mq <- 0.05 # robust quantiles for 'qCtr' (genefu.robust)

    # ---- No calibration ----
    if (calibration == "None") {
        message("No calibration")
        return(as.matrix(y))
    }

    # ---- Internal ----
    if (calibration == "Internal") {
        if (length(internal) == 0L || is.na(internal) || identical(internal, "medianCtr")) {
            y <- medianCtr(y)
        } else if (identical(internal, "meanCtr")) {
            y <- t(scale(t(y), center = TRUE, scale = TRUE))
        } else if (identical(internal, "qCtr")) {
            y <- t(apply(y, 1, function(x) (rescale(x, q = mq, na.rm = TRUE) - 0.5) * 2))
        } else if (!is.null(medians.all) && !is.null(colnames(medians.all)) && internal %in% colnames(medians.all)) {
            # use a cohort-derived *internal* medians column (e.g., "IHC.mdns")
            tm <- overlapSets(medians.all, y) # align genes
            med <- tm$x[, internal, drop = TRUE]
            y <- tm$y - med
        } else {
            stop(
                "Invalid internal calibration '", internal,
                "'. Use 'medianCtr', 'meanCtr', 'qCtr', or a column present in 'medians.all'."
            )
        }
        return(y)
    }


    # ---- External ----
    if (length(external) == 0L || is.na(external)) {
        stop("For calibration='External', please supply 'external' as a column name in medians.")
    }
    if (is.null(medians.all) || is.null(colnames(medians.all)) || !(external %in% colnames(medians.all))) {
        stop(
            "Column '", external, "' not found in medians table. Available: ",
            paste(colnames(medians.all), collapse = ", ")
        )
    }
    tm <- overlapSets(medians.all, y) # align genes
    y <- tm$y - tm$x[, external, drop = TRUE] # subtract EXTERNAL medians

    return(y)
}



#' Function for standardization
#' @param x gene expression matrix
#' @noRd
standardize <- function(x) {
    annAll <- dimnames(x)
    x <- scale(x)
    dimnames(x) <- annAll
    return(x)
}


#' Function for suffix of medians for gene centering
#' @noRd
getsuffix <- function(calibration,
    internal = NA,
    external = NA) {
    calibration <- match.arg(calibration, choices = c("None", "Internal", "External"))

    suffix <- switch(calibration,
        "None" = "None",
        "Internal" = internal,
        "External" = external
    )
    return(suffix)
}


#' Function for predicting PAM50 subtypes
#' @param x median train file
#' @param y gene expression matrix
#' @param std Logic
#' @param distm "spearman" (default), "euclidean", "correlation" or "pearson"
#' @param Subtype Logic. Please specify if it predicts Subtype-like subtype
#' @noRd
sspPredict <- function(x, y, std = FALSE, distm = "spearman", Subtype = TRUE) {
    distm <- match.arg(distm, choices = c("spearman", "euclidean", "correlation", "pearson"))

    dataMatrix <- x
    tdataMatrix <- y

    tmp <- overlapSets(dataMatrix, tdataMatrix)
    dataMatrix <- tmp$x
    tdataMatrix <- tmp$y

    dataMatrix_main <- dataMatrix
    tdataMatrix_main <- tdataMatrix

    if (nrow(dataMatrix) == 0L || ncol(tdataMatrix) == 0L) {
        # return an empty but well-formed object
        pred <- setNames(character(0), character(0))
        res <- list(
            predictions     = pred,
            testData        = as.matrix(tdataMatrix_main),
            distances       = matrix(numeric(0), nrow = 0, ncol = ncol(dataMatrix), dimnames = list(NULL, colnames(dataMatrix))),
            dist.RORSubtype = matrix(numeric(0), nrow = 0, ncol = min(4, ncol(dataMatrix)), dimnames = list(NULL, colnames(dataMatrix)[seq_len(min(4, ncol(dataMatrix)))])),
            centroids       = dataMatrix_main
        )
        if (Subtype) res$predictions.FourSubtype <- pred
        return(res)
    }

    sfeatureNames <- rownames(dataMatrix)

    if (std) {
        dataMatrix <- standardize(dataMatrix)
        tdataMatrix <- standardize(tdataMatrix)
    }

    nClasses <- ncol(dataMatrix)
    classLevels <- colnames(dataMatrix)

    distances <- matrix(NA_real_,
        ncol = nClasses, nrow = ncol(tdataMatrix),
        dimnames = list(colnames(tdataMatrix), colnames(dataMatrix))
    )

    for (j in seq_len(nClasses)) {
        if (distm == "euclidean") {
            combined <- cbind(dataMatrix[, j], tdataMatrix)
            dv <- dist(t(combined))
            distances[, j] <- dv[seq_len(ncol(tdataMatrix))]
        } else {
            mth <- if (distm == "spearman") "spearman" else "pearson"
            distances[, j] <- apply(tdataMatrix, 2, function(xx) {
                -cor(dataMatrix[, j], xx, method = mth, use = "pairwise.complete.obs")
            })
        }
    }

    # If a row is all NA (no finite distances), prediction := NA
    prediction <- rep(NA_character_, nrow(distances))
    names(prediction) <- rownames(distances)
    good <- rowSums(is.finite(distances)) > 0
    if (any(good)) {
        prediction[good] <- classLevels[apply(distances[good, , drop = FALSE], 1, which.min)]
    }

    if (Subtype) {
        nClasses4 <- min(4, nClasses)
        classLevels4 <- classLevels[seq_len(nClasses4)]
        dist4 <- matrix(NA_real_,
            ncol = nClasses4, nrow = ncol(tdataMatrix),
            dimnames = list(colnames(tdataMatrix), classLevels4)
        )
        for (j in seq_len(nClasses4)) {
            if (distm == "euclidean") {
                combined <- cbind(dataMatrix[, j], tdataMatrix)
                dv <- dist(t(combined))
                dist4[, j] <- dv[seq_len(ncol(tdataMatrix))]
            } else {
                mth <- if (distm == "spearman") "spearman" else "pearson"
                dist4[, j] <- apply(tdataMatrix, 2, function(xx) {
                    -cor(dataMatrix[, j], xx, method = mth, use = "pairwise.complete.obs")
                })
            }
        }
        prediction.Subtype <- rep(NA_character_, nrow(dist4))
        names(prediction.Subtype) <- rownames(dist4)
        good4 <- rowSums(is.finite(dist4)) > 0
        if (any(good4)) {
            prediction.Subtype[good4] <- classLevels4[apply(dist4[good4, , drop = FALSE], 1, which.min)]
        }
    }

    ## ROR distances: drop the 4 excluded genes first
    genes.ex <- c("BIRC5", "CCNB1", "GRB7", "MYBL2")
    t4 <- tdataMatrix_main[!(rownames(tdataMatrix_main) %in% genes.ex), , drop = FALSE]
    tmp4 <- overlapSets(dataMatrix, t4)
    cent4 <- tmp4$x[, seq_len(min(4, ncol(tmp4$x))), drop = FALSE]
    dat4 <- tmp4$y

    dist.RORSubtype <- matrix(NA_real_,
        ncol = ncol(cent4), nrow = ncol(dat4),
        dimnames = list(colnames(dat4), colnames(cent4))
    )
    if (ncol(dat4) > 0L) {
        for (j in seq_len(ncol(cent4))) {
            dist.RORSubtype[, j] <- apply(dat4, 2, function(xx) {
                -cor(cent4[, j], xx, method = "spearman", use = "pairwise.complete.obs")
            })
        }
    }

    if (Subtype) {
        res <- list(
            predictions = prediction,
            predictions.FourSubtype = prediction.Subtype,
            testData = as.matrix(tdataMatrix_main),
            distances = distances,
            dist.RORSubtype = dist.RORSubtype,
            centroids = dataMatrix_main
        )
    } else {
        res <- list(
            predictions            = prediction,
            testData               = as.matrix(tdataMatrix_main),
            distances              = distances,
            dist.RORSubtype        = dist.RORSubtype,
            centroids              = dataMatrix_main
        )
    }
    res
}


#' Function for risk calculation
#'
#' @param out The result of sspPredict() function.
#' @param out clinical table
#' @param Subtype Logic.
#' @param hasClinical Logic. Specify whether clinical information is included.
#'   For example, tumor size should be in the "T" column, and lymph node status
#' - "NODE": Lymph node status (0 = negative; \eqn{\ge 1}{>= 1} = positive).
#' @return ROR, ROR risk group and other indications
#' @noRd

RORgroup <- function(out,
    df.cln,
    Subtype = FALSE,
    hasClinical = FALSE) {
    sample <- data.frame(patientID = names(out$predictions))

    distance <- data.frame(out$distances, row.names = names(out$predictions))
    colnames(distance) <- c("Basal", "Her2", "LumA", "LumB", "Normal")

    Call <- data.frame(
        "Call" = out$predictions,
        row.names = names(out$predictions)
    )

    # providing proliferation signatures
    prolifG <- c(
        "CCNB1",
        "UBE2C",
        "BIRC5",
        "KNTC2",
        "CDC20",
        "PTTG1",
        "RRM2",
        "MKI67",
        "TYMS",
        "CEP55",
        "CDCA1"
    )
    prolifG.Subtype <- c(
        "ANLN",
        "CEP55",
        "ORC6L",
        "CCNE1",
        "EXO1",
        "PTTG1",
        "CDC20",
        "KIF2C",
        "RRM2",
        "CDC6",
        "KNTC2",
        "TYMS",
        "CDCA1",
        "MELK",
        "UBE2C",
        "CENPF",
        "MKI67",
        "UBE2T"
    )

    ###
    # some constants for ROR groups
    ###

    # for subtype only model
    glthreshold <- -0.15
    ghthreshold <- 0.1

    # for subtype + proliferation model
    gplthreshold <- -0.25
    gphthreshold <- 0.1

    # for combined model
    clthreshold <- -0.1
    chthreshold <- 0.2

    # for combined + proliferation model
    cplthreshold <- -0.2
    cphthreshold <- 0.2

    # thresholds for combined + proliferation + Subtype
    # node-negative
    cpl.Subtype.NODE.0 <- 40
    cph.Subtype.NODE.0 <- 60
    # node-positive (1–3 nodes)
    cpl.Subtype.NODE <- 15
    cph.Subtype.NODE <- 40
    ## >= 4 nodes would be assigned as high-risk

    ## ER and HER2 score
    # out$testData = mat

    # Safe extraction of ESR1 / ERBB2 scores; fill with NA if missing
    n_samp <- if (is.matrix(out$testData)) ncol(out$testData) else 0L
    samp_names <- if (n_samp) colnames(out$testData) else character(0)

    erScore <- if (!is.null(rownames(out$testData)) && "ESR1" %in% rownames(out$testData)) {
        out$testData["ESR1", ]
    } else {
        setNames(rep(NA_real_, n_samp), samp_names)
    }

    her2Score <- if (!is.null(rownames(out$testData)) && "ERBB2" %in% rownames(out$testData)) {
        out$testData["ERBB2", ]
    } else {
        setNames(rep(NA_real_, n_samp), samp_names)
    }


    er_her2 <- data.frame(
        "ER" = erScore,
        "HER2" = her2Score,
        row.names = colnames(out$testData)
    )


    # calculate the proliferation score
    # Filter test data for proliferation genes
    prolifG_data <- out$testData[rownames(out$testData) %in% prolifG, ]
    prolifScore <- apply(prolifG_data, 2, mean, na.rm = TRUE)

    # Filter test data for proliferation Subtype genes
    prolifG_Subtype <- out$testData[rownames(out$testData) %in% prolifG.Subtype, ]
    prolifScore_Subtype <- apply(prolifG_Subtype, 2, mean, na.rm = TRUE)

    ## confidence
    call.conf <- c()
    # for (j in seq_len(length(out$predictions))) {
    #     centroid_col <- which(colnames(out$centroids) == out$predictions[j])
    #     call.conf[j] <- 1 - cor.test(out$testData[, j],
    #         out$centroids[, centroid_col],
    #         method = "spearman",
    #         exact = FALSE
    #     )$p.value
    # }

    for (j in seq_len(length(out$predictions))) {
        centroid_col <- which(colnames(out$centroids) == out$predictions[j])
        x <- out$testData[, j, drop = TRUE]
        y <- out$centroids[, centroid_col, drop = TRUE]
        common <- intersect(names(x), rownames(out$centroids))
        # if names are NULL, fall back to positional intersection with rownames
        if (length(common) > 0L) {
            xi <- x[common]
            yi <- y[common]
        } else {
            # assume already aligned; trim to min length to avoid accidental mismatch
            n <- min(length(x), length(y))
            xi <- x[seq_len(n)]
            yi <- y[seq_len(n)]
        }
        ok <- is.finite(xi) & is.finite(yi)
        pv <- tryCatch(
            {
                if (sum(ok) >= 3L) cor.test(xi[ok], yi[ok], method = "spearman", exact = FALSE)$p.value else 1
            },
            error = function(e) 1
        )
        call.conf[j] <- 1 - pv
    }
    call.conf <- as.numeric(call.conf)
    call.conf[!is.finite(call.conf)] <- NA_real_
    call.conf <- round(call.conf, 2)
    call.conf <- data.frame(
        "Confidence" = call.conf,
        row.names = names(out$predictions)
    )


    ## calculate the risk scores
    ## genomic
    ## basal, her2, lumA, lumB in order
    genomic <- 0.04210193 * out$distances[, 1] +
        0.12466938 * out$distances[, 2] + -0.35235561 * out$distances[, 3] +
        0.14213283 * out$distances[, 4]

    genomicWprolif <- -0.0009299747 * out$distances[, 1] +
        0.0692289192 * out$distances[, 2] +
        -0.0951505484 * out$distances[, 3] +
        0.0493487685 * out$distances[, 4] +
        0.3385116381 * prolifScore


    # threshold the risk score
    griskgroups <- genomic
    griskgroups[genomic > ghthreshold] <- "high"
    griskgroups[genomic > glthreshold &
        genomic < ghthreshold] <- "med"
    griskgroups[genomic < glthreshold] <- "low"
    gpriskgroups <- genomicWprolif
    gpriskgroups[genomicWprolif > gphthreshold] <- "high"
    gpriskgroups[genomicWprolif > gplthreshold &
        genomicWprolif < gphthreshold] <- "med"
    gpriskgroups[genomicWprolif < gplthreshold] <- "low"

    genomic <- 100 * (genomic + 0.35) / 0.85
    genomicWprolif <- 100 * (genomicWprolif + 0.35) / 0.85

    ROR.genomic <- data.frame(
        "ROR-S (Subtype Only)" = genomic,
        "ROR-S Group (Subtype Only)" = griskgroups,
        "ROR-P (Subtype + Prolif)" = genomicWprolif,
        "ROR-P Group (Subtype + Prolif)" = gpriskgroups,
        check.names = FALSE
    )

    if (hasClinical) {
        Clinical <- df.cln[match(names(out$predictions), df.cln$PatientID), ]
        rownames(Clinical) <- Clinical$PatientID

        if ("TSIZE" %in% colnames(Clinical)) {
            xT <- as.numeric(as.vector(Clinical$TSIZE))

            combined <- 0.0442770 * out$distances[, 1] +
                0.1170297 * out$distances[, 2] +
                -0.2608388 * out$distances[, 3] +
                0.1055908 * out$distances[, 4] +
                0.1813751 * xT

            combinedWprolif <- -0.009383416 * out$distances[, 1] +
                0.073725503 * out$distances[, 2] +
                -0.090436516 * out$distances[, 3] +
                0.053013865 * out$distances[, 4] +
                0.131605960 * xT + 0.327259375 * prolifScore

            criskgroups <- combined
            criskgroups[combined > chthreshold] <- "high"
            criskgroups[combined > clthreshold &
                combined < chthreshold] <- "med"
            criskgroups[combined < clthreshold] <- "low"
            cpriskgroups <- combinedWprolif
            cpriskgroups[combinedWprolif > cphthreshold] <- "high"
            cpriskgroups[combinedWprolif > cplthreshold &
                combinedWprolif < cphthreshold] <- "med"
            cpriskgroups[combinedWprolif < cplthreshold] <- "low"

            combined <- 100 * (combined + 0.35) / 0.85
            combinedWprolif <- 100 * (combinedWprolif + 0.35) / 0.85

            ROR.combined <- data.frame(
                "ROR-C (Subtype + Clinic)" = combined,
                "ROR-C Group (Subtype + Clinic)" = cpriskgroups,
                "ROR-PC (Subtype + Clinic + Prolif)" = combinedWprolif,
                "ROR-PC Group (Subtype + Clinic + Prolif)" = cpriskgroups,
                check.names = FALSE
            )

            cmbWprolif.Subtype <- 54.7690 * (
                -0.0067 * out$dist.RORSubtype[, 1] + 0.4317 * out$dist.RORSubtype[, 2] -
                    0.3172 * out$dist.RORSubtype[, 3] + 0.4894 * out$dist.RORSubtype[, 4] +
                    0.1981 * prolifScore_Subtype + 0.1133 * xT + 0.8826
            )


            ## check NODE
            ## grouping by NODE
            if ("NODE" %in% colnames(Clinical)) {
                Clinical$NODE <- as.numeric(Clinical$NODE)

                cprskg.Subtype <- cmbWprolif.Subtype

                patients.NODE <- rownames(Clinical)[!is.na(Clinical$NODE)]

                if (length(patients.NODE) > 0) {
                    ## NODE >= 4
                    patients.NODE.4 <- rownames(Clinical)[Clinical$NODE >= 4 &
                        !is.na(Clinical$NODE)]

                    if (length(patients.NODE.4) > 0) {
                        cprskg.Subtype[patients.NODE.4] <- "high"
                    }

                    ## NODE < 4 & >0
                    pts.NODE.3 <- rownames(Clinical)[Clinical$NODE < 4 &
                        Clinical$NODE > 0 &
                        !is.na(Clinical$NODE)]

                    if (length(pts.NODE.3) > 0) {
                        temp <- cmbWprolif.Subtype[pts.NODE.3]
                        cprskg.Subtype[pts.NODE.3][temp > cph.Subtype.NODE] <- "high"
                        cprskg.Subtype[pts.NODE.3][temp > cpl.Subtype.NODE &
                            temp < cph.Subtype.NODE] <- "med"
                        cprskg.Subtype[pts.NODE.3][temp < cpl.Subtype.NODE] <- "low"
                    }

                    ## negative node
                    pts.NODE.0 <- rownames(Clinical)[Clinical$NODE == 0 &
                        !is.na(Clinical$NODE)]

                    if (length(pts.NODE.0) > 0) {
                        tmp <- cmbWprolif.Subtype[pts.NODE.0]
                        cprskg.Subtype[pts.NODE.0][tmp > cph.Subtype.NODE.0] <- "high"
                        which.med <- tmp > cpl.Subtype.NODE.0 & tmp < cph.Subtype.NODE.0
                        cprskg.Subtype[pts.NODE.0][which.med] <- "med"
                        cprskg.Subtype[pts.NODE.0][tmp < cpl.Subtype.NODE.0] <- "low"
                    }
                }

                ## check if NA NODE
                ## grouping by NA NODE
                patients.NA <- rownames(Clinical)[is.na(Clinical$NODE) |
                    is.na(Clinical$TSIZE)]

                if (length(patients.NA) > 0) {
                    cprskg.Subtype[patients.NA] <- NA
                }

                ROR.combined.Subtype <- data.frame(
                    "ROR-PC (Subtype + Clinic + Prolif.Subtype)" = cmbWprolif.Subtype,
                    "ROR-PC Group (Subtype + Clinic + Prolif.Subtype)" = cprskg.Subtype,
                    check.names = FALSE
                )
            } else {
                message("NODE info is missing.")
            }
        } else {
            message("TSIZE info is missing.")
        }


        outtable <- cbind(
            sample,
            distance,
            Call,
            call.conf,
            ROR.genomic,
            ROR.combined,
            ROR.combined.Subtype,
            er_her2
        )
    } else {
        outtable <- cbind(
            sample,
            distance,
            Call,
            call.conf,
            ROR.genomic,
            er_her2
        )
    }


    if (Subtype) {
        ## pass distances (just omitting normal)
        Call.Subtype <- data.frame(
            "Call.Subtype" = out$predictions.FourSubtype,
            row.names = names(out$predictions.FourSubtype)
        )
        outtable <- cbind(outtable, Call.Subtype)
        ord <- c(
            "patientID",
            colnames(distance),
            "Call",
            "Call.Subtype",
            setdiff(colnames(outtable), c("patientID", colnames(distance), "Call", "Call.Subtype"))
        )
        outtable <- outtable[, ord, drop = FALSE]
    }

    return(outtable)
}




#'
#' Functions for predicting PAM50 intrinsic subtypes and calculation of
#' proliferation
#'

#### functions to make calls using nearest-centroid (NC-based) strategies ####

#' Function for calling PAM50 subtypes by Parker et al.-based methods Here, we
#' integrated Parker et al.-based methods and genefu PAM50 implementation
#' @param mat gene expression matrix, log of normalized
#' @param df.cln clinical information table with PatientID and IHC column
#' @param calibration How to do calibration, "None"(default) means no
#'   calibration for gene expression matrix. When setting calibration =None, you
#'   dont need to set internal and external parameters.  "Internal" means
#'   calibration for gene expression matrix by itself. "External" means
#'   calibration by external cohort.
#' @param internal Specify the strategy for internal calibration,
#'   "medianCtr" (median centered, default), meanCtr and qCtr
#' @param external Specify the platform name(which column) of external medians
#'   calculated by train cohorts. When users want to use Medians prepared by
#'   user selves, this parameter should be "Given.mdns", not platform name.
#' @param medians If you specify "external" parameter as "Given.mdns", you
#'   should input matrix/table, 50 signatures in the first column and
#'   "Given.mdns" values in the second column.
#' @param Subtype Logic.
#' @param hasClinical Logic. Please specify if you prepared clinical
#'   information, like Tumore size as T column, lymphatic node status as NODE
#'   column.
#' @noRd
#'

makeCalls.parker <- function(mat,
    df.cln,
    calibration = c("None", "Internal", "External"),
    internal = NA,
    external = NA,
    medians = NULL,
    Subtype = FALSE,
    hasClinical = FALSE) {
    ## loading dataset
    data_env <- new.env(parent = emptyenv())
    data("BreastSubtypeRobj", envir = data_env, package = "BreastSubtypeR")
    BreastSubtypeRobj <- data_env[["BreastSubtypeRobj"]]


    calibration <- match.arg(calibration)

    fl.mdn <- BreastSubtypeRobj$medians

    if (calibration == "External" & external == "Given.mdns") {
        if (is.null(medians)) {
            stop(
                "Please input prepared medians as requires. "
            )
        } else {
            if (!is.data.frame(medians)) {
                medians <- as.data.frame(medians)
            }
            colnames(medians) <- c("X", "Given.mdns")

            medians.all <- merge(fl.mdn, medians, by = "X")
            rownames(medians.all) <- medians.all$X
            medians.all <- medians.all[, -1]

            suffix <- external
        }
    } else {
        medians.all <- fl.mdn

        rownames(medians.all) <- medians.all$X
        medians.all <- medians.all[, -1]
    }


    ## run
    # normalization
    mat <- docalibration(mat,
        medians.all,
        calibration,
        internal = internal,
        external = external
    )

    out <- sspPredict(
        BreastSubtypeRobj$centroid,
        mat,
        std = FALSE,
        distm = "spearman",
        Subtype = Subtype
    )

    if (Subtype) {
        Int.sbs <- data.frame(
            PatientID = names(out$predictions),
            BS = out$predictions,
            BS.Subtype = out$predictions.FourSubtype,
            row.names = NULL
        )
    } else {
        Int.sbs <- data.frame(
            PatientID = names(out$predictions),
            BS = out$predictions,
            row.names = NULL
        )
    }
    out$dist.RORSubtype <- -1 * out$dist.RORSubtype
    out$distances <- -1 * out$distances


    ## calculate and grouping
    out.ROR <- RORgroup(out, df.cln, Subtype = Subtype, hasClinical = hasClinical)

    return(list(
        BS.all = Int.sbs,
        score.ROR = out.ROR,
        mdns = medians.all,
        outList = out
    ))
}


#### function to form an ER-balanced subset and derive its median

#' Function to form an ER-balanced subset and derive its median
#' @param mat gene expression matrix
#' @param df.cln clinical information table with PatientID
#' @param calibration The calibration method to use. Options are "None",
#'   "Internal", or "External". If "Internal" is selected, see the "internal"
#'   parameter for further details. If "External" is selected, see the
#'   "external" parameter.
#' @param internal Specify the strategy for internal calibration. Options are
#'   median-centered ("medianCtr", default), mean-centered ("meanCtr"), or
#'   quantile-centered ("qCtr").
#' @param external Specify the platform name (i.e., the column name) for
#'   external medians, which are calculated by the training cohort. If you want
#'   to use user-provided medians, set this parameter to "Given.mdns" and
#'   provide the medians via the "medians" parameter.
#' @param medians If "Given.mdns" is specified for the "external" parameter,
#'   input a matrix/table where the first column contains 50 genes and the
#'   second column contains the corresponding "Given.mdns" values.
#' @param Subtype Logic.
#' @param hasClinical Logic. Specify whether clinical information is included.
#'   For example, tumor size should be in the "T" column, and lymph node status
#' - "NODE": Lymph node status (0 = negative; \eqn{\ge 1}{>= 1} = positive).
#' @param seed An integer value is used to set the random seed.
#' @noRd


makeCalls_ihc <- function(mat,
    df.cln,
    calibration = "Internal",
    internal = "IHC.mdns",
    external = NA,
    medians = NA,
    Subtype = FALSE,
    hasClinical = FALSE,
    seed = 118) {
    ## loading dataset
    data_env <- new.env(parent = emptyenv())
    data("BreastSubtypeRobj", envir = data_env, package = "BreastSubtypeR")
    BreastSubtypeRobj <- data_env[["BreastSubtypeRobj"]]


    ERN.ihc <- df.cln[which(df.cln$ER == "ER-"), ]
    dim(ERN.ihc) # [1] 153   9

    ERP.ihc <- df.cln[which(df.cln$ER == "ER+"), ]
    dim(ERP.ihc) # [1] 559   9

    # seed = 118
    if (dim(ERN.ihc)[1] > dim(ERP.ihc)[1]) {
        temp <- ERP.ihc
        ERP.ihc <- ERN.ihc
        ERN.ihc <- temp
    }
    # set.seed(seed);
    # i = sample(dim(ERP.ihc)[1],dim(ERN.ihc)[1])
    # take equal number of ER+ and ER- samples
    withr::with_seed(seed, {
        i <- sample(dim(ERP.ihc)[1], dim(ERN.ihc)[1])
        # take equal number of ER+ and ER- samples
    })

    length(ERP.ihc$PatientID[i]) # ER positive samples
    length(ERN.ihc$PatientID) # ER negative samples


    mbal.ihc <- mat[, c(ERP.ihc$PatientID[i], ERN.ihc$PatientID)]

    dim(mbal.ihc)

    # Find median
    suffix <- getsuffix(calibration = calibration, internal, external)

    mdns <- apply(mbal.ihc, 1, median, na.rm = TRUE)
    # compute median of each row i.e gene
    mdns.df <- as.data.frame(mdns)
    df.mdns <- data.frame(X = rownames(mdns.df), mdns.ihc = mdns.df$mdns)
    # ER-blanced set based on IHC status alone---

    colnames(df.mdns) <- c("X", suffix)

    ## merge mdns
    fl.mdn <- BreastSubtypeRobj$medians

    medians.all <- merge(fl.mdn, df.mdns, by = "X")
    rownames(medians.all) <- medians.all$X
    medians.all <- medians.all[, -1]

    ## centroids
    centroids <- BreastSubtypeRobj$centroid # pam50_centroids.txt

    ## normalization
    mat <- docalibration(mat, medians.all, calibration, internal)

    out <- sspPredict(centroids,
        mat,
        std = FALSE,
        distm = "spearman",
        Subtype = Subtype
    )

    if (Subtype) {
        Int.sbs <- data.frame(
            PatientID = names(out$predictions),
            BS = out$predictions,
            BS.Subtype = out$predictions.FourSubtype,
            row.names = names(out$predictions)
        )
    } else {
        Int.sbs <- data.frame(
            PatientID = names(out$predictions),
            BS = out$predictions,
            row.names = names(out$predictions)
        )
    }

    out$dist.RORSubtype <- -1 * out$dist.RORSubtype
    out$distances <- -1 * out$distances

    ## calculate and grouping
    out.ROR <- RORgroup(
        out, df.cln,
        hasClinical = hasClinical,
        Subtype = Subtype
    )

    return(list(
        BS.all = Int.sbs,
        score.ROR = out.ROR,
        mdns = medians.all,
        outList = out
    ))
}



#' Function for iterative ER-balanced subset gene centering
#' @param mat gene expression matrix
#' @param df.cln clinical information table with PatientID and IHC column
#' @param iterative Times to do iterative ER balanced procedure with certain
#'   ratio.
#' @param ratio The options are either 1:1 or 54 (ER+) : 64 (ER-) (default). The
#'   latter was ER ratio used for UNC230 train cohort.
#' @param calibration The calibration method to use. Options are "None",
#'   "Internal", or "External". If "Internal" is selected, see the "internal"
#'   parameter for further details. If "External" is selected, see the
#'   "external" parameter.
#' @param internal Specify the strategy for internal calibration. Options are
#'   median-centered ("medianCtr", default), mean-centered ("meanCtr"), or
#'   quantile-centered ("qCtr").
#' @param external Specify the platform name (i.e., the column name) for
#'   external medians, which are calculated by the training cohort. If you want
#'   to use user-provided medians, set this parameter to "Given.mdns" and
#'   provide the medians via the "medians" parameter.
#' @param medians If "Given.mdns" is specified for the "external" parameter,
#'   input a matrix/table where the first column contains 50 genes and the
#'   second column contains the corresponding "Given.mdns" values.
#' @param Subtype Logic. If `TRUE`, the function predicts four subtypes by
#'   excluding the Normal-like subtype.
#' @param hasClinical Logical. If `TRUE`, the function uses clinical data from
#'   the `pheno` table. Required columns include:
#' @param seed An integer value is used to set the random seed.
#' @noRd

makeCalls_ihc.iterative <- function(mat,
    df.cln,
    iteration = 100,
    ratio = 54 / 64,
    calibration = "Internal",
    internal = "ER.mdns",
    external = NA,
    medians = NA,
    Subtype = FALSE,
    hasClinical = FALSE,
    seed = 118) {
    ## loading dataset
    data_env <- new.env(parent = emptyenv())
    data("BreastSubtypeRobj", envir = data_env, package = "BreastSubtypeR")
    BreastSubtypeRobj <- data_env[["BreastSubtypeRobj"]]


    # load the published centroids for classifcation
    centroids <- BreastSubtypeRobj$centroid # pam50_centroids.txt

    ## preprocess the input matrix
    ### get ER- samples
    ERN.ihc <- df.cln[which(df.cln$ER == "ER-"), ]
    dim(ERN.ihc)

    ### get ER+ samples
    ERP.ihc <- df.cln[which(df.cln$ER == "ER+"), ]
    dim(ERP.ihc)

    ## check the ER composition
    if (dim(ERP.ihc)[1] < dim(ERN.ihc)[1]) {
        temp <- ERP.ihc
        ERP.ihc <- ERN.ihc
        ERN.ihc <- temp
    }

    ## check ratio
    ## make sure a reasonable ratio
    if (ratio > dim(ERP.ihc)[1] / dim(ERN.ihc)[1]) {
        n_ERP <- dim(ERP.ihc)[1]
        n_ERN <- dim(ERN.ihc)[1]
        stop(
            "please specify a ratio less than ",
            max(n_ERP, n_ERN) / min(n_ERP, n_ERN)
        )
    }

    withr::with_seed(seed, {
        res_ihc_iterative <- mapply(
            function(itr) {
                i <- sample(dim(ERP.ihc)[1], ceiling(ratio * dim(ERN.ihc)[1]))

                length(ERP.ihc$PatientID[i]) # ER positive samples
                length(ERN.ihc$PatientID) # ER negative samples

                mbal.ihc <- mat[, c(ERP.ihc$PatientID[i], ERN.ihc$PatientID)]

                suffix <- getsuffix(calibration = calibration, internal)

                # Calculate median
                mdns <- apply(mbal.ihc, 1, median, na.rm = TRUE)
                # compute median of each row i.e gene
                mdns.df <- as.data.frame(mdns)
                df.mdns <- data.frame(
                    X = rownames(mdns.df),
                    mdns.ihc = mdns.df$mdns
                )
                # ER-blanced set based on IHC status alone---
                colnames(df.mdns) <- c("X", suffix)

                ## integrate ihc.mdns
                fl.mdn <- BreastSubtypeRobj$medians

                medians.all <- merge(fl.mdn, df.mdns, by = "X")
                rownames(medians.all) <- medians.all$X
                medians.all <- medians.all[, -1]

                ## normalization
                mat <- docalibration(mat, medians.all, calibration, internal)

                out <- sspPredict(centroids,
                    mat,
                    std = FALSE,
                    distm = "spearman",
                    Subtype = Subtype
                )

                return(out)
            },
            paste0("itr.", seq(iteration)),
            SIMPLIFY = FALSE,
            USE.NAMES = TRUE
        )
    })

    ## get consensus intrinsic subtype
    Call_subtypes <- mapply(
        function(res_ihs) {
            res_ihs$predictions
        },
        res_ihc_iterative,
        SIMPLIFY = TRUE,
        USE.NAMES = TRUE
    )
    consensus_subtypes <- apply(Call_subtypes, 1, get_consensus_subtype)

    if (Subtype) {
        Call_subtypes.Subtype <- mapply(
            function(res_ihs) {
                res_ihs$predictions.FourSubtype
            },
            res_ihc_iterative,
            SIMPLIFY = TRUE,
            USE.NAMES = TRUE
        )
        cs_FourSubtype <- apply(Call_subtypes.Subtype, 1, get_consensus_subtype)

        Int.sbs <- data.frame(
            PatientID = names(consensus_subtypes),
            BS = consensus_subtypes,
            BS.Subtype = cs_FourSubtype,
            row.names = names(consensus_subtypes)
        )
    } else {
        Int.sbs <- data.frame(
            PatientID = names(consensus_subtypes),
            BS = consensus_subtypes,
            row.names = names(consensus_subtypes)
        )
    }

    mean_eve <- get_average_subtype(res_ihc_iterative, consensus_subtypes)

    if (Subtype) {
        out <- list(
            predictions = consensus_subtypes,
            predictions.FourSubtype = cs_FourSubtype,
            testData = mean_eve$testdata,
            distances = -1 * mean_eve$mean_distance,
            dist.RORSubtype = -1 * mean_eve$mean_distance.Subtype,
            centroids = centroids
        )
    } else {
        out <- list(
            predictions = consensus_subtypes,
            testData = mean_eve$testdata,
            distances = -1 * mean_eve$mean_distance,
            dist.RORSubtype = -1 * mean_eve$mean_distance.Subtype,
            centroids = centroids
        )
    }

    ## calculate and grouping
    out.ROR <- RORgroup(
        out, df.cln,
        hasClinical = hasClinical,
        Subtype = Subtype
    )

    if (Subtype) {
        res <- list(
            BS.all = Int.sbs,
            score.ROR = out.ROR,
            outList = out,
            BS.itr.keep = Call_subtypes,
            BS.itr.keep.Subtype = Call_subtypes.Subtype
        )
    } else {
        res <- list(
            BS.all = Int.sbs,
            score.ROR = out.ROR,
            outList = out,
            BS.itr.keep = Call_subtypes
        )
    }

    return(res)
}


#### form secondary ER-balanced set (refer to paper) leveraging PCA and
#### subsequent intermediate intrinsic subtypes

#' Function for the first step of PCA-PAM50 approach
#' @param mat gene expression matrix
#' @param df.cln clinical information table
#' @param calibration The calibration method to use, "Internal".
#' @param internal Specify the strategy for internal calibration, "PC1ihc.mdns".
#' @param external NA
#' @param medians NA
#' @param Subtype Logic
#' @param hasClinical Logic. Specify whether clinical information is included.
#'   For example, tumor size should be in the "T" column, and lymph node status
#' - "NODE": Lymph node status (0 = negative; \eqn{\ge 1}{>= 1} = positive).
#' @param seed An integer value is used to set the random seed.
#' @noRd

makeCalls.PC1ihc <- function(
        mat,
        df.cln,
        calibration = "Internal",
        internal = "PC1ihc.mdns",
        external = NA,
        medians = NA,
        Subtype = FALSE,
        hasClinical = FALSE,
        seed = 118) {
    ## loading dataset
    data_env <- new.env(parent = emptyenv())
    data("BreastSubtypeRobj", envir = data_env, package = "BreastSubtypeR")
    BreastSubtypeRobj <- data_env[["BreastSubtypeRobj"]]


    # Initial checks for 'df.cln' and 'mat'
    if (is.null(df.cln) ||
        !"PatientID" %in% colnames(df.cln) ||
        !"IHC" %in% colnames(df.cln)) {
        stop(
            "Clinical data 'df.cln' is missing or
            does not contain required columns 'PatientID' and 'IHC'.
            Refer to the vignette to create 'test.clinical' data."
        )
    }

    if (is.null(mat) || !is.matrix(mat)) {
        stop(
            "Gene expression matrix 'mat' is missing or not correctly formatted.
            Refer to the vignette to create 'test.matrix' data."
        )
    }

    # Pull the PCA components
    # rv      = rowVars(mat)
    # select  = order(rv, decreasing = TRUE)[seq_len(dim(mat)[1])]
    # the input is PAM50 matrix --50 genes -- get from dimension
    pca <- prcomp(t(mat)) # [select,]
    pc12 <- pca$x[, seq_len(2)] # get two principal
    df.pc1 <- data.frame(
        PatientID = rownames(pc12),
        PC1 = pc12[, 1],
        stringsAsFactors = FALSE
    )
    df.pca1 <- merge(df.cln, df.pc1, by = "PatientID")


    #--our function works best if majority of ER- cases
    # fall in the positive PC1 axis--check
    # Identify ER-negative cases
    er_negative <- !grepl("^L", df.pca1$IHC)

    # Determine if the majority of ER-negative cases fall in the negative axis
    # of PC1

    if (sum(df.pca1$PC1[er_negative] < 0) > sum(df.pca1$PC1[er_negative] >= 0)) {
        df.pca1$PC1 <- -df.pca1$PC1
    }

    # Ensure that IHC is not a factor or has all necessary levels defined
    if (is.factor(df.pca1$IHC)) {
        df.pca1$IHC <- as.character(df.pca1$IHC)
    }

    # Convert IHC column to uppercase to handle case insensitivity
    df.pca1$IHC <- toupper(df.pca1$IHC)

    # Function to count the number of misclassified cases
    # on a given PC1 point ---find the cutoff
    getno <- function(x) {
        p.rgt <- length(which(grepl("^L", df.pca1$IHC) &
            df.pca1$PC1 > x)) / length(which(grepl("^L", df.pca1$IHC)))
        n.lft <- length(which(!grepl("^L", df.pca1$IHC) &
            df.pca1$PC1 < x)) / length(which(!grepl("^L", df.pca1$IHC)))
        tot <- (p.rgt + n.lft) * 100
        return(list(PC1 = x, Mis = tot))
    }

    df.mis <- do.call(rbind.data.frame, lapply(seq(-20, 20, by = 0.1), getno))

    num.min <- df.mis$PC1[which(df.mis$Mis == min(df.mis$Mis))]

    # used mean to overcome situation where there are two minimum
    ERP.pc1ihc <- df.pca1[which(grepl("^L", df.pca1$IHC) &
        df.pca1$PC1 <= mean(num.min)), ]
    ERN.pc1ihc <- df.pca1[which(!grepl("^L", df.pca1$IHC) &
        df.pca1$PC1 > mean(num.min)), ]

    # dim(ERP.pc1ihc)
    # dim(ERN.pc1ihc)

    if (dim(ERP.pc1ihc)[1] < dim(ERN.pc1ihc)[1]) {
        temp <- ERN.pc1ihc
        ERN.pc1ihc <- ERP.pc1ihc
        ERP.pc1ihc <- temp
        rm(temp)
    }

    # set.seed(seed)
    # i = sample(dim(ERP.pc1ihc)[1], dim(ERN.pc1ihc)[1])
    # take equal number of ER+ and ER- samples
    withr::with_seed(seed, {
        i <- sample(dim(ERP.pc1ihc)[1], dim(ERN.pc1ihc)[1])
    })
    length(ERP.pc1ihc$PatientID[i])
    # ER positive samples
    length(ERN.pc1ihc$PatientID)
    # ER negative samples

    # subset pam50 matrix with the IDs corresponding to balanced ER+ and ER-
    mbal.pc1ihc <- mat[, c(ERP.pc1ihc$PatientID[i], ERN.pc1ihc$PatientID)]

    dim(mbal.pc1ihc)

    # Find median

    suffix <- getsuffix(calibration = calibration, internal)

    # compute median of each row i.e gene
    mdns <- apply(mbal.pc1ihc, 1, median, na.rm = TRUE)
    mdns.df <- as.data.frame(mdns)
    # ER-blanced set based on IHC status alone---
    df.mdns <- data.frame(X = rownames(mdns.df), mdns.pc1ihc = mdns.df$mdns)
    colnames(df.mdns) <- c("X", suffix)

    ## medians
    fl.mdn <- BreastSubtypeRobj$medians

    medians.all <- merge(fl.mdn, df.mdns, by = "X")
    rownames(medians.all) <- medians.all$X
    medians.all <- medians.all[, -1]


    # normalization
    mat <- docalibration(mat, medians.all, calibration, internal)

    out <- sspPredict(
        BreastSubtypeRobj$centroid,
        mat,
        std = FALSE,
        distm = "spearman",
        Subtype = Subtype
    )

    if (Subtype) {
        Int.sbs <- data.frame(
            PatientID = names(out$predictions),
            BS = out$predictions,
            BS.Subtype = out$predictions.FourSubtype,
            row.names = names(out$predictions)
        )
    } else {
        Int.sbs <- data.frame(
            PatientID = names(out$predictions),
            BS = out$predictions,
            row.names = names(out$predictions)
        )
    }

    out$dist.RORSubtype <- -1 * out$dist.RORSubtype
    out$distances <- -1 * out$distances

    ## calculate and grouping
    out.ROR <- RORgroup(
        out,
        df.cln,
        hasClinical = hasClinical,
        Subtype = Subtype
    )

    return(list(
        BS.all = Int.sbs,
        score.ROR = out.ROR,
        mdns = medians.all,
        outList = out
    ))
}

#' Function for the second step of PCA-PAM50 approach
#' @param mat gene expression matrix
#' @param df.pam clinical information table created using makeCalls.PC1ihc().
#' @param calibration The calibration method to use, "Internal".
#' @param internal Specify the strategy for internal calibration, "v1PAM.mdns".
#' @param external NA
#' @param medians NA
#' @param Subtype Logic.
#' @param hasClinical Logic. Specify whether clinical information is included.
#'   For example, tumor size should be in the "T" column, and lymph node status
#' - "NODE": Lymph node status (0 = negative; \eqn{\ge 1}{>= 1} = positive).
#' @param seed An integer value is used to set the random seed.
#' @noRd

makeCalls.v1PAM <- function(
        mat,
        df.pam,
        calibration = "Internal",
        internal = "v1PAM.mdns",
        external = NA,
        medians = NA,
        Subtype = FALSE,
        hasClinical = FALSE,
        seed = 118) {
    ## loading dataset
    data_env <- new.env(parent = emptyenv())
    data("BreastSubtypeRobj", envir = data_env, package = "BreastSubtypeR")
    BreastSubtypeRobj <- data_env[["BreastSubtypeRobj"]]


    ERN.pam <- df.pam[which(df.pam$PAM50 %in% c("Basal")), ]
    dim(ERN.pam)

    ERP.pam <- df.pam[which(df.pam$PAM50 %in% c("LumA")), ]
    dim(ERP.pam)

    # Determine the smaller size between ER+ and ER-
    sample_size <- min(dim(ERP.pam)[1], dim(ERN.pam)[1])

    # set.seed(seed)
    # i = sample(dim(ERP.pam)[1], dim(ERN.pam)[1])
    # take equal number of ER+ and ER- samples
    withr::with_seed(seed, {
        # i <- sample(dim(ERP.pam)[1], dim(ERN.pam)[1])
        if (dim(ERN.pam)[1] <= dim(ERP.pam)[1]) {
            i <- sample(dim(ERP.pam)[1], sample_size)
            ERP.pam_sampled <- ERP.pam[i, ]
            ERN.pam_sampled <- ERN.pam
        } else {
            i <- sample(dim(ERN.pam)[1], sample_size)
            ERN.pam_sampled <- ERN.pam[i, ]
            ERP.pam_sampled <- ERP.pam
        }
    })

    length(ERP.pam$PatientID[i])
    # ER positive samples
    length(ERN.pam$PatientID)
    # ER negative samples

    mbal.pam <- mat[, c(ERP.pam_sampled$PatientID, ERN.pam_sampled$PatientID)]
    dim(mbal.pam)

    # Find median
    suffix <- getsuffix(calibration = calibration, internal)

    mdns <- apply(mbal.pam, 1, median, na.rm = TRUE)
    # compute median of each row i.e gene
    mdns.df <- as.data.frame(mdns)
    df.mdns <- data.frame(X = rownames(mdns.df), mdns.pam = mdns.df$mdns)
    # ER-blanced set based on IHC status alone---
    colnames(df.mdns) <- c("X", suffix)

    ## median
    fl.mdn <- BreastSubtypeRobj$medians

    medians.all <- merge(fl.mdn, df.mdns, by = "X")
    rownames(medians.all) <- medians.all$X
    medians.all <- medians.all[, -1]


    ## normalization
    mat <- docalibration(mat, medians.all, calibration, internal)

    out <- sspPredict(
        BreastSubtypeRobj$centroid,
        mat,
        std = FALSE,
        distm = "spearman",
        Subtype = Subtype
    )

    if (Subtype) {
        Int.sbs <- data.frame(
            PatientID = names(out$predictions),
            BS = out$predictions,
            BS.Subtype = out$predictions.FourSubtype,
            row.names = names(out$predictions)
        )
    } else {
        Int.sbs <- data.frame(
            PatientID = names(out$predictions),
            BS = out$predictions,
            row.names = names(out$predictions)
        )
    }
    out$dist.RORSubtype <- -1 * out$dist.RORSubtype
    out$distances <- -1 * out$distances

    ## calculate and grouping
    out.ROR <- RORgroup(out,
        df.cln = df.pam,
        hasClinical = hasClinical,
        Subtype = Subtype
    )

    return(list(
        BS.all = Int.sbs,
        score.ROR = out.ROR,
        mdns.fl = medians.all,
        outList = out
    ))
}


#' Function for calling PAM50 subtypes by subgroup-specifc (ssBC) methods This
#' function is adapted from ssBC TNBC-BreastCancerRes2015 and subgroup specific
#' TNBC-JAMAOncol2024
#' @param mat gene expression matrix
#' @param df.cln clinical information table. The first column must be named
#'   "PatientID".
#' @param s Options are "ER" or "TN" or "ER.v2" or "TN.v2". Specify the medians
#'   you want. The original quantile is "ER" and "TN" of
#'   TNBC-BreastCancerRes2015.  If you choose "ER.v2" or "TN.v2", it means you
#'   choose quantile from TNBC-JAMAOncol2024.
#' @param Subtype Logic. Specify whether to predict Subtype-like subtyping.
#' @param hasClinical Logic. Specify whether clinical information is included.
#'   For example, tumor size should be in the "T" column, and lymph node status
#' - "NODE": Lymph node status (0 = negative; \eqn{\ge 1}{>= 1} = positive).
#' @noRd

makeCalls.ssBC <- function(
        mat,
        df.cln,
        s = c("ER", "TN", "ER.v2", "TN.v2"),
        Subtype = FALSE,
        hasClinical = FALSE) {
    ## loading dataset
    data_env <- new.env(parent = emptyenv())
    data("BreastSubtypeRobj", envir = data_env, package = "BreastSubtypeR")
    BreastSubtypeRobj <- data_env[["BreastSubtypeRobj"]]


    s <- match.arg(s)

    if (!(dim(mat)[2] == dim(df.cln)[1])) {
        stop(
            "Please input equal number of patient clinical information
            to the number of patient in gene expression matrix."
        )
    }

    gene.sigma <- BreastSubtypeRobj$ssBC.subgroupQuantile

    if (s == "ER") {
        ## use ER selected strategy

        ## if there is no sample in either of both, wont influence the code
        ERN_samples <- rownames(df.cln)[which(df.cln$ER == "ER-")]
        ERP_samples <- rownames(df.cln)[which(df.cln$ER == "ER+")]
        samples_selected <- list(ER_neg = ERN_samples, ER_pos = ERP_samples)
    } else if (s == "TN") {
        ## if there is no sample in either of both, wont influence the code
        TN_samples <- rownames(df.cln)[which(df.cln$TN == "TN")]
        nTN_samples <- rownames(df.cln)[which(df.cln$TN == "nonTN")]

        samples_selected <- list(TN = TN_samples, nonTN = nTN_samples)
    } else if (s == "ER.v2") {
        # TNBC-JAMAOncol2024

        ## if there is no sample in either of both, wont influence the code
        ERN_HER2N_samples <- rownames(df.cln)[which(df.cln$ER == "ER-" &
            df.cln$HER2 == "HER2-")]
        ERP_HER2N_samples <- rownames(df.cln)[which(df.cln$ER == "ER+" &
            df.cln$HER2 == "HER2-")]

        ERN_HER2P_samples <- rownames(df.cln)[which(df.cln$ER == "ER-" &
            df.cln$HER2 == "HER2+")]
        ERP_HER2P_samples <- rownames(df.cln)[which(df.cln$ER == "ER+" &
            df.cln$HER2 == "HER2+")]

        samples_selected <- list(
            ERneg_HER2neg = ERN_HER2N_samples,
            ERpos_HER2neg = ERP_HER2N_samples,
            HER2pos_ERneg = ERN_HER2P_samples,
            HER2pos_ERpos = ERP_HER2P_samples
        )
    } else if (s == "TN.v2") {
        ## selected cohort; TNBC-JAMAOncol2024
        ## if there is no sample in either of both, wont influence the code
        TN_samples <- rownames(df.cln)[which(df.cln$TN == "TN")]
        samples_selected <- list(TNBC = TN_samples)
    } else {
        stop("Please enter valid variable for s ")
    }

    res <- mapply(
        function(element) {
            # restrict to PAM50 genes first
            pam50 <- rownames(BreastSubtypeRobj$centroid)
            samp <- samples_selected[[element]]
            if (length(samp) == 0L) {
                # return a well-formed empty result for this subgroup
                pred <- setNames(character(0), character(0))
                empty <- list(
                    predictions = pred,
                    testData = matrix(numeric(0), nrow = length(pam50), ncol = 0, dimnames = list(pam50, NULL)),
                    distances = matrix(numeric(0),
                        nrow = 0, ncol = ncol(BreastSubtypeRobj$centroid),
                        dimnames = list(NULL, colnames(BreastSubtypeRobj$centroid))
                    ),
                    dist.RORSubtype = matrix(numeric(0),
                        nrow = 0, ncol = 4,
                        dimnames = list(NULL, colnames(BreastSubtypeRobj$centroid)[seq_len(4)])
                    )
                )
                if (Subtype) empty$predictions.FourSubtype <- pred
                return(empty)
            }

            x.m <- mat[rownames(mat) %in% pam50, samp, drop = FALSE]
            if (nrow(x.m) == 0L) {
                # no overlapping PAM50 genes in this subgroup
                pred <- setNames(character(0), character(0))
                empty <- list(
                    predictions = pred,
                    testData = matrix(numeric(0), nrow = length(pam50), ncol = 0, dimnames = list(pam50, NULL)),
                    distances = matrix(numeric(0),
                        nrow = 0, ncol = ncol(BreastSubtypeRobj$centroid),
                        dimnames = list(NULL, colnames(BreastSubtypeRobj$centroid))
                    ),
                    dist.RORSubtype = matrix(numeric(0),
                        nrow = 0, ncol = 4,
                        dimnames = list(NULL, colnames(BreastSubtypeRobj$centroid)[seq_len(4)])
                    )
                )
                if (Subtype) empty$predictions.FourSubtype <- pred
                return(empty)
            }

            # gene-specific quantiles: align safely
            gene.sigma.o <- gene.sigma[rownames(x.m), element, drop = TRUE]

            # per-row quantile with guards
            x.sigma <- vapply(seq_len(nrow(x.m)), function(i) {
                p <- gene.sigma.o[i]
                if (!is.finite(p)) {
                    return(NA_real_)
                }
                tryCatch(quantile(x.m[i, ], probs = p, na.rm = TRUE, names = FALSE, type = 7),
                    error = function(e) NA_real_
                )
            }, numeric(1))

            x.m <- sweep(x.m, 1, x.sigma, FUN = "-")

            sspPredict(
                BreastSubtypeRobj$centroid,
                x.m,
                std = FALSE,
                distm = "spearman",
                Subtype = Subtype
            )
        },
        names(samples_selected),
        SIMPLIFY = FALSE,
        USE.NAMES = TRUE
    )


    ## prepare data for ROR

    if (s == "ER") {
        ## use ER selected strategy

        predictions <- c(res$ER_neg$predictions, res$ER_pos$predictions)
        distances <- rbind(res$ER_neg$distances, res$ER_pos$distances)
        testData <- .align_and_cbind_by_genes(res$ER_neg$testData, res$ER_pos$testData)
        dist.RORSubtype <- rbind(
            res$ER_neg$dist.RORSubtype,
            res$ER_pos$dist.RORSubtype
        )

        if (Subtype) {
            predictions.FourSubtype <- c(
                res$ER_neg$predictions.FourSubtype,
                res$ER_pos$predictions.FourSubtype
            )
        }
    } else if (s == "ER.v2") {
        predictions <- c(
            res$ERneg_HER2neg$predictions,
            res$ERpos_HER2neg$predictions,
            res$HER2pos_ERneg$predictions,
            res$HER2pos_ERpos$predictions
        )
        distances <- rbind(
            res$ERneg_HER2neg$distances,
            res$ERpos_HER2neg$distances,
            res$HER2pos_ERneg$distances,
            res$HER2pos_ERpos$distances
        )
        testData <- .align_and_cbind_by_genes(
            res$ERneg_HER2neg$testData,
            res$ERpos_HER2neg$testData,
            res$HER2pos_ERneg$testData,
            res$HER2pos_ERpos$testData
        )
        dist.RORSubtype <- rbind(
            res$ERneg_HER2neg$dist.RORSubtype,
            res$ERpos_HER2neg$dist.RORSubtype,
            res$HER2pos_ERneg$dist.RORSubtype,
            res$HER2pos_ERpos$dist.RORSubtype
        )

        if (Subtype) {
            predictions.FourSubtype <- c(
                res$ERneg_HER2neg$predictions.FourSubtype,
                res$ERpos_HER2neg$predictions.FourSubtype,
                res$HER2pos_ERneg$predictions.FourSubtype,
                res$HER2pos_ERpos$predictions.FourSubtype
            )
        }
    } else if (s == "TN") {
        predictions <- c(
            res$TN$predictions,
            res$nonTN$predictions
        )
        distances <- rbind(
            res$TN$distances,
            res$nonTN$distances
        )
        testData <- .align_and_cbind_by_genes(
            res$TN$testData,
            res$nonTN$testData
        )
        dist.RORSubtype <- rbind(
            res$TN$dist.RORSubtype,
            res$nonTN$dist.RORSubtype
        )

        if (Subtype) {
            predictions.FourSubtype <- c(
                res$TN$predictions.FourSubtype,
                res$nonTN$predictions.FourSubtype
            )
        }
    } else if (s == "TN.v2") {
        predictions <- res$TNBC$predictions
        distances <- res$TNBC$distances
        testData <- res$TNBC$testData
        dist.RORSubtype <- res$TNBC$dist.RORSubtype

        if (Subtype) {
            predictions.FourSubtype <- res$TNBC$predictions.FourSubtype
        }
    }

    if (Subtype) {
        Int.sbs <- data.frame(
            PatientID = names(predictions),
            BS = predictions,
            BS.Subtype = predictions.FourSubtype,
            row.names = names(predictions),
            check.names = FALSE
        )

        ## out list
        out <- list(
            predictions = predictions,
            predictions.FourSubtype = predictions.FourSubtype,
            testData = testData,
            distances = distances,
            dist.RORSubtype = dist.RORSubtype,
            centroids = BreastSubtypeRobj$centroid
        )
    } else {
        ## subtype table
        Int.sbs <- data.frame(
            PatientID = names(predictions),
            BS = predictions,
            row.names = names(predictions)
        )

        ## out list
        out <- list(
            predictions = predictions,
            testData = testData,
            distances = distances,
            dist.RORSubtype = dist.RORSubtype,
            centroids = BreastSubtypeRobj$centroid
        )
    }

    out$dist.RORSubtype <- -1 * out$dist.RORSubtype
    out$distances <- -1 * out$distances

    ## calculate and grouping ROR score
    out.ROR <- RORgroup(
        out, df.cln,
        hasClinical = hasClinical,
        Subtype = Subtype
    )

    ## supplementing NA samples
    if (length(out$predictions) != ncol(mat)) {
        sample_ex <- setdiff(colnames(mat), Int.sbs$PatientID)

        predictions_ex <- rep(NA, length(sample_ex))
        names(predictions_ex) <- sample_ex
        out$predictions <- c(out$predictions, predictions_ex)

        testData_ex <- matrix(
            data = NA, nrow = nrow(testData), ncol = length(sample_ex),
            dimnames = list(rownames(testData), sample_ex)
        )
        out$testData <- cbind(out$testData, testData_ex)

        distances_ex <- matrix(
            data = NA, nrow = length(sample_ex), ncol = 5,
            dimnames = list(sample_ex, colnames(distances))
        )
        out$distances <- rbind(out$distances, distances_ex)

        dist.RORSubtype_ex <- matrix(
            data = NA, nrow = length(sample_ex), ncol = 4,
            dimnames = list(sample_ex, colnames(dist.RORSubtype))
        )
        out$dist.RORSubtype <- rbind(out$dist.RORSubtype, dist.RORSubtype_ex)

        out.ROR_ex <- matrix(
            data = NA, nrow = length(sample_ex), ncol = ncol(out.ROR),
            dimnames = list(sample_ex, colnames(out.ROR))
        )
        out.ROR <- rbind(out.ROR, out.ROR_ex)

        if (Subtype) {
            predictions.FourSubtype_ex <- rep(NA, length(sample_ex))
            names(predictions.FourSubtype_ex) <- sample_ex
            out$predictions.FourSubtype <- c(out$predictions.FourSubtype, predictions.FourSubtype_ex)

            Int.sbs_ex <- data.frame(
                PatientID = sample_ex,
                BS = NA,
                BS.Subtype = NA,
                row.names = sample_ex
            )
            Int.sbs <- rbind(Int.sbs, Int.sbs_ex)
        } else {
            Int.sbs_ex <- data.frame(
                PatientID = sample_ex,
                BS = NA,
                row.names = sample_ex
            )
            Int.sbs <- rbind(Int.sbs, Int.sbs_ex)
        }
    }

    ## reorder
    Int.sbs <- Int.sbs[colnames(mat), ]
    out.ROR <- out.ROR[colnames(mat), ]
    out$dist.RORSubtype <- out$dist.RORSubtype[colnames(mat), ]
    out$distances <- out$distances[colnames(mat), ]

    out$testData <- out$testData[, colnames(mat)]
    out$predictions <- out$predictions[colnames(mat)]

    return(list(
        BS.all = Int.sbs,
        score.ROR = out.ROR,
        mdns = gene.sigma,
        outList = out
    ))
}
