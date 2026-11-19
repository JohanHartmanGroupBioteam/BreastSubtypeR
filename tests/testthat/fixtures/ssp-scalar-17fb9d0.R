# Verbatim definitions from R/NC_functions.R at
# 17fb9d0358db41b5d7ed777a0c0434db0b4262db.
# Frozen scalar reference for the whole sspPredict() call.

overlapSets <- function(x, y) {
    # subset the two lists to have a commonly ordered gene list
    x <- x[dimnames(x)[[1]] %in% dimnames(y)[[1]], ]
    y <- y[dimnames(y)[[1]] %in% dimnames(x)[[1]], ]

    # and sort such that thing are in the correct order
    x <- x[sort.list(row.names(x)), ]
    y <- y[sort.list(row.names(y)), ]

    return(list(x = x, y = y))
}

standardize <- function(x) {
    annAll <- dimnames(x)
    x <- scale(x)
    dimnames(x) <- annAll
    return(x)
}

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
