# The site labels describe the expected main, four-subtype and ROR paths when
# Subtype is TRUE. A missing site is unreachable due to an error or early exit.
ssp_predict_cases <- function() {
    model <- new.env(parent = emptyenv())
    utils::data("BreastSubtypeRobj", package = "BreastSubtypeR", envir = model)
    centres <- model$BreastSubtypeRobj$centroid
    genes <- rownames(centres)
    excluded <- c("BIRC5", "CCNB1", "GRB7", "MYBL2")
    retained <- setdiff(genes, excluded)

    # A reproducible untied panel with nonidentical sample ranks. Generate it
    # without consuming RNG, so cases can also run with no .Random.seed.
    expression <- outer(seq_along(genes), seq_len(100L), function(i, j)
        sin(i * 1.7 + j * 0.31) + cos(i * j * 0.17))
    dimnames(expression) <- list(genes, paste0("sample", seq_len(ncol(expression))))
    small <- expression[, 1:4]

    batched <- c(main = "B", four = "B", ror = "B")
    scalar <- c(main = "S", four = "S", ror = "S")
    ror_only <- c(main = "S", four = "S", ror = "B")
    scalar_ror <- c(main = "B", four = "B", ror = "S")
    fifth_only <- c(main = "S", four = "B", ror = "B")
    spec <- function(x = centres, y = small, std = FALSE, distm = "spearman",
                     sites = batched) {
        list(args = list(x = x, y = y, std = std, distm = distm), sites = sites)
    }
    change_sample <- function(rows, value) {
        sample <- small
        sample[rows, 1L] <- value
        sample
    }
    change_centre <- function(rows, column, value) {
        x <- centres
        x[rows, column] <- value
        x
    }

    ties <- round(2 * small)
    storage.mode(ties) <- "integer"
    tied_centres <- as.matrix(round(4 * centres))
    storage.mode(tied_centres) <- "integer"
    reversed <- small[rev(genes), c(4L, 1L)]
    unnamed <- duplicate_names <- unusual_names <- reversed
    colnames(unnamed) <- NULL
    colnames(duplicate_names) <- rep("replicate", 2L)
    colnames(unusual_names) <- c(NA_character_, "")
    constant <- change_sample(genes, 1)
    ror_constant <- change_sample(retained, 1)
    one_ror_gene <- c(excluded, retained[1L])
    no_overlap <- as.data.frame(small)
    rownames(no_overlap) <- paste0("other", seq_along(genes))

    duration <- small
    class(duration) <- "difftime"
    attr(duration, "units") <- "secs"
    # These ranks have exactly zero correlation. Negation in the admitted
    # path produces -0, observable with the strict full-object comparator.
    zero_centres <- data.frame(first = 1:4, second = 4:1, row.names = retained[1:4])
    zero_samples <- cbind(a = c(2, 4, 1, 3), b = 1:4)
    rownames(zero_samples) <- retained[1:4]
    nested_matrix <- nested_frame <- centres
    nested_matrix[[1L]] <- I(cbind(centres[[1L]], centres[[2L]]))
    nested_frame[[1L]] <- centres[, 1:2, drop = FALSE]

    list(
        complete = spec(y = expression),
        integer_ties = spec(x = tied_centres, y = ties),
        partial_panel = spec(y = small[seq_len(37L), ]),
        four_classes = spec(x = centres[, 1:4]),
        two_classes = spec(x = as.matrix(centres[, 1:2])),
        zero_distance = spec(x = zero_centres, y = zero_samples),
        reversed = spec(y = reversed),
        unnamed = spec(y = unnamed),
        duplicate_names = spec(y = duplicate_names),
        missing_and_empty_names = spec(y = unusual_names),
        all_constant_sample = spec(y = constant, sites = scalar),
        ror_constant_sample = spec(y = ror_constant, sites = scalar_ror),
        excluded_na = spec(y = change_sample(excluded[1L], NA_real_), sites = ror_only),
        retained_na = spec(y = change_sample(retained[1L], NA_real_), sites = scalar),
        retained_nan = spec(y = change_sample(retained[1L], NaN), sites = scalar),
        retained_inf = spec(y = change_sample(retained[1:2], c(Inf, -Inf)), sites = scalar),
        excluded_inf = spec(y = change_sample(excluded[1:2], c(Inf, -Inf)), sites = ror_only),
        standardized = spec(std = TRUE),
        standardized_constant = spec(y = constant, std = TRUE, sites = scalar),
        standardized_overflow = spec(y = change_sample(genes, small[, 1L] * 1e200),
            std = TRUE, sites = ror_only),
        pearson = spec(distm = "pearson", sites = ror_only),
        correlation = spec(distm = "correlation", sites = ror_only),
        euclidean = spec(distm = "euclidean", sites = c(main = "E", four = "E", ror = "B")),
        euclidean_constant = spec(y = constant, distm = "euclidean",
            sites = c(main = "E", four = "E", ror = "S")),
        pearson_constant = spec(y = constant, distm = "pearson", sites = scalar),
        class_constant = spec(x = change_centre(genes, 2L, 0), sites = scalar),
        class_na = spec(x = change_centre(retained[1L], 2L, NA_real_), sites = scalar),
        fifth_constant = spec(x = change_centre(genes, 5L, 1), sites = fifth_only),
        fifth_na = spec(x = change_centre(retained[1L], 5L, NA_real_), sites = fifth_only),
        excluded_class_na = spec(x = change_centre(excluded[1L], 2L, NA_real_), sites = ror_only),
        frame_input = spec(y = as.data.frame(small), sites = scalar),
        duplicated_gene = spec(y = rbind(small, small[1L, , drop = FALSE]),
            sites = c(main = "S")),
        one_sample = spec(y = small[, 1L, drop = FALSE], sites = character()),
        one_gene_matrix = spec(y = small[retained[1L], , drop = FALSE], sites = character()),
        one_gene_frame = spec(y = as.data.frame(small)[retained[1L], , drop = FALSE],
            sites = scalar),
        late_alignment_error = spec(y = constant[one_ror_gene, ], sites = scalar[1:2]),
        late_alignment_frame = spec(y = as.data.frame(constant[one_ror_gene, ]), sites = scalar),
        empty_overlap = spec(y = no_overlap, sites = character()),
        zero_samples = spec(y = small[, integer(), drop = FALSE], sites = character()),
        classed_matrix = spec(y = duration, sites = scalar),
        matrix_column = spec(x = nested_matrix, sites = c(main = "S")),
        frame_column = spec(x = nested_frame, sites = c(main = "S")),
        signed_zero = spec(y = change_sample(genes, rep(c(-0, 0), length.out = length(genes))),
            sites = scalar)
    )
}
