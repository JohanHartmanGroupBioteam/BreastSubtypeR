visual_example <- function() {
    env <- new.env()
    data("OSLO2EMIT0obj", package = "BreastSubtypeR", envir = env)
    obj <- env$OSLO2EMIT0obj
    calls <- obj$res$results$genefu.robust$BS.all
    list(obj = obj, x = SummarizedExperiment::assay(obj$data_input$se_NC),
         out = data.frame(PatientID = calls$PatientID, Subtype = calls$BS))
}

visual_draw <- function(p) {
    path <- tempfile(fileext = ".pdf")
    grDevices::pdf(path, width = 10, height = 7.4)
    on.exit({ grDevices::dev.off(); unlink(path) })
    ComplexHeatmap::draw(p)@ht_list[[1L]]
}

visual_reference_heatmap <- function(x, out) {
    # Scientific specification before the styling refresh: no styling helper.
    z <- t(scale(t(x[, out$PatientID])))
    ComplexHeatmap::Heatmap(
        z, clustering_distance_rows = "spearman",
        clustering_method_rows = "average", show_row_dend = FALSE,
        cluster_column_slices = TRUE, column_split = out$Subtype,
        clustering_distance_columns = "spearman",
        clustering_method_columns = "average", show_column_names = FALSE
    )
}

expect_visual_tree_identity <- function(reference, candidate) {
    for (field in c("matrix", "row_order", "column_order", "row_order_list",
                    "column_order_list", "row_dend_list", "column_dend_list",
                    "row_dend_slice", "column_dend_slice")) {
        expect_identical(methods::slot(candidate, field),
                         methods::slot(reference, field), info = field)
    }
    expect_identical(candidate@matrix_param$column_split,
                     reference@matrix_param$column_split)
    expect_identical(candidate@matrix_param$row_split,
                     reference@matrix_param$row_split)
}

test_that("heatmap styling preserves scaled values, trees, order and NA slices", {
    ex <- visual_example()
    before <- serialize(ex, NULL)
    colors <- c(Basal = "red", Her2 = "hotpink", LumA = "darkblue",
                LumB = "skyblue", Normal = "green")
    out_na <- ex$out
    out_na$Subtype <- ex$obj$res$results$ssBC.v2$BS.all$BS
    expect_identical(sum(is.na(out_na$Subtype)), 8L)
    for (out in list(ex$out, out_na)) {
        p <- Vis_heatmap(ex$x, out)
        expect_s4_class(p, "Heatmap")
        reference <- visual_draw(visual_reference_heatmap(ex$x, out))
        candidate <- visual_draw(p)
        expect_visual_tree_identity(reference, candidate)
        z <- t(scale(t(ex$x[, out$PatientID])))
        expect_identical(p@matrix, z)
        anchors <- c(min(z), 0, max(z))
        expected_map <- circlize::colorRamp2(
            anchors, c("#2166AC", "#FFFFFF", "#B2182B"))
        expect_identical(attr(p@matrix_color_mapping@col_fun, "breaks"), anchors)
        expect_identical(p@matrix_color_mapping@col_fun(as.vector(z)),
                         expected_map(as.vector(z)))
        annotation <- p@top_annotation@anno_list$Subtype
        expect_identical(get("value", annotation@fun@var_env),
                         setNames(out$Subtype, out$PatientID))
        expected_annotation <- ComplexHeatmap::HeatmapAnnotation(
            Subtype = out$Subtype, col = list(Subtype = colors))
        expect_identical(annotation@color_mapping@colors,
                         expected_annotation@anno_list$Subtype@color_mapping@colors)
        expect_identical(annotation@color_mapping@na_col, "#BEBEBEFF")
        expect_true(p@row_dend_param$show)
        expect_identical(p@row_dend_param$distance, "spearman")
        expect_identical(p@row_dend_param$method, "average")
    }
    expect_identical(serialize(ex, NULL), before)
})

test_that("multi-method plot keeps the two sorts, calls, categories and entropy", {
    ex <- visual_example()
    data <- ex$obj$res$res_subtypes
    expected <- data[order(data[, ncol(data) - 1L]), ]
    expected <- expected[order(expected[, ncol(expected)], decreasing = FALSE), ]
    p <- Vis_Multi(data)
    expect_s4_class(p, "Heatmap")
    expect_identical(p@matrix, t(as.matrix(expected[, -ncol(expected)])))
    expect_identical(sum(is.na(p@matrix)), 8L)
    expect_identical(as.vector(get("value", p@top_annotation@anno_list$Entropy@fun@var_env)),
                     expected$entropy)
    h <- visual_draw(p)
    expect_identical(h@column_order, seq_len(nrow(data)))
    expect_identical(h@row_order, seq_len(ncol(data) - 1L))
    category <- get("value", p@right_annotation@anno_list$Category@fun@var_env)
    expect_identical(as.character(category), c(rep("NC-based", 6), rep("SSP-based", 2)))
    expect_identical(p@matrix_color_mapping@na_col, "#BEBEBEFF")
    expect_identical(p@matrix_color_mapping@colors,
        c(Basal = "#FF0000FF", Her2 = "#FF69B4FF", LumA = "#00008BFF",
          LumB = "#87CEEBFF", Normal = "#00FF00FF"))
})

test_that("pie counts, denominators, rounding and polar geometry are preserved", {
    ex <- visual_example()
    out_na <- ex$out
    out_na$Subtype <- ex$obj$res$results$ssBC.v2$BS.all$BS
    for (out in list(ex$out, out_na)) {
        p <- Vis_pie(out)
        expect_s3_class(p, "ggplot")
        tab <- data.frame(table(out$Subtype))
        expect_identical(p$data$Freq, tab$Freq)
        expect_identical(p$data$Var1, tab$Var1)
        expect_identical(p$data$perc, round(tab$Freq / sum(tab$Freq) * 100, 2))
        expect_identical(sum(p$data$Freq), sum(!is.na(out$Subtype)))
        expect_match(p$labels$caption,
            sprintf("%d assigned calls; %d missing calls excluded", sum(tab$Freq),
                    sum(is.na(out$Subtype))), fixed = TRUE)
        expect_identical(p$scales$get_scales("fill")$labels,
                         sprintf("%s   %d (%.2f%%)", tab$Var1, tab$Freq, p$data$perc))
        reference <- ggplot2::ggplot(p$data, ggplot2::aes(x = "", y = perc, fill = Var1)) +
            ggplot2::geom_bar(stat = "identity") +
            ggplot2::scale_fill_manual(values = c(Basal = "red", Her2 = "hotpink",
                LumA = "darkblue", LumB = "skyblue", Normal = "green")) +
            ggplot2::coord_polar("y", direction = -1)
        expect_identical(ggplot2::ggplot_build(p)$data[[1]],
                         ggplot2::ggplot_build(reference)$data[[1]])
        expect_identical(p$coordinates$direction, -1)
    }
})

test_that("PCA and neutral boxes retain their numerical quantities", {
    ex <- visual_example()
    p <- Vis_PCA(ex$x, ex$out)
    pca <- stats::prcomp(t(ex$x), center = TRUE, scale. = TRUE)
    expect_identical(p$data[, colnames(pca$x)], as.data.frame(pca$x))
    expect_identical(p$data$PatientID, rownames(pca$x))
    expect_identical(p$data$Subtype, ex$out$Subtype)
    scree <- Vis_PCA(ex$x, ex$out, Eigen = TRUE)
    expect_s3_class(scree, "ggplot")
    expect_identical(scree$data$Variance, (pca$sdev^2 / sum(pca$sdev^2) * 100)[1:10])
    correlations <- ex$obj$res$results$genefu.robust$outList$distances
    box <- Vis_boxplot(ex$out, correlations)
    expected <- data.frame(predictions = ex$out$Subtype,
                           cor = apply(correlations, 1, max))
    reference <- ggplot2::ggplot(expected, ggplot2::aes(x = predictions, y = cor)) +
        ggplot2::geom_boxplot()
    expect_identical(box$data, expected)
    expect_identical(ggplot2::ggplot_build(box)$data,
                     ggplot2::ggplot_build(reference)$data)
    expect_identical(box$labels$y, "Maximum centroid correlation")
    expect_null(box$mapping$fill)
})

test_that("true four-class API results plot without relabeling five-class calls", {
    ex <- visual_example()
    withr::local_seed(20260831)
    withr::local_options(BreastSubtypeR.verbose = FALSE)
    fresh <- suppressMessages(BS_Multi(ex$obj$data_input,
        methods = c("genefu.robust", "AIMS"), Subtype = TRUE, hasClinical = FALSE))
    calls <- fresh$results$genefu.robust$BS.all
    expect_true("BS.Subtype" %in% names(calls))
    expect_true(any(calls$BS == "Normal"))
    expect_true(all(calls$BS.Subtype %in% c("Basal", "Her2", "LumA", "LumB")))
    out <- data.frame(PatientID = calls$PatientID, Subtype = calls$BS.Subtype)
    h <- visual_draw(Vis_heatmap(ex$x, out))
    expect_visual_tree_identity(visual_draw(visual_reference_heatmap(ex$x, out)), h)
    pie <- Vis_pie(out)
    expect_identical(sum(pie$data$Freq), nrow(out))
    expect_false("Normal" %in% as.character(pie$data$Var1))
    multi <- Vis_Multi(fresh$res_subtypes.Subtype)
    expect_false("Normal" %in% as.vector(multi@matrix))
    expect_identical(sum(is.na(multi@matrix)), sum(is.na(fresh$res_subtypes.Subtype)))
    expect_s4_class(visual_draw(multi), "Heatmap")
})

test_that("existing all-no-call and unusual-input behavior is not repaired", {
    ex <- visual_example()
    none <- ex$out
    none$Subtype[] <- NA_character_
    expect_s4_class(visual_draw(Vis_heatmap(ex$x, none)), "Heatmap")
    no_methods <- ex$obj$res$res_subtypes
    no_methods[, -ncol(no_methods)] <- NA_character_
    no_methods$entropy <- NA_real_
    h <- visual_draw(Vis_Multi(no_methods))
    expect_true(all(is.na(h@matrix)))
    expect_true(all(is.na(get("value", h@top_annotation@anno_list$Entropy@fun@var_env))))
    # These are pre-existing failures, not invitations to change the policy.
    expect_error(ggplot2::ggplot_build(Vis_pie(none)), "Var1")
    constant <- ex$x
    constant[1, ] <- 1
    expect_error(Vis_PCA(constant, ex$out), "constant/zero")
    expect_error(Vis_heatmap(constant, ex$out))
    unknown <- ex$out
    unknown$Subtype[1] <- "unknown"
    expect_error(visual_draw(Vis_heatmap(ex$x, unknown)), "cannot map colors")
    shuffled <- ex$out[rev(seq_len(nrow(ex$out))), ]
    correlations <- ex$obj$res$results$genefu.robust$outList$distances
    expect_identical(Vis_boxplot(shuffled, correlations)$data$predictions,
                     shuffled$Subtype)
    expect_identical(Vis_boxplot(shuffled, correlations)$data$cor,
                     unname(apply(correlations, 1, max)))
})
