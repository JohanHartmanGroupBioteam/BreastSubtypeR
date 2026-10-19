# Shared presentation settings; these helpers do not transform plot data.
.bsr_plot_theme <- function() {
    ggplot2::theme_classic(base_size = 11, base_family = "sans") +
        ggplot2::theme(
            axis.text = ggplot2::element_text(size = 10, colour = "#28323C"),
            axis.title = ggplot2::element_text(size = 11),
            plot.title = ggplot2::element_text(
                size = 13, face = "bold", margin = ggplot2::margin(b = 7)),
            plot.subtitle = ggplot2::element_text(
                size = 10, colour = "#4B5563", margin = ggplot2::margin(b = 10)),
            plot.caption = ggplot2::element_text(
                size = 9, colour = "#4B5563", hjust = 0,
                margin = ggplot2::margin(t = 12)),
            legend.title = ggplot2::element_text(size = 10, face = "bold"),
            legend.text = ggplot2::element_text(size = 10),
            legend.key.size = grid::unit(4, "mm"),
            legend.spacing.y = grid::unit(2, "mm"),
            plot.margin = ggplot2::margin(12, 14, 10, 12)
        )
}

.bsr_legend_style <- function(title, ...) {
    list(title = title,
         title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
         labels_gp = grid::gpar(fontsize = 10), border = "white", ...)
}

.bsr_missing_key <- function() {
    .bsr_legend_style("Missing calls", labels = "NA",
                      legend_gp = grid::gpar(fill = "grey", col = NA))
}

# Reserve a footer so an ordinary Heatmap object prints with the approved
# horizontal keys. Use supported constructors and a post-draw callback;
# do not mutate Heatmap slots or global draw options.
.bsr_legend_footer <- function(name, legends) {
    force(name)
    force(legends)
    annotation <- do.call(ComplexHeatmap::HeatmapAnnotation, c(
        setNames(list(ComplexHeatmap::anno_empty(
            border = FALSE, height = grid::unit(22, "mm"))), name),
        list(show_annotation_name = FALSE)
    ))
    draw_footer <- function(ht) {
        previous_viewport <- grid::current.viewport()$name
        on.exit({
            if (previous_viewport == "ROOT") grid::upViewport(0)
            else grid::seekViewport(previous_viewport)
        }, add = TRUE)
        # Keep the position relative to the device. Converting annotation
        # locations to inches freezes their positions during RStudio replay.
        # The footer reserves 22 mm above the default 2 mm plot padding.
        grid::upViewport(0)
        grid::pushViewport(grid::viewport(
            x = grid::unit(0.5, "npc"), y = grid::unit(13, "mm"),
            width = grid::unit(1, "npc") - grid::unit(4, "mm"),
            height = grid::unit(22, "mm")
        ))
        on.exit(grid::upViewport(), add = TRUE, after = FALSE)
        keys <- lapply(legends, function(x) do.call(ComplexHeatmap::Legend, x))
        # packLegend measures on its own temporary device, so resolve this
        # width on the active plot device before handing it to the packer.
        available_width <- grid::convertWidth(grid::unit(1, "npc"), "mm")
        packed <- ComplexHeatmap::packLegend(
            list = keys, direction = "horizontal",
            max_width = available_width, gap = grid::unit(4, "mm"))
        ComplexHeatmap::draw(packed, x = grid::unit(0.5, "npc"),
                            y = grid::unit(0.5, "npc"))
    }
    list(annotation = annotation, draw = draw_footer)
}
