# Calibration choices are available in the initial HTML and on each rebind.
.shiny_pam50_defaults <- list(calibration = "Internal", internal = "medianCtr",
                             external = "Given.mdns")

.shiny_pam50_controls <- function(ns = shiny::NS(NULL), defaults = .shiny_pam50_defaults,
                                  allow_custom = TRUE) {
  bslib::card(
    style = "overflow: visible;",
    bslib::layout_columns(
      col_widths = c(4, 8),

      shiny::selectizeInput(
        ns("calibration"), "Calibration strategy",
        choices = list("None" = "None", "External" = "External", "Internal" = "Internal"),
        selected = defaults$calibration, width = "100%",
        options = list(openOnFocus = TRUE, dropdownParent = "body")
      ),

      shiny::div(
        style = "position: relative; overflow: visible;",

        shiny::conditionalPanel(
          ns = ns,
          condition = "input.calibration == 'External'",
          shiny::tagList(
            shiny::selectizeInput(
              ns("external"),
              "Reference medians",
              choices = c(
                if (allow_custom) stats::setNames("Given.mdns", "Custom (upload file\u2026)"),
                "Built-in: nCounter"                                     = "nCounter",
                "Built-in: RNA-seq (Freeze 2012-09-07)"                  = "RNAseq.Freeze.20120907",
                "Built-in: totalRNA FFPE (2015-11-11)"                   = "totalRNA.FFPE.20151111",
                "Built-in: RNA-seq V2"                                   = "RNAseq.V2",
                "Built-in: RNA-seq V1"                                   = "RNAseq.V1",
                stats::setNames("GC.4x44Kcustom", "Built-in: Agilent custom 4\u00d744K (GC)"),
                "Built-in: Agilent 244K"                                 = "Agilent_244K",
                stats::setNames("commercial_1x44k_postMeanCollapse_WashU", "Built-in: 1\u00d744K (WashU, post-collapse)"),
                stats::setNames("commercial_4x44k_postMeanCollapse_WashU_v2", "Built-in: 4\u00d744K (WashU v2, post-collapse)"),
                "Built-in: htp1.5 (WU update)"                           = "htp1.5_WU_update",
                "Built-in: arrayTrain (post-collapse)"                   = "arrayTrain_postMeanCollapse"
              ),
              selected = if (!allow_custom && identical(defaults$external, "Given.mdns"))
                  "nCounter" else defaults$external,
              width = "100%",
              options = list(
                placeholder = "Choose a reference medians set\u2026",
                openOnFocus = TRUE,
                dropdownParent = "body"
              )
            ),
            if (allow_custom) shiny::conditionalPanel(
              ns = ns,
              condition = "input.external == 'Given.mdns'",
              shiny::tagList(
                shiny::fileInput(
                  ns("medians"),
                  "Upload custom medians (.csv / .txt)",
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv", ".txt")
                ),
                shiny::helpText(shiny::HTML(
                  "Expected format: two columns \u2014 ",
                  "<code>X</code> (gene symbol) and <code>Given.mdns</code> (median)."
                ))
              )
            )
          )
        ),

        shiny::conditionalPanel(
          ns = ns,
          condition = "input.calibration == 'Internal'",
          shiny::selectizeInput(
            ns("internal"), "Internal calibration method",
            choices = list(
              "medianCtr (parker.original)" = "medianCtr",
              "meanCtr (genefu.scale)"      = "meanCtr",
              "qCtr (genefu.robust)"        = "qCtr"
            ),
            selected = defaults$internal, width = "100%",
            options = list(openOnFocus = TRUE, dropdownParent = "body")
          )
        )
      )
    )
  )
}
