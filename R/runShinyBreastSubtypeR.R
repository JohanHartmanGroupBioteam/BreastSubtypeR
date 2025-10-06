#' Launch the iBreastSubtypeR Shiny app
#'
#' @description Starts the Shiny UI bundled with the BreastSubtypeR package.
#' The launcher can (optionally) attach Shiny/Bslib so UI/server can use
#' unqualified functions like `tags`, `icon`, `fileInput`, etc.
#'
#' @param attach Character vector of packages to attach before launch.
#'        Defaults to c("shiny","bslib"). Set to character(0) to skip attaching.
#' @param attach_tidyverse Logical; if TRUE and tidyverse is installed, it will
#'        be attached quietly for the session (purely optional).
#' @param max_upload_mb Numeric; Shiny upload size limit (in MB). Default 1000.
#' @return Opens the app; returns the value of `shiny::runApp()`.
#'
#' @examples
#' if (interactive()) {
#'   iBreastSubtypeR()
#'   iBreastSubtypeR(attach = character(0))
#' }
#' 
#' @export
iBreastSubtypeR <- function(attach = c("shiny", "bslib"),
    attach_tidyverse = FALSE,
    max_upload_mb = 1000) {
    # helper: check install + attach quietly
    .attach_if <- function(pkgs) {
        pkgs <- unique(pkgs)
        if (!length(pkgs)) {
            return(invisible(NULL))
        }
        installed <- rownames(utils::installed.packages())
        missing <- setdiff(pkgs, installed)
        if (length(missing)) {
            stop(
                sprintf(
                    "Please install required package(s) before launching the app: %s",
                    paste(missing, collapse = ", ")
                ),
                call. = FALSE
            )
        }
        # Attach to search path so ui.R/server.R can use unqualified calls
        for (p in pkgs) {
            suppressPackageStartupMessages(
                requireNamespace(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
            )
        }
        invisible(NULL)
    }

    # Attach Shiny/Bslib (and optionally tidyverse) for this R session
    .attach_if(attach)
    if (isTRUE(attach_tidyverse) && "tidyverse" %in% rownames(utils::installed.packages())) {
        suppressPackageStartupMessages(
            requireNamespace("tidyverse", quietly = TRUE, warn.conflicts = FALSE, character.only = TRUE)
        )
    }

    # App directory shipped inside the package
    appDir <- system.file("ShinyBreastSubtypeR", package = "BreastSubtypeR")
    if (identical(appDir, "") || !dir.exists(appDir)) {
        stop("Shiny app directory 'inst/ShinyBreastSubtypeR' not found in the package.", call. = FALSE)
    }
    
    # Increase upload limit (MB -> bytes)
    options(shiny.maxRequestSize = max_upload_mb * 1024^2)

    shiny::runApp(appDir, display.mode = "normal")
}

#' @title (Deprecated) Run iBreastSubtypeR
#' @description Back-compat wrapper; use [iBreastSubtypeR()] instead.
#' @param ... Arguments passed on to [iBreastSubtypeR()].
#' @export
runShinyBreastSubtypeR <- function(...) {
    .Deprecated("iBreastSubtypeR",
        package = "BreastSubtypeR",
        msg = "runShinyBreastSubtypeR() is deprecated; use iBreastSubtypeR()."
    )
    iBreastSubtypeR(...)
}
