# --- helpers ---
safeCssUnit <- function(x, fallback = "auto") {
  if (is.null(x) || (length(x) == 0)) return(fallback)
  tryCatch(shiny::validateCssUnit(x),
           error = function(e) if (is.numeric(x) && is.finite(x)) paste0(x, "px") else fallback)
}

bs_callout <- function(title, body, color = c("info","success","warning","danger","secondary")) {
  color <- match.arg(color)
  if (requireNamespace("bslib", quietly = TRUE)) {
    ns <- asNamespace("bslib")
    if (exists("callout", envir = ns, inherits = FALSE)) {
      return(get("callout", envir = ns)(title = title, body, color = color))
    }
  }
  cl <- switch(color, "success"="alert alert-success","warning"="alert alert-warning",
               "danger"="alert alert-danger","secondary"="alert alert-secondary","info"="alert alert-info")
  shiny::tags$div(class = cl, role = "alert", shiny::tags$div(shiny::tags$b(title)), body)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

.upload_separator <- function(name) {
  if (tolower(tools::file_ext(name)) == "csv") "," else "\t"
}

# Preserve matching keys and TSIZE labels needed for clinical scoring.
.read_table <- function(inFile, key) {
  read_metadata <- getFromNamespace(".read_shiny_metadata", "BreastSubtypeR")
  read_metadata(inFile$datapath, .upload_separator(inFile$name), key)
}

# Small helper to print tables like "ER+=42, ER-=18"
.tab_str <- function(tab) {
  if (is.null(tab)) return("n/a")
  paste(sprintf("%s=%s", names(tab), as.integer(tab)), collapse = " · ")
}

# Helper: summarize cohort depending on what columns exist. The cohort-aware
# modes call get_methods(), which requires both ER and HER2 columns (TN is
# optional), so readiness is FALSE whenever one of them is absent even if the
# present columns are validly coded.
.summarize_cohort <- function(ph) {
  out <- .summarize_cohort_columns(ph)
  missing_cols <- setdiff(c("ER", "HER2"), names(ph))
  if (isTRUE(out$ok) && length(missing_cols)) {
    out$ok  <- FALSE
    out$msg <- sprintf("Cohort-aware analysis requires both 'ER' and 'HER2' columns; missing: %s",
                       paste(missing_cols, collapse = ", "))
  }
  out
}

.summarize_cohort_columns <- function(ph) {
  out <- list(kind = "none", ok = FALSE, msg = "No cohort columns found", stats = NULL)
  
  tn_present <- FALSE
  if ("TN" %in% names(ph)) {
    vals <- na.omit(as.character(ph$TN))
    bad  <- setdiff(unique(vals), c("TN","nonTN"))
    if (length(bad)) {
      out$kind <- "TN"; out$msg <- paste("Invalid TN values:", paste(bad, collapse=", "))
      return(out)
    }
    tn_present <- TRUE
    out$kind  <- "TN"
    out$ok    <- TRUE
    out$msg   <- "OK"
    out$stats <- list(
      nTN    = sum(ph$TN == "TN",    na.rm = TRUE),
      nNonTN = sum(ph$TN == "nonTN", na.rm = TRUE)
    )
  }
  
  haveER   <- "ER"   %in% names(ph)
  haveHER2 <- "HER2" %in% names(ph)
  if (haveER || haveHER2) {
    msg <- character()
    if (haveER) {
      er <- na.omit(as.character(ph$ER))
      badER <- setdiff(unique(er), c("ER+","ER-"))
      if (length(badER)) msg <- c(msg, paste("Invalid ER:", paste(badER, collapse=", ")))
    }
    if (haveHER2) {
      h2 <- na.omit(as.character(ph$HER2))
      badH <- setdiff(unique(h2), c("HER2+","HER2-"))
      if (length(badH)) msg <- c(msg, paste("Invalid HER2:", paste(badH, collapse=", ")))
    }
    
    stats_here <- list()
    if (haveER)   { stats_here$nERpos <- sum(ph$ER=="ER+",   na.rm=TRUE)
    stats_here$nERneg <- sum(ph$ER=="ER-",   na.rm=TRUE) }
    if (haveHER2) { stats_here$nHpos  <- sum(ph$HER2=="HER2+",na.rm=TRUE)
    stats_here$nHneg  <- sum(ph$HER2=="HER2-",na.rm=TRUE) }
    if (haveER && haveHER2) {
      stats_here$n_ERpos_HER2neg <- sum(ph$ER=="ER+" & ph$HER2=="HER2-", na.rm=TRUE)
      stats_here$n_ERpos_HER2pos <- sum(ph$ER=="ER+" & ph$HER2=="HER2+", na.rm=TRUE)
      stats_here$n_ERneg_HER2neg <- sum(ph$ER=="ER-" & ph$HER2=="HER2-", na.rm=TRUE)
      stats_here$n_ERneg_HER2pos <- sum(ph$ER=="ER-" & ph$HER2=="HER2+", na.rm=TRUE)
    }
    
    out$stats <- c(out$stats, stats_here)
    # keep kind = "TN" if TN present; otherwise ER/HER2 kind
    if (!tn_present) {
      out$kind <- if (haveER && haveHER2) "ERHER2" else if (haveER) "ER" else "HER2"
      out$ok   <- length(msg) == 0
      out$msg  <- if (out$ok) "OK" else paste(msg, collapse="; ")
    }
    return(out)
  }
  out
}

.read_expression <- getFromNamespace(".read_shiny_expression", "BreastSubtypeR")

# tiny DOI helper
.doi <- function(doi) sprintf("<a href='https://doi.org/%s' target='_blank' rel='noopener noreferrer'>doi</a>", doi)

# Build Calls-only and Full tables for single-method outputs
.build_single_exports <- function(res, want_4) {
  stopifnot(!is.null(res$BS.all), "PatientID" %in% names(res$BS.all))
  calls_col <- if (want_4 && "BS.Subtype" %in% names(res$BS.all)) "BS.Subtype"
  else if ("BS" %in% names(res$BS.all)) "BS"
  else setdiff(names(res$BS.all), "PatientID")[1]
  nm <- if (want_4) "Call_4class" else "Call_5class"
  
  calls_tbl <- data.frame(PatientID = res$BS.all$PatientID,
                          stringsAsFactors = FALSE, check.names = FALSE)
  calls_tbl[[nm]] <- res$BS.all[[calls_col]]
  
  full_tbl <- as.data.frame(res$BS.all, stringsAsFactors = FALSE, check.names = FALSE)
  if ("BS" %in% names(full_tbl))        { full_tbl$Call_5class <- full_tbl$BS;         full_tbl$BS <- NULL }
  if ("BS.Subtype" %in% names(full_tbl)){ full_tbl$Call_4class <- full_tbl$BS.Subtype; full_tbl$BS.Subtype <- NULL }
  
  if (!is.null(res$score.ROR) && nrow(res$score.ROR) > 0) {
    sr <- as.data.frame(res$score.ROR, stringsAsFactors = FALSE, check.names = FALSE)
    if (!"PatientID" %in% names(sr)) {
      if ("sample" %in% names(sr)) names(sr)[names(sr) == "sample"] <- "PatientID" else sr$PatientID <- rownames(sr)
    }
    keep <- setdiff(names(sr), names(full_tbl))
    if (length(keep)) {
      full_tbl <- merge(full_tbl, sr[, c("PatientID", keep), drop = FALSE],
                        by = "PatientID", all.x = TRUE, sort = FALSE)
    }
  }
  list(calls = calls_tbl, full = full_tbl)
}

# Posterior rows and rule columns in the native result have positional sample
# axes. Verify them against the actual model input before aligning to calls.
.ssp_export_metrics <- function(result, model, input_ids, call_ids) {
  is_ror <- identical(model, "ssp.ror")
  k <- if (identical(model, "ssp.pam50")) "24" else "21"
  labels <- if (is_ror) sprintf("c%03d", seq(5, 90, 5)) else c("Basal", "Her2", "LumA", "LumB")
  if (identical(model, "ssp.pam50")) labels <- c(labels, "Normal")

  ids <- rownames(result$cl)
  post <- result$all.probs[[k]]
  rules <- result$rules.matrix
  if (!length(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids) ||
      !identical(ids, input_ids) || length(call_ids) != length(ids) ||
      anyDuplicated(call_ids) || !setequal(call_ids, ids) ||
      !identical(colnames(result$data.used), ids) ||
      !identical(colnames(result$cl), k) || !is.matrix(result$prob) ||
      !identical(dimnames(result$prob), dimnames(result$cl)) ||
      !is.matrix(post) || nrow(post) != length(ids) || !identical(colnames(post), labels) ||
      (!is.null(rownames(post)) && !identical(rownames(post), ids)) ||
      !is.matrix(rules) || !nrow(rules) || ncol(rules) != length(ids) ||
      (!is.null(colnames(rules)) && !identical(colnames(rules), ids))) {
    stop(sprintf("Cannot align %s metrics. Check sample IDs and native output dimensions.", model),
         call. = FALSE)
  }

  prefix <- if (is_ror) "SSP_ROR_" else "SSP_"
  metrics <- data.frame(row.names = ids)
  if (is_ror) metrics$SSP_ROR_Bin <- result$cl[, k]
  metrics[[if (is_ror) "SSP_ROR_Bin_Prob" else "SSP_Call_Prob"]] <- result$prob[, k]
  metrics[paste0(prefix, "Prob_", labels)] <- as.data.frame(post)
  metrics[[paste0(prefix, "Rules_Evaluable")]] <- as.integer(colSums(!is.na(rules)))
  metrics[[paste0(prefix, "Rules_Total")]] <- nrow(rules)
  metrics[match(call_ids, ids), , drop = FALSE]
}

.ssp_ror_coverage_note <- function(full) {
  incomplete <- full$SSP_ROR_Rules_Evaluable < full$SSP_ROR_Rules_Total
  if (!any(incomplete)) return(character())
  note <- sprintf(paste("Incomplete SSP.ROR rule coverage in %d of %d samples.",
                        "Inspect SSP_ROR_Rules_Evaluable out of %d rules per sample."),
                  sum(incomplete), length(incomplete), full$SSP_ROR_Rules_Total[1L])
  zero <- sum(full$SSP_ROR_Rules_Evaluable == 0L)
  if (zero) {
    note <- paste(note, sprintf("Rows with SSP_ROR_Rules_Evaluable = 0 (%d) have zero evaluable rules.", zero),
                  "Their native bins and posteriors have no expression support.",
                  "Check probe-to-Entrez annotation and missing expression for these samples.")
  }
  note
}

# Human-readable labels travel in comment metadata, never as renamed columns.
.experimental_ror_pc_export_note <- function(df) {
  has_score <- getFromNamespace(".has_experimental_ror_pc", "BreastSubtypeR")
  if (!has_score(df)) return(character())
  fields <- getFromNamespace(".experimental_ror_pc_fields", "BreastSubtypeR")()
  notice <- getFromNamespace(".experimental_ror_pc_notice", "BreastSubtypeR")()
  paste(notice, "Output labels:",
        paste(sprintf("%s = %s", names(fields), fields), collapse = "; "),
        "The original Parker ROR-PC fields retain their existing labels.")
}

server <- function(input, output, session) {
  
  # --- small reactives ---
  has_clinical <- reactive({ isTRUE((input$hasClinical %||% FALSE)) })
  k_is_4class  <- reactive({ identical(input$k_subtypes, "4") })
  
  # --- citation header for exports ---
  PKG_VER      <- as.character(utils::packageVersion("BreastSubtypeR"))
  CITATION_DOI <- "10.1093/nargab/lqaf131"
  CITATION_NOTE <- sprintf("BreastSubtypeR v%s - NAR Genomics and Bioinformatics (2025) DOI: %s",
                           PKG_VER, CITATION_DOI)
  
  .write_tsv_with_header <- function(path, df, header_lines = CITATION_NOTE, quote_ids = FALSE) {
    con <- file(path, open = "wt", encoding = "UTF-8")
    on.exit(close(con), add = TRUE)
    if (length(header_lines)) writeLines(paste0("# ", header_lines), con)
    utils::write.table(df, con, sep = "\t", quote = if (quote_ids) 1L else FALSE,
                       qmethod = "double", row.names = FALSE)
  }
  
  .pam50_variant_from_internal <- function(internal) {
    switch(internal,
           "medianCtr" = "parker.original",
           "meanCtr"   = "genefu.scale",
           "qCtr"      = "genefu.robust",
           "parker.original")
  }
  
  .method_header_line <- function(params, want_4, use_clin) {
    method_label <- params$method %||% "Unknown"
    if (identical(method_label, "PAM50")) {
      cal <- params$calibration %||% "None"
      if (identical(cal, "Internal")) {
        internal <- params$internal %||% "medianCtr"
        variant  <- .pam50_variant_from_internal(internal)
        method_label <- sprintf("PAM50 (variant: %s, calibration: Internal (%s))", variant, internal)
      } else if (identical(cal, "External")) {
        ext <- params$external %||% "Unknown"
        ext_label <- if (identical(ext, "Given.mdns")) "Custom (Given.mdns)" else ext
        method_label <- sprintf("PAM50 (variant: parker.original, calibration: External (%s))", ext_label)
      } else {
        method_label <- "PAM50 (variant: parker.original, calibration: None)"
      }
    }
    sprintf("Method: %s | Classes: %s | ROR: %s",
            method_label, if (want_4) "4" else "5", if (use_clin) "ON" else "OFF")
  }
  
  # Raw inputs, mapping and predictions belong to this session only.
  reactive_files <- shiny::reactiveValues(
    GEX = NULL, clinic = NULL, anno = NULL, source = NULL,
    RawCounts = NULL, example_id = NULL, input_generation = 0L,
    data_input = NULL, mapped_key = NULL, prediction = NULL, error = NULL
  )

  examples <- getFromNamespace(".shiny_examples", "BreastSubtypeR")
  read_example <- getFromNamespace(".read_shiny_example", "BreastSubtypeR")
  manual_mode <- reactiveVal("norm")
  raw_mode <- reactive({
    id <- reactive_files$example_id
    if (is.null(id)) identical(manual_mode(), "raw") else examples[[id]]$raw
  })

  # New file bindings on replacement also discard uploads still in flight.
  upload_generation <- reactiveVal(0L)
  upload_id <- function(id) {
    generation <- upload_generation()
    if (generation == 0L) id else paste(id, generation, sep = "_")
  }
  upload_labels <- c(GEX = "Upload expression matrix", clinic = "Upload clinical table",
                     anno = "Upload annotation table")
  for (key in names(upload_labels)) local({
    id <- key
    output[[paste0("upload_", id)]] <- renderUI({
      fileInput(upload_id(id), upload_labels[[id]],
                accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv", ".txt"))
    })
  })
  # Retire the previous file immediately when its replacement starts uploading.
  # Shiny supplies a new temporary path only after the transfer succeeds.
  retired_files <- reactiveVal(list())
  uploads <- reactive({
    files <- stats::setNames(vector("list", length(upload_labels)), names(upload_labels))
    retired <- retired_files()
    for (id in names(files)) {
      file <- input[[upload_id(id)]]
      if (!identical(file, retired[[id]])) files[id] <- list(file)
    }
    files
  })
  input_key <- reactive({
    list(generation = reactive_files$input_generation, source = reactive_files$source,
         example = reactive_files$example_id, raw = raw_mode(),
         uploads = if (identical(reactive_files$source, "upload")) uploads())
  })

  # Rebind calibration on each primary input change. Primary upload controls
  # retain their own lifetime so the other two slots survive a manual edit.
  default_calibration <- getFromNamespace(".shiny_pam50_defaults", "BreastSubtypeR")
  pam50_controls <- getFromNamespace(".shiny_pam50_controls", "BreastSubtypeR")
  calibration_defaults <- reactiveVal(default_calibration)
  calibration_ns <- reactive({
    generation <- reactive_files$input_generation
    NS(if (generation > 0L) paste0("data", generation) else NULL)
  })
  output$pam50_calibration <- renderUI({
    pam50_controls(calibration_ns(), calibration_defaults())
  })
  calibration_choice <- function(id) {
    input[[calibration_ns()(id)]] %||% calibration_defaults()[[id]]
  }
  # Shiny assigns each completed file upload a fresh temporary path. Remember
  # its owning generation even if the custom panel is subsequently hidden.
  reference_owners <- new.env(parent = emptyenv())
  retiring_reference <- reactiveVal(NULL)
  reference_file <- reactive({
    file <- input[[calibration_ns()("medians")]]
    if (!is.data.frame(file) || nrow(file) != 1L ||
        !is.character(file$datapath) || length(file$datapath) != 1L ||
        is.na(file$datapath) || !nzchar(file$datapath)) return(NULL)
    path <- normalizePath(file$datapath, winslash = "/", mustWork = FALSE)
    if (identical(path, retiring_reference())) return(NULL)
    generation <- reactive_files$input_generation
    owner <- reference_owners[[path]]
    if (!is.null(owner) && !identical(owner, generation)) return(NULL)
    reference_owners[[path]] <- generation
    file
  })
  observeEvent(input[[calibration_ns()("medians")]], reference_file(), priority = 115)
  advance_dataset <- function() {
    reference_file() # claim any outgoing file even if its panel is hidden
    choices <- default_calibration
    for (id in names(choices)) choices[[id]] <- calibration_choice(id)
    if (identical(choices$calibration, "External") &&
        identical(choices$external, "Given.mdns")) choices <- default_calibration
    calibration_defaults(choices)
    reactive_files$input_generation <- reactive_files$input_generation + 1L
    retiring_reference(NULL)
  }

  prediction_parameters <- reactive({
    method <- input$BSmethod
    params <- list(method = method,
                   classes = if (identical(method, "AIMS")) "5" else input$k_subtypes)
    if (!is.null(method) && method %in% c("PAM50", "cIHC", "cIHC.itr", "PCAPAM50", "ssBC"))
      params$clinical <- isTRUE(input$hasClinical)
    if (identical(method, "PAM50")) {
      params$calibration <- calibration_choice("calibration")
      if (identical(params$calibration, "Internal")) params$internal <- calibration_choice("internal")
      if (identical(params$calibration, "External")) {
        params$external <- calibration_choice("external")
        if (identical(params$external, "Given.mdns")) params$medians <- reference_file()
      }
    }
    if (identical(method, "cIHC.itr")) {
      params$ratio <- input$ratio
      params$iteration <- input$iteration
    }
    if (identical(method, "ssBC")) params$subgroup <- input$s
    params
  })

  invalidate_prediction <- function() {
    reactive_files$prediction <- NULL
    reactive_files$error <- NULL
    removeNotification("analysis-error")
  }
  invalidate_mapping <- function() {
    reactive_files$data_input <- NULL
    reactive_files$mapped_key <- NULL
    invalidate_prediction()
  }
  report_error <- function(e) {
    invalidate_prediction()
    # req() and need(..., FALSE) abort silently; named validation stays visible.
    if (inherits(e, "shiny.silent.error") && !nzchar(conditionMessage(e)))
      return(invisible(NULL))
    reactive_files$error <- conditionMessage(e)
    showNotification(reactive_files$error, type = "error", duration = NULL,
                     id = "analysis-error")
  }
  output$analysis_error <- renderUI({
    if (is.null(reactive_files$error)) return(NULL)
    bs_callout("Analysis needs attention", reactive_files$error, "danger")
  })

  clear_inputs <- function() {
    invalidate_mapping()
    reactive_files$example_id <- reactive_files$source <- NULL
    reactive_files$GEX <- reactive_files$clinic <- reactive_files$anno <- NULL
    reactive_files$RawCounts <- NULL
  }
  sync_mode <- function() {
    session$sendCustomMessage("bsrDataType", list(
      locked = !is.null(reactive_files$example_id),
      value = if (raw_mode()) "raw" else "norm"))
  }
  observeEvent(input$is_raw_counts, {
    if (!is.null(reactive_files$example_id)) {
      sync_mode()
    } else if (length(input$is_raw_counts) == 1L && input$is_raw_counts %in% c("raw", "norm") &&
               !identical(input$is_raw_counts, manual_mode())) {
      if (identical(reactive_files$source, "upload")) advance_dataset()
      manual_mode(input$is_raw_counts)
      reactive_files$RawCounts <- NULL
      invalidate_mapping()
    }
  }, priority = 110)

  activate_upload <- function() {
    advance_dataset()
    from_example <- !is.null(reactive_files$example_id)
    clear_inputs()
    reactive_files$source <- "upload"
    if (from_example) sync_mode()
  }
  observeEvent(input$upload_started, {
    if (identical(input$upload_started, calibration_ns()("medians"))) {
      file <- reference_file()
      if (!is.null(file)) {
        path <- normalizePath(file$datapath, winslash = "/", mustWork = FALSE)
        reference_owners[[path]] <- -1L
        retiring_reference(path)
      }
      return()
    }
    ids <- vapply(names(upload_labels), upload_id, character(1))
    slot <- names(ids)[ids == input$upload_started]
    if (length(slot) != 1L) return()
    retired <- retired_files()
    retired[slot] <- list(input[[ids[[slot]]]])
    retired_files(retired)
    activate_upload()
  }, priority = 105)
  observeEvent(uploads(), {
    if (any(lengths(uploads()) > 0L) || identical(reactive_files$source, "upload"))
      activate_upload()
  }, priority = 100)

  observeEvent(prediction_parameters(), {
    # Rebinding the calibrated controls can echo their existing defaults.
    p <- reactive_files$prediction
    if (mapping_current() && (is.null(p) || !identical(p$parameters, prediction_parameters())))
      invalidate_prediction()
  }, ignoreInit = TRUE, priority = 100)

  mapping_current <- reactive({
    !is.null(reactive_files$data_input) &&
      identical(reactive_files$mapped_key, input_key())
  })
  current_input_coverage <- reactive({
    if (!mapping_current()) return(NULL)
    tryCatch(
      list(data = BreastSubtypeR::InputCoverage(reactive_files$data_input),
           error = NULL),
      error = function(e) list(data = NULL, error = conditionMessage(e))
    )
  })
  output$input_coverage <- renderUI({
    coverage <- current_input_coverage()
    if (is.null(coverage)) return(NULL)
    if (!is.null(coverage$error)) {
      return(bslib::card(
        bslib::card_header("Input coverage"),
        tags$p(paste("Input coverage is unavailable:", coverage$error))
      ))
    }
    tab <- coverage$data
    rules <- vapply(seq_len(nrow(tab)), function(i) {
      if (is.na(tab$rules_total[i])) return("N/A")
      if (tab$rules_evaluable_min[i] == tab$rules_evaluable_max[i]) {
        sprintf("%s/%s", tab$rules_evaluable_min[i], tab$rules_total[i])
      } else {
        sprintf("%s-%s/%s per sample", tab$rules_evaluable_min[i],
                tab$rules_evaluable_max[i], tab$rules_total[i])
      }
    }, character(1))
    bslib::card(
      bslib::card_header("Input coverage"),
      tags$div(class = "table-responsive", tabindex = "0", role = "region",
        `aria-label` = "Input coverage by method and model",
        tags$table(class = "table table-sm",
          tags$thead(tags$tr(
            tags$th(scope = "col", "Method/model"),
            tags$th(scope = "col", "Genes represented"),
            tags$th(scope = "col", "Genes absent"),
            tags$th(scope = "col", "Rules evaluable")
          )),
          tags$tbody(lapply(seq_len(nrow(tab)), function(i) {
            tags$tr(tags$th(scope = "row", tab$model[i]),
                    tags$td(sprintf("%s/%s", tab$genes_present[i], tab$genes_total[i])),
                    tags$td(tab$genes_absent[i]), tags$td(rules[i]))
          }))
        )
      ),
      tags$p("NC-based methods share the 50-gene PAM50 subtype reference; AIMS and both SSPBC models are SSP-based. These rows describe input coverage, not which methods were selected or run."),
      tags$p("Descriptive coverage of processed expression inputs. Coverage does not establish prediction reliability or change method eligibility.")
    )
  })
  current_prediction <- reactive({
    p <- reactive_files$prediction
    if (is.null(p) || !mapping_current() ||
        !identical(p$inputs, input_key()) ||
        !identical(p$parameters, prediction_parameters())) return(NULL)
    p
  })

  # Use all aligned probes in both mappings, including for raw library sizes.
  aligned_input <- function() {
    samples <- intersect(colnames(reactive_files$GEX), reactive_files$clinic$PatientID)
    validate(need(length(samples) > 0, "No overlapping sample IDs between GEX columns and clinical ‘PatientID’."))
    probeID <- intersect(rownames(reactive_files$GEX), rownames(reactive_files$anno))
    validate(need(length(probeID) > 0, "No overlapping probes between GEX and annotation."))

    if (isTRUE(reactive_files$RawCounts)) {
      validate(need("Length" %in% names(as.data.frame(reactive_files$anno)),
                    "Raw counts selected but annotation lacks ‘Length’ (bp)."))
    }

    SummarizedExperiment::SummarizedExperiment(
      assays  = list(counts = reactive_files$GEX[probeID, samples, drop = FALSE]),
      rowData = S4Vectors::DataFrame(reactive_files$anno[probeID, , drop = FALSE]),
      colData = S4Vectors::DataFrame(reactive_files$clinic[samples, , drop = FALSE])
    )
  }

  # A previous run's full radio value can still be in flight when renderUI
  # replaces the control. Only the current SSP run can opt in to its mapping.
  export_kind <- reactive({
    choice <- input$export_kind %||% "calls"
    if (identical(choice, "calls")) return("calls")
    full_choice <- current_prediction()$ssp$full_choice %||% "full"
    if (identical(choice, full_choice)) "full" else NULL
  })

  # This reactive caches success or failure for the current prediction. Keep
  # export_kind outside it so switching back to calls does not evict the cache.
  ssp_full <- reactive({
    p <- current_prediction()
    req(p$ssp)
    tryCatch({
      subtype <- .ssp_export_metrics(p$ssp$native, p$ssp$model, p$ssp$input_ids, p$calls$PatientID)
      mapped <- BreastSubtypeR::Mapping(se_obj = aligned_input(), RawCounts = reactive_files$RawCounts,
                                        impute = FALSE, include_ssp_ror = TRUE, verbose = TRUE)
      ror_input <- mapped$se_SSP_ROR
      ror <- BreastSubtypeR::BS_sspbc(se_obj = ror_input, ssp.name = "ssp.ror")
      metrics <- .ssp_export_metrics(ror, "ssp.ror", colnames(ror_input), p$calls$PatientID)
      full <- cbind(p$calls, subtype, metrics)
      list(full = full, note = .ssp_ror_coverage_note(full))
    }, error = function(e) list(error = paste("SSP.ROR full export unavailable.", conditionMessage(e))))
  })

  selected_export <- reactive({
    p <- current_prediction()
    if (!is.null(p$ssp) && identical(export_kind(), "full")) {
      result <- ssp_full()
      p$full <- result$full
      p$full_error <- result$error
      p$full_note <- result$note
      p$header <- sprintf("Method: %s | Classes: %d | ROR: categorical ssp.ror | Probabilities: fitted class/bin membership",
                          p$ssp$model, if (p$want_4) 4L else 5L)
    }
    p
  })

  supports_export <- function(p, kind) {
    !is.null(p) && length(kind) == 1L && kind %in% c("calls", "full") &&
      !is.null(p[[kind]]) && nrow(p[[kind]]) > 0L
  }

  # Check validity inside the callback, including direct requests to old URLs.
  export_result <- function(kind) {
    p <- selected_export()
    if (is.null(p)) stop("No current prediction. Preprocess if needed, then run subtyping.", call. = FALSE)
    if (identical(kind, "full") && !is.null(p$full_error)) stop(p$full_error, call. = FALSE)
    if (!supports_export(p, kind))
      stop("The current prediction does not provide this export. Choose calls only.", call. = FALSE)
    p
  }
  download_filename <- function() {
    kind <- export_kind()
    p <- export_result(kind)
    sprintf("results-%s-%s-%s.txt", p$method,
            if (p$want_4) "4class" else "5class", kind)
  }
  write_download <- function(file) {
    kind <- export_kind()
    p <- export_result(kind)
    hdr <- c(CITATION_NOTE, p$header)
    if (!identical(p$method, "AUTO")) {
      hdr <- c(hdr, sprintf("Export: %s | Date: %s | R: %s", kind,
                           format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), getRversion()))
    }
    if (identical(kind, "full")) {
      note <- .experimental_ror_pc_export_note(p[[kind]])
      if (length(note)) hdr[length(hdr)] <- paste(hdr[length(hdr)], note, sep = " | ")
    }
    # Protect SSP sample IDs containing TSV delimiters, quotes or the comment marker.
    quote_ids <- !is.null(p$ssp) && any(grepl('[\t\r\n"#]', p$calls$PatientID))
    .write_tsv_with_header(file, p[[kind]], header_lines = c(hdr, p$full_note), quote_ids = quote_ids)
  }
  output$download <- downloadHandler(filename = download_filename, content = write_download)

  output$plotSection <- renderUI({
    p <- current_prediction()
    if (is.null(p)) return(NULL)
    if (identical(p$method, "AUTO")) {
      bslib::card(
        bslib::card_header("Method comparison and entropy"),
        tags$p(class = "plot-scroll-hint", "Swipe horizontally to view the full plot."),
        tags$div(class = "plot-scroll", tabindex = "0", role = "region",
          `aria-label` = "Method comparison plot; scroll horizontally on narrow screens",
          plotOutput("multi_plot", height = 560)))
    } else {
      bslib::layout_columns(
        col_widths = bslib::breakpoints(sm = c(12, 12), lg = c(5, 7)),
        bslib::card(bslib::card_header("Subtype calls"),
          tags$div(class = "plot-scroll", tabindex = "0", role = "region",
            `aria-label` = "Subtype proportion plot; scroll horizontally on narrow screens",
            plotOutput("pie1"))),
        bslib::card(bslib::card_header("Expression and subtype calls"),
          tags$p(class = "plot-scroll-hint", "Swipe horizontally to view the full plot."),
          tags$div(class = "plot-scroll", tabindex = "0", role = "region",
            `aria-label` = "Expression heatmap; scroll horizontally on narrow screens",
            plotOutput("heat2", height = 720))))
    }
  })
  output$multi_plot <- renderPlot({
    p <- current_prediction(); req(p, identical(p$method, "AUTO"))
    BreastSubtypeR::Vis_Multi(p$plot)
  }, res = 96, alt = "Subtype calls across selected methods, with per-sample entropy and subtype legend.")
  output$pie1 <- renderPlot({
    p <- current_prediction(); req(p, !identical(p$method, "AUTO"))
    BreastSubtypeR::Vis_pie(p$plot)
  }, res = 96, alt = "Proportions of subtype calls for the selected method.")
  output$heat2 <- renderPlot({
    p <- current_prediction(); req(p, !identical(p$method, "AUTO"))
    BreastSubtypeR::Vis_heatmap(as.matrix(p$matrix), out = p$plot)
  }, res = 96, alt = "Expression heatmap annotated with subtype calls for the selected method.")

  # ---- UI dynamic help ----
  output$gex_help <- renderUI({
    if (raw_mode()) {
      div(class = "method-help",
          tags$b("Expression matrix requirements"),
          tags$ul(
            tags$li(tags$b("Format:"), " genes × samples ", tags$em("integer raw counts (RNA-seq only)"), "."),
            tags$li(tags$b("Rows:"), " gene IDs matching ", tags$code("annotation$probe"), "."),
            tags$li(tags$b("Columns:"), " sample IDs matching ", tags$code("clinical$PatientID"), "."),
            tags$li(tags$b("Processing:"), " NC → log2-CPM (upper quartile); SSP → FPKM (linear) using gene Length.")
          ))
    } else {
      div(class = "method-help",
          tags$b("Expression matrix requirements"),
          tags$ul(
            tags$li(tags$b("Format:"), " genes × samples ", tags$em("log2-normalized expression"), "."),
            tags$li(tags$b("Examples:"), " ", tags$code("log2(FPKM+1)"), ", microarray/nCounter log2."),
            tags$li(tags$b("Rows:"), " gene IDs matching ", tags$code("annotation$probe"), "."),
            tags$li(tags$b("Columns:"), " sample IDs matching ", tags$code("clinical$PatientID"), "."),
            tags$li(tags$b("Processing:"), " NC → uses log2 values; SSP → back-transforms to linear.")
          ))
    }
  })
  
  output$clin_help <- renderUI({
    if (has_clinical()) {
      div(class = "method-help",
          tags$b("Clinical data requirements"),
          HTML("<ul>
              <li><b>Minimum:</b> <code>PatientID</code>.</li>
              <li><b>Method-specific:</b> cIHC/itr/PCAPAM50 need <code>ER</code>. ssBC needs <code>ER</code> (and <code>HER2</code> for ER.v2) or <code>TN</code>.</li>
              <li><b>For clinical ROR:</b> Supply both <code>TSIZE</code> and <code>NODE</code>.
              Missing clinical observations or nonfinite NODE disable clinical ROR for the cohort. In uploaded TSIZE fields, blank, NA and NaN denote missing values.
              With complete annotations, an invalid TSIZE gives a warning and missing size-dependent research ROR scores and groups only for the affected patients, including those with four or more positive nodes. Subtype calls and genomic ROR remain available.</li>
              <li><b>Clinical coding:</b> <code>TSIZE</code> is 0 (at most 2 cm) or 1 (more than 2 cm).
              <code>NODE</code> is the number of positive lymph nodes. The additional subtype-specific ROR-PC risk group requires finite nonnegative integer counts;
              negative or fractional values give a warning and a missing group for that sample, while numerical ROR scores remain available.
              Binary 0/1 is accepted but cannot distinguish 1&ndash;3 from 4 or more nodes for the additional subtype-specific ROR-PC risk group.
              NODE does not enter numerical ROR scores.</li>
            </ul>"))
    } else {
      div(class = "method-help",
          tags$b("Clinical data requirements"),
          HTML("<ul>
              <li><b>Minimum:</b> <code>PatientID</code>.</li>
              <li><b>Optional for NC methods (clinical ROR):</b> binary tumor size <code>TSIZE</code> and positive-node count <code>NODE</code>. Use actual counts when available; binary 0/1 cannot distinguish 1&ndash;3 from 4 or more nodes for the additional subtype-specific ROR-PC risk group.</li>
            </ul>"))
    }
  })
  
  output$anno_help <- renderUI({
    if (raw_mode()) {
      div(class = "method-help",
          tags$b("Feature annotation requirements"),
          HTML("<ul>
              <li><code>probe</code> — matches GEX row names.</li>
              <li><code>ENTREZID</code> — mandatory.</li>
              <li><code>Length</code> — gene length (bp), mandatory for raw counts.</li>
            </ul>"))
    } else {
      div(class = "method-help",
          tags$b("Feature annotation requirements"),
          HTML("<ul>
              <li><code>probe</code> — matches GEX row names.</li>
              <li><code>ENTREZID</code> — mandatory.</li>
            </ul>"))
    }
  })
  
  # Method help + AUTO chip + lock AIMS to 5-class unchanged (shortened)
  output$method_help <- renderUI({
    m <- input$BSmethod
    if (is.null(m) || length(m) < 1) return(NULL)
    m <- as.character(m)[1]  # ensure length-1 scalar
    
    .p <- function(title, items) {
      lis <- Map(function(k, v) tags$li(HTML(paste0("<b>", k, ":</b> ", v))), names(items), items)
      div(class = "method-help", tags$b(title), tags$br(), tags$ul(lis))
    }
    
    switch(m,
           "AIMS"   = .p("AIMS", list("Category"="SSP-based","Description"="Absolute Intrinsic Molecular Subtyping.","Key refs"="Paquet & Hallett, 2015 JNCI")),
           "sspbc"  = .p("SSPBC", list("Category"="SSP-based","Description"="SCAN-B SSP models (4 or 5 classes).","Key refs"="Staaf et al., 2022 npj Breast Cancer")),
           "PAM50"  = .p("PAM50", list("Variants"="parker.original | genefu.scale | genefu.robust")),
           NULL)
  })
  
  output$auto_chip <- renderUI({
    if (!identical(input$BSmethod, "AUTO Mode")) return(NULL)
    if (!mapping_current()) {
      return(tags$span(class = "chip auto pending",
        tags$b("AUTO: awaiting preprocessing"),
        tags$span(class = "chip-note", "Run Step 1 (Preprocess & map) to continue.")))
    }
    chk <- auto_requirements()
    if (!isTRUE(chk$show)) return(NULL)
    cls <- if (isTRUE(chk$ready)) "chip auto ready" else "chip auto blocked"
    
    stats <- ""
    if (isTRUE(chk$ready)) {
      er  <- if (!is.null(chk$ER_tab))   paste(names(chk$ER_tab), as.integer(chk$ER_tab), sep="=", collapse=" · ") else "ER n/a"
      h2  <- if (!is.null(chk$HER2_tab)) paste(names(chk$HER2_tab), as.integer(chk$HER2_tab), sep="=", collapse=" · ") else "HER2 n/a"
      stats <- sprintf(" &nbsp; <span class='chip-note'>(n=%d · %s | %s)</span>", chk$n, er, h2)
    } else if (!is.null(chk$msg)) {
      stats <- sprintf(" &nbsp; <span class='chip-note'>— %s</span>", chk$msg)
    }
    
    htmltools::tags$span(
      class = cls,
      icon("magic"),
      htmltools::HTML("<b>AUTO:</b> cohort-aware selection"),
      htmltools::HTML(stats)
    )
  })
  observeEvent(input$BSmethod, {
    if (identical(input$BSmethod, "AIMS")) {
      updateRadioButtons(session, "k_subtypes",
                         choices  = list("5 classes (includes Normal-like)" = "5"),
                         selected = "5")
    } else {
      prev <- isolate(input$k_subtypes); if (!prev %in% c("4","5")) prev <- "5"
      updateRadioButtons(session, "k_subtypes",
                         choices  = list("5 classes (includes Normal-like)" = "5",
                                         "4 classes (excludes Normal-like)" = "4"),
                         selected = prev)
    }
  }, ignoreInit = FALSE)
  
  output$export_selector <- renderUI({
    p <- current_prediction()
    if (is.null(p)) return(helpText("Run subtyping before downloading results."))
    has_full <- !is.null(p$ssp) || supports_export(p, "full")
    choices <- list("Calls only" = "calls")
    if (has_full) choices[["Full metrics (incl. ROR)"]] <- p$ssp$full_choice %||% "full"
    selected <- isolate(input$export_kind) %||% "calls"
    if (!selected %in% unlist(choices)) selected <- "calls"
    tagList(
      radioButtons("export_kind", "Export content", choices = choices,
                   inline = TRUE, selected = selected),
      if (!has_full) helpText("Full metrics are unavailable for this method. Download calls only.")
    )
  })
  output$experimental_ror_pc_notice <- renderUI({
    p <- current_prediction()
    if (!identical(input$export_kind %||% "calls", "full") ||
        !supports_export(p, "full")) return(NULL)
    has_score <- getFromNamespace(".has_experimental_ror_pc", "BreastSubtypeR")
    if (!has_score(p$full)) return(NULL)
    fields <- getFromNamespace(".experimental_ror_pc_fields", "BreastSubtypeR")()
    notice <- getFromNamespace(".experimental_ror_pc_notice", "BreastSubtypeR")()
    bs_callout("Experimental 18-gene ROR-PC", tagList(
      tags$p(notice),
      tags$dl(lapply(seq_along(fields), function(i) tagList(
        tags$dt(unname(fields[i])), tags$dd(tags$code(names(fields)[i]))))),
      tags$p("The original Parker ROR-PC fields retain their existing labels.")
    ), "warning")
  })
  output$ssp_ror_notice <- renderUI({
    p <- selected_export()
    if (!is.null(p$full_error))
      return(bs_callout("SSP.ROR export needs attention",
                        paste(p$full_error, "Choose Calls only to download the subtype results."), "danger"))
    if (length(p$full_note)) bs_callout("SSP.ROR rule coverage", p$full_note, "warning")
  })
  output$download_control <- renderUI({
    p <- selected_export()
    if (is.null(p)) return(NULL)
    kind <- export_kind()
    if (!supports_export(p, kind))
      return(NULL)
    downloadButton(
      "download", HTML("Download results<br>(TSV)"),
      class = "btn-primary",
      style = "max-width: 100%; width: 220px; padding: 10px 16px; white-space: normal; line-height: 1.4;"
    )
  })
  
  # --- AUTO preflight panel ---
  auto_requirements <- reactive({
    if (!identical(input$BSmethod, "AUTO Mode"))
      return(list(show = FALSE))
    
    di <- reactive_files$data_input
    if (is.null(di) || is.null(di$se_NC))
      return(list(show = TRUE, ready = FALSE, msg = "Run Step 1 (Preprocess & map) first."))
    
    ph <- as.data.frame(SummarizedExperiment::colData(di$se_NC))
    summ <- .summarize_cohort(ph)
    
    # Ready only when coding is valid for the detected cohort-kind
    list(
      show  = TRUE,
      ready = isTRUE(summ$ok),
      msg   = if (summ$ok) "OK" else summ$msg,
      kind  = summ$kind,
      stats = summ$stats
    )
  })
  
  output$auto_preflight <- renderUI({
    if (!mapping_current()) return(NULL)
    chk <- auto_requirements()
    if (!isTRUE(chk$show)) return(NULL)
    
    # Compact lines
    er_line <- her2_line <- subgroup_line <- tn_line <- NULL
    s <- chk$stats %||% list()
    
    if (!is.null(s$nERpos) || !is.null(s$nERneg))
      er_line <- sprintf("ER+: %s · ER-: %s", s$nERpos %||% 0, s$nERneg %||% 0)
    if (!is.null(s$nHpos)  || !is.null(s$nHneg))
      her2_line <- sprintf("HER2+: %s · HER2-: %s", s$nHpos %||% 0, s$nHneg %||% 0)
    if (!is.null(s$n_ERpos_HER2neg)) {
      subgroup_line <- paste(
        sprintf("ER+/HER2-: %s", s$n_ERpos_HER2neg %||% 0),
        sprintf("ER-/HER2-: %s", s$n_ERneg_HER2neg %||% 0),
        sprintf("ER+/HER2+: %s", s$n_ERpos_HER2pos %||% 0),
        sprintf("ER-/HER2+: %s", s$n_ERneg_HER2pos %||% 0),
        sep = " · "
      )
    }
    if (!is.null(s$nTN))
      tn_line <- sprintf("TN: %s · nonTN: %s", s$nTN %||% 0, s$nNonTN %||% 0)
    
    body <- tagList(
      HTML(sprintf("Cohort fields detected: <b>%s</b>", chk$kind)),
      if (!is.null(er_line))    tagList(tags$br(), HTML(er_line))    else NULL,
      if (!is.null(her2_line))  tagList(tags$br(), HTML(her2_line))  else NULL,
      if (!is.null(subgroup_line)) tagList(tags$br(), HTML(subgroup_line)) else NULL,
      if (!is.null(tn_line))    tagList(tags$br(), HTML(tn_line))    else NULL,
      if (identical(chk$kind, "TN")) " (ER/HER2 present; TN takes precedence)" else "",
      tags$br(),
      tags$small(if (isTRUE(chk$ready)) "Coding OK. You can run AUTO." else sprintf("Issue: %s", chk$msg))
    )
    
    if (isTRUE(chk$ready)) {
      bs_callout("AUTO preflight: ready", body, "success")
    } else {
      bs_callout("AUTO preflight: action needed", body, "danger")
    }
  })
  
  # ---- About modal ----
  observeEvent(input$about, {
    showModal(modalDialog(
      title = "About & Citation", easyClose = TRUE, footer = NULL,
      HTML(paste0(
        "<p><b>BreastSubtypeR</b> — intrinsic molecular subtyping for breast cancer (R/Bioconductor).</p>",
        "<p><b>Please cite:</b><br/>",
        "Yang Q., Hartman J., Sifakis E.G. <em>BreastSubtypeR</em>. ",
        "<em>NAR Genomics and Bioinformatics</em> (2025). ",
        "<a href='https://doi.org/10.1093/nargab/lqaf131' target='_blank' rel='noopener noreferrer'>https://doi.org/10.1093/nargab/lqaf131</a></p>",
        "<p><a href='https://academic.oup.com/nargab/pages/editors-choice-articles?login=false#:~:text=BreastSubtypeR' target='_blank' rel='noopener noreferrer'>Editor's Choice (2025)</a> - ",
        "<a href='https://doi.org/10.1093/nargab/lqaf131' target='_blank' rel='noopener noreferrer'>NAR Genomics and Bioinformatics</a>.</p>"
      ))
    ))
  })
  
  mode_label <- function(raw) if (raw) "Raw counts" else "Normalized (log2)"
  dimensions_label <- function(spec) {
    sprintf("%s %s and %d tumors", format(spec$dimensions[1L], big.mark = ","),
            if (spec$raw) "expression features" else "genes", spec$dimensions[2L])
  }
  output$active_source <- renderText({
    id <- reactive_files$example_id
    if (!is.null(id)) {
      spec <- examples[[id]]
      return(paste(spec$name, dimensions_label(spec), mode_label(spec$raw), sep = " | "))
    }
    if (identical(reactive_files$source, "upload"))
      return(sprintf("Manual upload | %s | %d of 3 files supplied",
                     mode_label(raw_mode()), sum(lengths(uploads()) > 0L)))
    "No active dataset"
  })
  open_examples <- function() {
    choices <- lapply(examples, function(spec) {
      tags$span(tags$strong(spec$name), tags$br(),
        dimensions_label(spec), tags$br(),
        if (spec$raw) "RNA-seq coverage sums" else "Log2-normalized RNA-seq expression",
        tags$br(), spec$units, tags$br(),
        paste("Data type:", mode_label(spec$raw)), tags$br(),
        if (spec$raw) "Gene length annotation: Included" else "Gene length annotation: Not required")
    })
    showModal(modalDialog(
      title = "Load example data",
      p("Select a packaged breast cancer expression dataset to explore the BreastSubtypeR workflow."),
      radioButtons("example_choice", "Dataset", choiceNames = unname(choices),
                   choiceValues = names(examples), width = "100%", selected = reactive_files$example_id %||% "oslo"),
      p("TCGA-BRCA uses the raw-count workflow with base-pair coverage sums. Conventional sequencing read counts measure reads rather than covered bases."),
      footer = tagList(modalButton("Cancel"), actionButton("confirm_example", "Load selected dataset",
                                                          class = "btn-primary"))))
  }
  observeEvent(input$load_example, open_examples())
  observeEvent(input$load_example_step1, open_examples())
  observeEvent(input$confirm_example, {
    id <- input$example_choice
    removeModal()
    advance_dataset()
    clear_inputs()
    retired_files(list())
    upload_generation(upload_generation() + 1L)
    tryCatch({
      staged <- read_example(id)
      reactive_files$GEX <- staged$GEX
      reactive_files$anno <- staged$anno
      reactive_files$clinic <- staged$clinic
      reactive_files$example_id <- id
      reactive_files$RawCounts <- examples[[id]]$raw
      reactive_files$source <- "example"
      showNotification("Example loaded. Click Preprocess & map.", type = "message")
    }, error = report_error)
    sync_mode()
  }, priority = 120)

  # --- Preprocess & map: works for uploads OR example ---
  observeEvent(input$map, tryCatch({
    invalidate_mapping()
    files <- uploads()
    have_uploads <- identical(reactive_files$source, "upload") &&
      all(lengths(files) > 0L)
    have_example <- identical(reactive_files$source, "example") &&
      !is.null(reactive_files$GEX) && !is.null(reactive_files$clinic) && !is.null(reactive_files$anno)
    if (!have_uploads && !have_example)
      stop("Load the example or upload expression, clinical and annotation files before preprocessing.")
    
    if (have_example) sync_mode()
    withProgress(message = "Running Mapping...", value = 0, {
      incProgress(0.15, detail = "Reading inputs...")
      
      if (have_uploads) {
        reactive_files$GEX <- reactive_files$clinic <- reactive_files$anno <- NULL
        clinic_df <- .read_table(files$clinic, "PatientID")
        stopifnot("PatientID" %in% names(clinic_df))
        rownames(clinic_df) <- clinic_df$PatientID
        anno_df <- .read_table(files$anno, "probe")
        stopifnot(all(c("probe","ENTREZID") %in% names(anno_df)))
        rownames(anno_df) <- anno_df$probe
        reactive_files$GEX <- .read_expression(files$GEX$datapath,
            .upload_separator(files$GEX$name), clinic_df$PatientID, anno_df$probe, raw_mode())
        reactive_files$clinic <- clinic_df
        reactive_files$anno <- anno_df
        reactive_files$RawCounts <- raw_mode()
        showNotification(sprintf("Data type: %s", if (reactive_files$RawCounts) "Raw counts" else "Normalized (log2)"),
                         type = "message")
      }
      
      incProgress(0.5, detail = "Aligning...")
      se_obj <- aligned_input()

      incProgress(0.8, detail = "Mapping()...")
      tryCatch({
        output_text <- capture.output({
          message("Mapping() -> RawCounts = ", reactive_files$RawCounts)
          data_input <- BreastSubtypeR::Mapping(se_obj = se_obj, RawCounts = raw_mode(),
                                impute = FALSE, verbose = TRUE)
          cat("Step 1 is complete.\nYou may now proceed to Step 2.\n")
        })
        reactive_files$data_input <- data_input
        reactive_files$mapped_key <- input_key()
        showNotification(HTML(paste(output_text, collapse = "<br>")),
                         type = "message", duration = 4)
      }, error = report_error)
      
      incProgress(1, detail = "Completed.")
    })
  }, error = report_error))
  
  # ---- RUN subtyping (shortened: same logic you already had) ----
  .validate_has_clinical <- function(se, want) {
    out <- list(use = FALSE, se = se, msg = NULL)
    if (!want || is.null(se)) return(out)
    cd <- as.data.frame(SummarizedExperiment::colData(se))
    missing <- setdiff(c("TSIZE","NODE"), names(cd))
    if (length(missing) > 0) { out$msg <- paste("missing column(s):", paste(missing, collapse = ", ")); return(out) }
    if (!is.numeric(cd$NODE)) cd$NODE <- suppressWarnings(as.numeric(cd$NODE))
    SummarizedExperiment::colData(se)$NODE <- cd$NODE
    # Pass supplied TSIZE values, including factor labels, to the shared decoder.
    if (anyNA(cd$TSIZE) || any(!is.finite(cd$NODE))) {
      out$msg <- "missing TSIZE or nonfinite NODE"
      return(out)
    }
    out$use <- TRUE; out$se <- se; out
  }
  
  observeEvent(input$run, tryCatch(withCallingHandlers({
    invalidate_prediction()
    if (!mapping_current()) stop("Run Step 1 (Preprocess & map) before subtyping.")
    provenance <- list(inputs = input_key(), parameters = prediction_parameters())
    se_NC  <- reactive_files$data_input$x_NC  %||% reactive_files$data_input$se_NC
    se_SSP <- reactive_files$data_input$x_SSP %||% reactive_files$data_input$se_SSP
    want_4 <- k_is_4class()
    if (identical(input$BSmethod, "AIMS")) want_4 <- FALSE
    
    # ---------- AUTO ----------
    if (identical(input$BSmethod, "AUTO Mode")) {
      chk <- auto_requirements()
      if (!isTRUE(chk$ready)) {
        stop("AUTO needs ER and HER2 (ER+/ER-, HER2+/HER2-).")
      }
      withProgress(message = "Performing analysis (AUTO Mode)...", value = 0.5, {
        result_auto <- BreastSubtypeR::BS_Multi(
          data_input  = reactive_files$data_input,
          methods     = "AUTO",
          Subtype     = want_4,
          hasClinical = FALSE
        )
        tab <- if (isTRUE(want_4)) result_auto$res_subtypes.Subtype else result_auto$res_subtypes
        calls_tbl_auto <- data.frame(PatientID = rownames(tab), tab, check.names = FALSE)
        reactive_files$prediction <- c(provenance, list(
          method = "AUTO", want_4 = want_4,
          header = sprintf("Method: AUTO | Classes: %s | ROR: OFF", if (want_4) "4" else "5"),
          calls = calls_tbl_auto, full = NULL, plot = tab))
      })
      return(invisible(NULL))
    }
    
    # ---------- Single methods ----------
    is_nc_method <- input$BSmethod %in% c("PAM50","cIHC","cIHC.itr","PCAPAM50","ssBC")
    want_clin <- isTRUE(input$hasClinical)
    chkclin <- .validate_has_clinical(if (is_nc_method) se_NC else NULL, want_clin)
    se_NC    <- chkclin$se
    use_clin <- want_clin && chkclin$use
    if (want_clin && !use_clin && is_nc_method) {
      showNotification(paste("Clinical variables disabled —", chkclin$msg), type = "warning", duration = 5)
    }
    
    withProgress(message = "Performing analysis...", value = 0.55, {
      res <- NULL
      ssp <- NULL
      
      if (input$BSmethod == "PAM50") {
        params <- provenance$parameters
        cal <- params$calibration
        args <- list(se_obj = se_NC, hasClinical = use_clin, Subtype = want_4)
        if (identical(cal, "Internal")) {
          args$calibration <- "Internal"; args$internal <- params$internal
        } else if (identical(cal, "External")) {
          args$calibration <- "External"
          if (identical(params$external, "Given.mdns")) {
            if (is.null(params$medians)) stop("Upload custom medians before running PAM50 with Given.mdns.")
            medians <- if (tolower(tools::file_ext(params$medians$name)) == "csv") {
              read.csv(params$medians$datapath, header = TRUE, row.names = NULL, sep = ",", check.names = FALSE)
            } else {
              read.table(params$medians$datapath, header = TRUE, row.names = NULL, sep = "\t",
                         fill = TRUE, stringsAsFactors = FALSE, check.names = FALSE, quote = "", comment.char = "")
            }
            args$external <- "Given.mdns"; args$medians <- medians
          } else args$external <- params$external
        } else args$calibration <- "None"
        res <- do.call(BreastSubtypeR::BS_parker, args)
      }
      
      if (input$BSmethod == "cIHC") {
        res <- BreastSubtypeR::BS_cIHC(se_obj = se_NC, Subtype = want_4, hasClinical = use_clin)
      }
      if (input$BSmethod == "cIHC.itr") {
        if (!is.numeric(input$ratio) || length(input$ratio) != 1L ||
            !is.finite(input$ratio) || input$ratio <= 0)
          stop("ER+ : ER- target ratio must be a single finite positive number (default 54 / 64).")
        if (!is.numeric(input$iteration) || length(input$iteration) != 1L ||
            !is.finite(input$iteration) || input$iteration < 10 ||
            input$iteration != floor(input$iteration) || input$iteration > .Machine$integer.max)
          stop("Iterations must be a whole number from 10 to 2147483647.")
        res <- BreastSubtypeR::BS_cIHC.itr(se_obj = se_NC, iteration = input$iteration, ratio = input$ratio,
                                           Subtype = want_4, hasClinical = use_clin)
      }
      if (input$BSmethod == "PCAPAM50") {
        res <- BreastSubtypeR::BS_PCAPAM50(se_obj = se_NC, Subtype = want_4, hasClinical = use_clin)
      }
      if (input$BSmethod == "ssBC") {
        req(se_NC)
        cd <- as.data.frame(SummarizedExperiment::colData(se_NC))
        opt <- as.character(input$s %||% "")[1]
        need_cols <- switch(opt,
                            "ER"    = c("ER"),
                            "ER.v2" = c("ER","HER2"),
                            "TN"    = c("TN"),
                            "TN.v2" = c("TN"),
                            character(0)
        )
        miss <- setdiff(need_cols, names(cd))
        if (length(miss)) {
          showNotification(paste("ssBC requires:", paste(need_cols, collapse=", "),
                                 "| missing:", paste(miss, collapse=", ")), type="error", duration=7)
          return(invisible(NULL))
        }
        withCallingHandlers({
          res <- BreastSubtypeR::BS_ssBC(se_obj = se_NC, s = input$s, Subtype = want_4, hasClinical = use_clin)
        }, warning = function(w) {
          # The experimental warning is shown persistently with its full
          # export; keep unrelated warnings visible as before.
          if (!inherits(w, "bsr_experimental_ror_pc_warning"))
            showNotification(conditionMessage(w), type = "warning", duration = 15)
          invokeRestart("muffleWarning")
        })
      }
      
      if (input$BSmethod == "AIMS") {
        req(se_SSP)
        res_AIMS <- BreastSubtypeR::BS_AIMS(se_obj = se_SSP)
        if (!is.null(res_AIMS$cl) && is.matrix(res_AIMS$cl) && ncol(res_AIMS$cl) >= 1) {
          res$BS.all <- data.frame(
            PatientID = rownames(res_AIMS$cl),
            BS = res_AIMS$cl[, 1, drop = TRUE],
            row.names = rownames(res_AIMS$cl),
            check.names = FALSE
          )
        }
      }
      
      if (input$BSmethod == "sspbc") {
        req(se_SSP)
        model <- if (isTRUE(want_4)) "ssp.subtype" else "ssp.pam50"
        res_sspbc <- BreastSubtypeR::BS_sspbc(se_obj = se_SSP, ssp.name = model)
        ssp <- list(model = model, native = res_sspbc, input_ids = colnames(se_SSP),
                    full_choice = paste0("full-", input$run))
        if (!is.null(res_sspbc$cl) && is.matrix(res_sspbc$cl) && ncol(res_sspbc$cl) >= 1) {
          res$BS.all <- data.frame(
          PatientID = rownames(res_sspbc$cl),
          BS = res_sspbc$cl[, 1, drop = TRUE],
          row.names = rownames(res_sspbc$cl),
          check.names = FALSE
        )
        }
      }
      
      if (is.null(res) || is.null(res$BS.all) || nrow(res$BS.all) == 0L)
        stop("No per-sample subtype results were returned. Check the selected method and inputs.")
      ex <- .build_single_exports(res, want_4)

      # Plots
      out <- NULL
      if (!is.null(res) && !is.null(res$BS.all) && nrow(res$BS.all) > 0) {
        labcol <- if (isTRUE(want_4) && "BS.Subtype" %in% names(res$BS.all)) "BS.Subtype"
        else if ("BS" %in% names(res$BS.all)) "BS"
        else setdiff(names(res$BS.all), "PatientID")[1]
        out <- data.frame(PatientID = res$BS.all$PatientID,
                          Subtype   = res$BS.all[[labcol]],
                          check.names = FALSE)
      } else if (!is.null(res) && !is.null(res$cl) && is.matrix(res$cl)) {
        out <- data.frame(PatientID = rownames(res$cl),
                          Subtype   = res$cl[, 1, drop = TRUE],
                          check.names = FALSE)
      }
      
      if (is.null(out) || nrow(out) == 0) {
        stop("No per-sample subtype table returned. Check the selected method and inputs.")
      } else {
        mat <- if (input$BSmethod %in% c("sspbc", "AIMS")) {
          req(se_SSP); log2(SummarizedExperiment::assay(se_SSP) + 1)
        } else {
          req(se_NC);  SummarizedExperiment::assay(se_NC)
        }
        reactive_files$prediction <- c(provenance, list(
          method = input$BSmethod, want_4 = want_4,
          header = .method_header_line(provenance$parameters, want_4, use_clin),
          calls = ex$calls,
          full = if (input$BSmethod %in% c("AIMS", "sspbc")) NULL else ex$full,
          plot = out, matrix = mat), if (!is.null(ssp)) list(ssp = ssp))
      }
    })
  }, bsr_invalid_tsize_warning = function(w) {
    showNotification(conditionMessage(w), type = "warning", duration = 15)
    invokeRestart("muffleWarning")
  }), error = report_error))
}
