# Define UI for iBreastSubtypeR
hosted <- getFromNamespace(".shiny_hosted_profile", "BreastSubtypeR")()
app_theme <- bslib::bs_theme(
  version = 5,
  base_font = bslib::font_collection("system-ui", "Segoe UI", "sans-serif"),
  primary = "#9C245F"
)

ui <- bslib::page_fluid(
  theme = app_theme,
  title = "iBreastSubtypeR",
  
  # --- Global CSS fixes ---
  tags$head(tags$style(HTML("
    .bslib-card, .card, .card-body { overflow: visible !important; }
    .selectize-control, .selectize-dropdown, .dropdown-menu { z-index: 1040 !important; }
    .modal-backdrop { z-index: 1050 !important; }
    .modal { z-index: 1060 !important; }
    body { line-height: 1.5; }
    .container-fluid { max-width: 1440px; padding: 20px clamp(12px, 3vw, 36px); }
    .app-header { display: flex; gap: 18px; align-items: center; padding-bottom: 18px; border-bottom: 1px solid var(--bs-border-color); }
    .app-logo { width: 84px; height: 96px; flex-shrink: 0; }
    .app-header h1 { font-size: clamp(1.6rem, 3vw, 2rem); margin: 0 0 4px; }
    .app-header p { margin: 0 0 6px; }
    .header-meta { display: flex; gap: 8px 16px; flex-wrap: wrap; align-items: center; font-size: .9rem; }
    .research-label { font-weight: 600; color: var(--bs-emphasis-color); }
    .workflow-nav { display: flex; gap: 12px 24px; flex-wrap: wrap; margin: 14px 0 24px; }
    .section-heading { font-size: 1.4rem; margin: 24px 0 12px; scroll-margin-top: 16px; }
    .section-toolbar { display: flex; justify-content: space-between; gap: 12px; flex-wrap: wrap; align-items: center; margin-bottom: 12px; }
    .section-toolbar .section-heading { margin: 0; }
    .card-header { font-weight: 600; }
    .bslib-card .card-body { gap: 12px; }
    .input-essential { font-size: .9rem; margin-bottom: 8px; }
    .help-details { margin: 4px 0 12px; }
    .help-details summary { cursor: pointer; font-weight: 600; padding: 6px 0; }
    .method-help { margin: 6px 0 12px; padding: 10px 12px; border-left: 3px solid var(--bs-border-color); background: var(--bs-tertiary-bg); border-radius: 6px; }
    .method-help ul { margin-bottom: 0; }
    .chip.auto { border: 1px solid; padding: 6px 10px; border-radius: 6px; display: inline-flex; gap: 8px; flex-wrap: wrap; align-items: center; font-size: .9rem; margin-bottom: 12px; }
    .chip.auto.pending { color: var(--bs-body-color); background: var(--bs-tertiary-bg); border-color: var(--bs-border-color); }
    .chip.auto.ready { color: var(--bs-success-text-emphasis); background: var(--bs-success-bg-subtle); border-color: var(--bs-success-border-subtle); }
    .chip.auto.blocked { color: var(--bs-danger-text-emphasis); background: var(--bs-danger-bg-subtle); border-color: var(--bs-danger-border-subtle); }
    .chip-note { font-size: .85rem; }
    .map-actions { display: flex; gap: 12px 20px; align-items: center; flex-wrap: wrap; }
    .map-actions .help-block { margin: 0; }
    .plot-scroll { overflow-x: auto; }
    .plot-scroll-hint { display: none; margin: 0; font-size: .85rem; color: var(--bs-secondary-color); }
    #multi_plot, #heat2 { min-width: 760px; }
    #pie1 { min-width: 360px; }
    :is(a, button, input, select, summary, [tabindex]):focus-visible,
    .selectize-input.focus, .btn-file:focus-within { outline: 3px solid var(--bs-link-color); outline-offset: 3px; }
    .table-responsive:focus-visible, .plot-scroll:focus-visible { outline-offset: -3px; }
    footer { margin: 20px 0 4px; color: var(--bs-secondary-color); }
    @media (max-width: 575px) {
      .plot-scroll-hint { display: block; }
      .app-header { gap: 12px; align-items: flex-start; }
      .app-logo { width: 56px; height: 64px; }
      .header-meta { gap: 6px 12px; }
      .workflow-nav { gap: 8px 16px; }
      .section-toolbar .btn, .map-actions .btn { width: 100%; }
    }
  "))),
  
  # Preserve scroll position helper
  tags$head(
    tags$script(HTML("
      (function () {
        // Only presentation changes while locked. Sending an example selection
        // through the input binding would queue it as a later manual edit.
        Shiny.addCustomMessageHandler('bsrDataType', function(mode) {
          var radios = document.querySelectorAll('input[name=is_raw_counts]');
          radios.forEach(function(radio) {
            radio.checked = radio.value === mode.value;
            radio.disabled = mode.locked;
          });
          if (!mode.locked) Shiny.setInputValue('is_raw_counts', mode.value);
        });
        // Shiny handles drops with a synthetic change. Capture both native
        // entry points before the binding starts its asynchronous transfer.
        function startUpload(input, files) {
          if (input && files && files.length &&
              input.closest('#upload_GEX, #upload_clinic, #upload_anno, #pam50_calibration')) {
            Shiny.setInputValue('upload_started', input.id, {priority: 'event'});
          }
        }
        document.addEventListener('change', function(event) {
          startUpload(event.target, event.target.files);
        }, true);
        document.addEventListener('drop', function(event) {
          var zone = event.target.closest('div.input-group');
          if (zone) startUpload(zone.querySelector('input[type=file]'),
                                event.dataTransfer && event.dataTransfer.files);
        }, true);
        var lastY = 0, restoreOnIdle = false;
        function saveScroll(){ lastY = window.scrollY || window.pageYOffset || 0; }
        function restoreScroll(){ window.scrollTo(0, lastY || 0); }
        function scrollToMap(){ var el = document.getElementById('map'); if (el && el.scrollIntoView) el.scrollIntoView({ behavior:'smooth', block:'center' }); }
        ['GEX','clinic','anno'].forEach(function(id) {
          var field = '#upload_' + id + ' input[type=file]';
          $(document).on('mousedown click focusin', '#upload_' + id + ' label, ' + field, function(){ saveScroll(); setTimeout(restoreScroll, 0); });
          $(document).on('change', field, function(){ setTimeout(function(){ restoreScroll(); if (id==='anno') scrollToMap(); }, 50); });
          $(document).on('focusout', field, function(){ setTimeout(restoreScroll, 0); });
        });
        $(document).on('mousedown click', '#map', function(){ saveScroll(); restoreOnIdle = true; });
        $(document).on('shiny:idle', function(){ if (restoreOnIdle){ setTimeout(function(){ restoreScroll(); scrollToMap(); restoreOnIdle=false; }, 50); }});
      })();
    "))
  ),
  
  tags$header(class = "app-header",
    bslib::card_image("logo.svg", alt = "BreastSubtypeR logo",
      class = "app-logo", border_radius = "none"),
    div(
      h1("iBreastSubtypeR"),
      p("Intrinsic molecular subtyping with published classifiers."),
      div(class = "header-meta",
        span(class = "research-label", "For research use"),
        if (hosted) span("Remote processing on the host. Choose a packaged example. All file uploads are disabled.")
        else span("The application runs locally in your R session."),
        actionLink("load_example", "Load example data..."),
        actionLink("about", "About & Citation")
      )
    )
  ),
  tags$nav(class = "workflow-nav", `aria-label` = "Workflow",
    tags$a(href = "#step1", "Data and preprocessing"),
    tags$a(href = "#step2", "Analysis settings"),
    tags$a(href = "#results", "Results and downloads")
  ),
  
  #### Step 1
  div(class = "section-toolbar",
    h2(id = "step1", class = "section-heading", "1. Data and preprocessing"),
    actionButton("load_example_step1", "Load example data...", icon = icon("database"))
  ),
  div(role = "status", textOutput("active_source")),
  bslib::layout_columns(
    col_widths = bslib::breakpoints(sm = c(12, 12, 12), lg = c(4, 4, 4)),
    
    # GEX
    bslib::card(
      bslib::card_header("1) Gene expression (GEX)"),
      uiOutput("upload_GEX"),
      radioButtons("is_raw_counts", "Data type",
                   choices  = list("normalized (log2)" = "norm", "Raw counts (RNA-seq only)" = "raw"),
                   selected = "norm", inline = TRUE),
      p(class = "input-essential", "Genes in rows; samples in columns. Gene IDs must match annotation probe IDs; sample IDs must match clinical PatientID."),
      tags$details(class = "help-details",
        tags$summary("Expression format and processing"), uiOutput("gex_help"))
    ),
    
    # Clinical
    bslib::card(
      bslib::card_header("2) Clinical metadata"),
      uiOutput("upload_clinic"),
      p(class = "input-essential", "Required: PatientID. AUTO requires ER and HER2, coded ER+/ER- and HER2+/HER2-. Other requirements depend on the method."),
      tags$details(class = "help-details",
        tags$summary("Clinical fields and research ROR"), uiOutput("clin_help"))
    ),
    
    # Annotation
    bslib::card(
      bslib::card_header("3) Feature annotation"),
      uiOutput("upload_anno"),
      p(class = "input-essential", "Required: probe and ENTREZID. For raw counts, also supply gene Length in base pairs."),
      tags$details(class = "help-details",
        tags$summary("Annotation format"), uiOutput("anno_help"))
    )
  ),

  bslib::card(
    div(class = "map-actions",
      actionButton("map", "Preprocess & map", icon = icon("sliders"), class = "btn-primary"),
      helpText("Runs Mapping() to align IDs and normalize (as needed). When complete, continue to Step 2.")
    )
  ),
  uiOutput("input_coverage"),
  
  #### Step 2
  h2(id = "step2", class = "section-heading", "2. Analysis settings"),
  
  bslib::card(
    style = "overflow: visible;",
    
    bslib::layout_columns(col_widths = bslib::breakpoints(sm = c(12, 12), lg = c(7, 5)),
    selectizeInput(
      "BSmethod", "Subtyping method",
      choices = list(
        "AUTO - select methods for this cohort" = "AUTO Mode",
        "PAM50 (parker.original | genefu.scale | genefu.robust)" = "PAM50",
        "cIHC" = "cIHC",
        "cIHC.itr" = "cIHC.itr",
        "PCAPAM50" = "PCAPAM50",
        "ssBC" = "ssBC",
        "AIMS" = "AIMS",
        "SSPBC" = "sspbc"
      ),
      selected = "AUTO Mode", width = "100%",
      options = list(openOnFocus = TRUE, dropdownParent = "body")
    ),
    
    radioButtons(
      "k_subtypes", "Subtype classes",
      choices = list("5 classes (includes Normal-like)" = "5",
                     "4 classes (excludes Normal-like)" = "4"),
      selected = "5", inline = TRUE
    )),
    uiOutput("auto_chip"),
    uiOutput("auto_preflight"),
    
    uiOutput("method_help"),
    tags$details(class = "help-details",
      tags$summary("Method guide"),
      tags$ul(
        tags$li(tags$b("AUTO: "), "selects methods according to implemented receptor-annotation and subgroup-size rules. It does not establish complete model-input adequacy or clinical validity."),
        tags$li(tags$b("NC-based: "), "PAM50 (parker.original, genefu.scale, genefu.robust), cIHC, cIHC.itr, PCAPAM50, ssBC and ssBC.v2."),
        tags$li(tags$b("SSP-based: "), "AIMS and SSPBC. AIMS supports five classes only."),
        tags$li("Choose four classes to exclude Normal-like, or five classes to include it. Supported NC methods provide research ROR scores. SSPBC full exports include categorical SSP.ROR bins and their model probabilities. Incomplete SSP.ROR rule coverage is reported for the full export.")
      )
    ),
    
    conditionalPanel(
      condition = "input.BSmethod == 'PAM50' || input.BSmethod == 'cIHC' || input.BSmethod == 'cIHC.itr' || input.BSmethod == 'PCAPAM50' || input.BSmethod == 'ssBC'",
      checkboxInput("hasClinical", "Use clinical variables (ROR)", value = FALSE)
    ),
    
    conditionalPanel(
      condition = "input.BSmethod == 'PAM50'",
      htmltools::tagAppendChild(uiOutput("pam50_calibration"),
        getFromNamespace(".shiny_pam50_controls", "BreastSubtypeR")(allow_custom = !hosted))
    ),
    
    uiOutput("calib_help"),
    
    conditionalPanel(
      condition = "input.BSmethod == 'cIHC.itr'",
      bslib::layout_columns(
        col_widths = c(6, 6),
        numericInput("iteration", label = "Iterations", value = 100, min = 10, step = 10,
                     max = if (hosted) getFromNamespace(".shiny_iteration_limit", "BreastSubtypeR")() else NA),
        tagList(
          numericInput("ratio", label = "ER+ : ER- target ratio", value = 54/64,
                       min = 5e-324, step = "any"),
          helpText("Target number of ER-positive samples per ER-negative sample. Use a finite positive number; values above 1 are allowed. Default: 54 / 64.")
        )
      )
    ),
    
    conditionalPanel(
      condition = "input.BSmethod == 'ssBC'",
      bslib::layout_column_wrap(
        selectInput("s", "Subgroup",
                    choices = list("ER" = "ER", "ER.v2" = "ER.v2", "TN" = "TN", "TN.v2" = "TN.v2"),
                    selected = "ER.v2")
      )
    ),
    
    actionButton("run", "Run subtyping", icon = icon("play-circle"), class = "btn-primary")
  ),
  
  uiOutput("analysis_error"),

  h2(id = "results", class = "section-heading", "Results and downloads"),
  p("Compare subtype calls in the plots below. Download current results as tab-separated text (TSV)."),
  # Visualization
  uiOutput("plotSection"),
  
  # Export
  bslib::card(
    bslib::layout_columns(
      col_widths = bslib::breakpoints(sm = c(12, 12), lg = c(8, 4)),
      uiOutput("export_selector"),
      div(style = "text-align:center;",
          uiOutput("download_control")
      )
    ),
    uiOutput("experimental_ror_pc_notice"),
    uiOutput("ssp_ror_notice")
  ),
  
  tags$hr(),
  tags$footer(
    class = "text-center",
    tags$small(
      HTML("&copy; Karolinska Institutet — BreastSubtypeR. DOI: ",
           "<a href='https://doi.org/10.1093/nargab/lqaf131' target='_blank' rel='noopener noreferrer'>10.1093/nargab/lqaf131</a>")
    )
  )
)
