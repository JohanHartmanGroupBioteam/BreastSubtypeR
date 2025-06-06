# Define UI for iBreastSubtypeR

ui <- bslib::page_fluid(
    titlePanel(
        "iBreastSubtypeR: Interactive Breast Cancer Intrinsic Subtyping"
    ),
    bslib::card(
        tags$div(
            style = "display: flex; align-items: center;",
            # Flexbox container
            tags$div(
                style = "flex: 1; text-align: left; margin-right: 20px;",
                bslib::card_image("logo.svg", height = "180px")
            ),
            tags$div(
                style = "flex: 3; text-align: center;",
                # Right-aligned text
                tags$div(
                    style = "margin-bottom: 10px;",
                    bslib::card_header("Welcome to iBreastSubtypeR!")
                ),
                "This tool integrates intrinsic molecular subtyping methods
        for breast cancer,including nearest-centroid (NC-based)
        and single-sample predictor (SSP-based) approaches.
        It employs standardized input and output formats,
        providing a unified framework
        that is highly compatible with other R packages
        in the gene expression profiling field.",
                tags$div(
                    style = "margin-top: 10px;",
                    bslib::card_footer("Enjoy your subtyping journey!")
                )
            )
        )
    ),

    #### input your data ####
    # Main page section with a title
    h3("Step 1: Please enter your data"),
    # Add a title here

    ## main page, input gene expression etc.
    bslib::layout_column_wrap(
        col_width = 3,
        ## gene expression
        bslib::card(
            bslib::card_header("Input 1: Gene expression table"),
            fileInput(
                "GEX",
                "choose file",
                accept = c(
                    "text/csv",
                    "text/comma-separated-values,text/plain",
                    ".csv",
                    ".txt"
                )
            )
        ),

        ## clinic information
        ## csv/text
        bslib::card(
            bslib::card_header("Input 2: Clinical table"),
            fileInput(
                "clinic",
                "choose file",
                accept = c(
                    "text/csv",
                    "text/comma-separated-values,text/plain",
                    ".csv",
                    ".txt"
                )
            ),
        ),

        ## feature information
        ## csv/text
        bslib::card(
            bslib::card_header("Input 3: Feature table"),
            fileInput(
                "anno",
                "choose file",
                accept = c(
                    "text/csv",
                    "text/comma-separated-values,text/plain",
                    ".csv",
                    ".txt"
                )
            ),
        )
    ),

    #### Map button ####
    bslib::card(actionButton("map", "Map Now", icon = icon("map"))),
    ## layout 3

    #### Please select your method
    # Main page section with a title

    h3("Step 2: Please select method & parameters for subtyping"),
    # Add a title here

    ## options
    bslib::card(
        selectInput(
            "BSmethod",
            "Select subtyping method",
            choices = list(
                "PAM50.parker" = "PAM50.parker",
                "cIHC" = "cIHC",
                "cIHC.itr" = "cIHC.itr",
                "PCAPAM50" = "PCAPAM50",
                "ssBC" = "ssBC",
                "AIMS" = "AIMS",
                "sspbc" = "sspbc"
            ),
            selected = "PAM50.parker"
        ),

        # Conditional panels for each choice
        conditionalPanel(
            condition = "input.BSmethod == 'PAM50.parker' ",
            ## parker

            bslib::layout_column_wrap(
                selectInput(
                    "hasClinical",
                    "Has Clinical",
                    choices = list("FALSE" = "FALSE", "TRUE" = "TRUE"),
                    selected = "FALSE"
                ),
                selectInput(
                    "calibration",
                    "Select calibration",
                    choices = list(
                        "None" = "None",
                        "External" = "External",
                        "Internal" = "Internal"
                    ),
                    selected = "Internal"
                )
            ),
            bslib::card(
                conditionalPanel(
                    condition = "input.calibration == 'External'",
                    bslib::layout_column_wrap(selectInput(
                        "external",
                        "External",
                        choices = list(
                            "Given.mdns" = "Given.mdns",
                            "RNAseq.V2" = "RNAseq.V2",
                            "RNAseq.V1" = "RNAseq.V1",
                            "GC.4x44Kcustom" = "GC.4x44Kcustom",
                            "Agilent_244K" = "Agilent_244K",
                            "commercial_1x44k_postMeanCollapse_WashU" =
                                "commercial_1x44k_postMeanCollapse_WashU",
                            "commercial_4x44k_postMeanCollapse_WashU_v2" =
                                "commercial_4x44k_postMeanCollapse_WashU_v2",
                            "htp1.5_WU_update" = "htp1.5_WU_update",
                            "arrayTrain_postMeanCollapse" =
                                "arrayTrain_postMeanCollapse"
                        ),
                        selected = "RNAseq.V2"
                    ))
                ),

                # Show a new file input field when calibration is "External"
                # and external is "Given.mdns"
                conditionalPanel(
                    condition = "input.calibration == 'External' &&
                    input.external == 'Given.mdns'",
                    fileInput(
                        "medians",
                        "Upload Given.mdns File",
                        accept = c(
                            "text/csv",
                            "text/comma-separated-values,text/plain",
                            ".csv",
                            ".txt"
                        )
                    )
                ),
                conditionalPanel(
                    condition = "input.calibration == 'Internal'",
                    selectInput(
                        "internal",
                        "Internal",
                        choices = list(
                            "parker.original" = "medianCtr",
                            "genefu.scale" = "meanCtr",
                            "genefu.robust" = "qCtr"
                        ),
                        selected = "parker.original"
                    )
                )
            )
        ),
        conditionalPanel(
            condition = "input.BSmethod == 'cIHC'",
            bslib::layout_column_wrap(selectInput(
                "hasClinical",
                "Has Clinical",
                choices = list("FALSE" = "FALSE", "TRUE" = "TRUE"),
                selected = "FALSE"
            ))
        ),
        conditionalPanel(
            condition = "input.BSmethod == 'cIHC.itr'",
            bslib::layout_column_wrap(
                numericInput("iteration", label = "Iteration", value = 100),
                numericInput("ratio", label = "Ratio", value = 54 / 64),
                selectInput(
                    "hasClinical",
                    "Has Clinical",
                    choices = list("FALSE" = "FALSE", "TRUE" = "TRUE"),
                    selected = "FALSE"
                )
            )
        ),
        conditionalPanel(
            condition = "input.BSmethod == 'PCAPAM50'",
            bslib::layout_column_wrap(selectInput(
                "hasClinical",
                "Has Clinical",
                choices = list("FALSE" = "FALSE", "TRUE" = "TRUE"),
                selected = "FALSE"
            ))
        ),
        conditionalPanel(
            condition = "input.BSmethod == 'ssBC'",
            bslib::layout_column_wrap(
                selectInput(
                    "s",
                    "Subgroup",
                    choices = list(
                        "ER" = "ER",
                        "TN" = "TN",
                        "ER.v2" = "ER.v2",
                        "HER2+" = "HER2+",
                        "TNBC" = "TNBC"
                    ),
                    selected = "ER.v2"
                ),
                selectInput(
                    "hasClinical",
                    "Has Clinical",
                    choices = list("FALSE" = "FALSE", "TRUE" = "TRUE"),
                    selected = "FALSE"
                )
            )
        ),
        conditionalPanel(condition = "input.BSmethod == 'AIMS'", ),
        conditionalPanel(
            condition = "input.BSmethod == 'sspbc'",
            bslib::layout_column_wrap(
                selectInput(
                    "Subtype",
                    "Subtype method (TRUE: 4-subtypes; FALSE: 5-subtypes)",
                    choices = list("FALSE" = "FALSE", "TRUE" = "TRUE"),
                    selected = "FALSE"
                )
            )
        ),
        bslib::card(actionButton("run", "Subtype Now", icon = icon("cog")))
    ),
    ## card option

    # Visualization layout placed separately and conditionally rendered
    uiOutput("plotSection"),

    # Footer with download and reset buttons
    fluidRow(
        column(12,
            align = "center",
            div(
                style = "margin-top: 20px; margin-bottom: 20px;",
                downloadButton(
                    "download",
                    "Download Results",
                    style =
                        "width: 220px;
                       height: 40px;
                       background-color: #FF69B4;
                       color: white;
                       border: none;
                       text-align: center;
                       line-height: 20px;"
                )
            )
        )
    )
)
