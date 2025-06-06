# Define server logic

server <- function(input, output) {
    output$logo <- renderImage(
        {
            list(src = "logo.svg", height = "100%", inline = FALSE)
        },
        deleteFile = FALSE
    )

    #### upload data ####

    # Reactive values to store uploaded data and results
    reactive_files <- shiny::reactiveValues(GEX = NULL, clinic = NULL, anno = NULL, data_input = NULL, output_res = NULL)

    shiny::observeEvent(input$map, {
        req(input$GEX, input$clinic, input$anno) # Ensure files are uploaded

        withProgress(message = "Running Mapping...", value = 0, {
            # Update progress bar for each step
            incProgress(0.2, detail = "Loading gene expression data...")
            inFile <- input$GEX
            fileExt <- tools::file_ext(inFile$name)

            reactive_files$GEX <- if (fileExt == "csv") {
                read.csv(inFile$datapath, header = TRUE, row.names = 1, sep = ",")
            } else if (fileExt == "txt") {
                read.table(inFile$datapath, header = TRUE, row.names = 1, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
            } else {
                stop("Unsupported file type")
            }

            incProgress(0.4, detail = "Loading clinical data...")
            inFile <- input$clinic
            fileExt <- tools::file_ext(inFile$name)

            reactive_files$clinic <- if (fileExt == "csv") {
                read.csv(inFile$datapath, header = TRUE, row.names = NULL, sep = ",")
            } else if (fileExt == "txt") {
                read.table(inFile$datapath, header = TRUE, row.names = NULL, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
            } else {
                stop("Unsupported file type")
            }
            rownames(reactive_files$clinic) <- reactive_files$clinic$PatientID

            incProgress(0.6, detail = "Loading feature annotation data...")
            inFile <- input$anno
            fileExt <- tools::file_ext(inFile$name)

            reactive_files$anno <- if (fileExt == "csv") {
                read.csv(inFile$datapath, header = TRUE, row.names = NULL, sep = ",")
            } else if (fileExt == "txt") {
                read.table(inFile$datapath, header = TRUE, row.names = NULL, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
            } else {
                stop("Unsupported file type")
            }

            ## Run Mapping function
            samples <- intersect(colnames(reactive_files$GEX), reactive_files$clinic$PatientID)
            probeID <- intersect(rownames(reactive_files$GEX), reactive_files$anno$probe)
            rownames(reactive_files$anno) = reactive_files$anno$probe
            
            incProgress(0.8, detail = "Running Mapping function...")
            tryCatch(
                {
                    # Redirect console output to a variable
                    se_obj = SummarizedExperiment::SummarizedExperiment(assays = list(counts = reactive_files$GEX[probeID, samples]),
                                                  rowData = reactive_files$anno[probeID,],
                                                  colData = reactive_files$clinic[samples,]
                                                  )
                    
                    output_text <- capture.output({
                        data_input <- Mapping(se_obj = se_obj, impute = TRUE, verbose = TRUE)

                        cat("Step 1 is complete.", sep = "\n")
                        cat("You may now proceed to Step 2.", sep = "\n")
                    })

                    # Store results
                    reactive_files$data_input <- data_input

                    # Display output text as notification
                    showNotification(HTML(paste(output_text, collapse = "<br>")), type = "message", duration = NULL)
                },
                error = function(e) {
                    showNotification(paste("Error:", e$message), type = "error", duration = NULL)
                }
            )

            incProgress(1, detail = "Completed.")
        })
    })

    #### perform analysis ####

    shiny::observeEvent(input$run, {
        # Check if data has been mapped
        req(reactive_files$data_input)

        withProgress(message = "Performing analysis...", value = 0, {
            incProgress(0.2, detail = "Initializing the method...")

            # Execute the selected method
            if (input$BSmethod == "PAM50.parker") {

                args <- list(
                    se_obj = reactive_files$data_input$se_NC,
                    calibration = input$calibration,
                    hasClinical = input$hasClinical
                )

                # Modify the list of arguments based on conditions
                if (input$calibration == "Internal") {
                    args$internal <- input$internal
                } else if (input$calibration == "External" & input$external == "Given.mdns") {
                    req(input$medians)
                    inFile <- input$medians
                    fileExt <- tools::file_ext(inFile$name)

                    # Read the file based on its extension
                    medians <- if (fileExt == "csv") {
                        read.csv(inFile$datapath, header = TRUE, row.names = NULL, sep = ",")
                    } else if (fileExt == "txt") {
                        read.table(inFile$datapath, header = TRUE, row.names = NULL, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
                    } else {
                        stop("Unsupported file type")
                    }

                    args$external <- input$external
                    args$medians <- medians
                } else if (input$calibration == "External" & input$external != "Given.mdns") {
                    args$external <- input$external
                }


                incProgress(0.5, detail = "Running PAM50.parker method...")
                res <- do.call(BreastSubtypeR::BS_parker, args)

                output_res <- res$score.ROR
            }

            if (input$BSmethod == "cIHC") {
                incProgress(0.5, detail = "Running cIHC method...")
                res <- BreastSubtypeR::BS_cIHC(se_obj = reactive_files$data_input$se_NC, hasClinical = input$hasClinical)

                output_res <- res$score.ROR
            }

            if (input$BSmethod == "cIHC.itr") {
                incProgress(0.5, detail = "Running cIHC.itr method...")
                res <- BreastSubtypeR::BS_cIHC.itr(se_obj = reactive_files$data_input$se_NC, iteration = input$iteration, ratio = input$ratio, hasClinical = input$hasClinical)

                output_res <- res$score.ROR
            }

            if (input$BSmethod == "PCAPAM50") {
                incProgress(0.5, detail = "Running PCAPAM50 method...")
                res <- BreastSubtypeR::BS_PCAPAM50(se_obj = reactive_files$data_input$se_NC, hasClinical = input$hasClinical)

                output_res <- res$score.ROR
            }

            if (input$BSmethod == "ssBC") {
                incProgress(0.5, detail = "Running ssBC method...")
                res <- BreastSubtypeR::BS_ssBC(se_obj = reactive_files$data_input$se_NC, s = input$s, hasClinical = input$hasClinical)

                output_res <- res$score.ROR
            }

            if (input$BSmethod == "AIMS") {
                incProgress(0.5, detail = "Running AIMS method...")

                data("BreastSubtypeRobj", package = "BreastSubtypeR")
                res <- BreastSubtypeR::BS_AIMS(se_obj = reactive_files$data_input$se_SSP)

                res$BS.all <- data.frame(
                    PatientID = rownames(res$cl),
                    BS = res$cl[, 1]
                )

                output_res <- res$BS.all
            }

            if (input$BSmethod == "sspbc") {
                incProgress(0.5, detail = "Running sspbc method...")

                res_sspbc <- BreastSubtypeR::BS_sspbc(reactive_files$data_input$se_SSP, ssp.name = "ssp.pam50")

                BS.all <- data.frame(PatientID = rownames(res_sspbc), BS = res_sspbc[, 1], row.names = rownames(res_sspbc))

                if (input$Subtype == "TRUE") {
                    res_sspbc.Subtype <- BreastSubtypeR::BS_sspbc(reactive_files$data_input$se_SSP, ssp.name = "ssp.subtype")

                    BS.all$BS.Subtype <- res_sspbc.Subtype[, 1]
                }
                res <- list()
                res$BS.all <- BS.all
                output_res <- BS.all
            }


            incProgress(0.9, detail = "Finalizing the results...")
            ## Allow downloading the result
            reactive_files$output_res <- output_res
            output$download <- downloadHandler(
                filename = function() {
                    paste("results-", input$BSmethod, ".txt", sep = "")
                },
                content = function(file) {
                    write.table(reactive_files$output_res, file, row.names = F, sep = "\t", quote = F)
                }
            )

            ## Allow visualization
            out <- data.frame(
                PatientID = res$BS.all$PatientID,
                Subtype = if (input$BSmethod == "sspbc" && input$Subtype == "TRUE") {
                    res$BS.all$BS.Subtype
                } else {
                    res$BS.all$BS
                }
            )
            output$pie1 <- renderPlot(BreastSubtypeR::Vis_pie(out))

            matrix <- if (input$BSmethod == "sspbc" || input$BSmethod == "AIMS") {
                log2(SummarizedExperiment::assay(reactive_files$data_input$se_SSP) + 1) %>% as.matrix()
            } else {
                SummarizedExperiment::assay(reactive_files$data_input$se_NC) %>% as.matrix()
            }

            if (input$BSmethod == "sspbc" || input$BSmethod == "AIMS") {
                rownames(matrix) <- NULL
            }

            output$heat2 <- renderPlot(BreastSubtypeR::Vis_heatmap(matrix, out = out))

            output$plotSection <- renderUI({
                req(input$run) # Wait until the run button is clicked

                # Once 'Run' is clicked, show the plot layout
                bslib::layout_columns(
                    col_width = 2,
                    bslib::card(plotOutput("pie1")),
                    bslib::card(plotOutput("heat2"))
                )
            })

            incProgress(1, detail = "Analysis is complete.")
        })

        # Show a notification when the analysis is complete
        # showNotification(  paste( "Analysis complete" ,"You can exit or continue to run another subtyping method.", sep = "\n") , type = "message", duration = 5)
        showNotification(HTML(paste(c(
            "Analysis is complete.",
            "Please review your results."
        ), collapse = "<br>")), type = "message", duration = 5)
    })
}
