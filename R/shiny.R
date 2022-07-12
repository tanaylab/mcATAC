app_ui <- function(request) {
    css <- "
            .checkbox_inline .shiny-input-container {
            display: inline-block;
            width: auto;
            }
            #inline label{ display: table-cell; text-align: center; vertical-align: middle; }
            #inline .form-group { display: table-row;}
            "
    if (!is.null(shiny_mct@id)) {
        title <- glue("mcATAC viewer: {shiny_mct@id}")
    } else {
        title <- "mcATAC viewer"
    }
    tagList(
        shinyjs::useShinyjs(), # Set up shinyjs
        fluidPage(
            tags$style(type = "text/css", css),
            titlePanel(title),
            fluidRow(
                column(
                    12,
                    wellPanel(
                        fluidRow(
                            align = "left",
                            column(
                                11,
                                tagList(
                                    div(
                                        tags$b("Move:"),
                                        shinyWidgets::actionGroupButtons(
                                            inputIds = c(paste0("shift_left_", rev(c(10, 47.5, 95))), paste0("shift_right_", c(10, 47.5, 95))),
                                            labels = c("<<<", "<<", "<", ">", ">>", ">>>")
                                        ),
                                        tags$b("Zoom in:"),
                                        shinyWidgets::actionGroupButtons(
                                            inputIds = paste0("zoom_in_", c(1.5, 3, 10)),
                                            labels = c("1.5x", "3x", "10x")
                                        ),
                                        tags$b("Zoom out:"),
                                        shinyWidgets::actionGroupButtons(
                                            inputIds = paste0("zoom_out_", c(1.5, 3, 10)),
                                            labels = c("1.5x", "3x", "10x")
                                        ),
                                        span(class = "checkbox_inline", checkboxInput("detect_dca", "Detect DCA", value = TRUE)),
                                        tags$div(numericInput("dca_peak_lf_thresh", "Peak thr:", value = 1.2, min = 0, max = 10, step = 0.1), id = "inline", style = "display:inline-block;vertical-align: middle;margin-left:50px;"),
                                        tags$div(numericInput("dca_trough_lf_thresh", "Trough thr:", value = -1, max = 0, min = -10, step = 0.1), id = "inline", style = "display:inline-block;vertical-align: middle;"),
                                        tags$div(numericInput("dca_sz_frac_for_peak", "Max:", value = 0.25, max = 1, min = 0, step = 0.05), id = "inline", style = "display:inline-block;vertical-align: middle;")
                                    ),
                                    fluidRow(
                                        column(
                                            1,
                                            tags$div(
                                                actionButton("back", "", icon = icon("undo"), style = "margin-top:25px;"),
                                                actionButton("forward", "", icon = icon("redo"), style = "margin-top:25px;margin-left:0px;")
                                            )
                                        ),
                                        column(
                                            7,
                                            shinyWidgets::searchInput(
                                                inputId = "coords",
                                                label = "Enter coordinates:",
                                                placeholder = "chrom:start-end",
                                                btnSearch = icon("circle"),
                                                btnReset = icon("remove"),
                                                width = "100%"
                                            )
                                        ),
                                        column(3, shinyWidgets::virtualSelectInput(
                                            "genes",
                                            "Select gene:",
                                            choices = NULL,
                                            multiple = FALSE,
                                            search = TRUE,
                                            inline = TRUE
                                        )),
                                        column(1, actionButton(
                                            inputId = "update_gene_coord",
                                            label = "Go to gene",
                                            style = "margin-top: 25px; margin-bottom: 7.5px; margin-left: 1px;"
                                        ))
                                    ),
                                    em(textOutput("current_coords", inline = TRUE))
                                )
                            )
                        )
                    )
                )
            ),
            fluidRow(
                shinycssloaders::withSpinner(
                    plotOutput(
                        "region_plot",
                        height = "50vh",
                        brush = brushOpts(
                            id = "region_brush",
                            direction = "x",
                            resetOnNew = TRUE,
                            delay = 1000
                        )
                    )
                )
            )
        )
    )
}

app_server <- function(input, output, session) {
    intervals <- reactiveVal()
    globals <- reactiveValues()

    if (has_rna(shiny_mct) && has_cell_type(shiny_mct) && has_cell_type_colors(shiny_mct)) {
        hc <- mc_hclust_rna(shiny_mct, force_cell_type = TRUE)
    } else {
        hc <- NULL
    }

    promoters <- get_promoters(upstream = 5e4, downstream = 5e4) %>%
        mutate(
            coords = glue("{chrom}:{start}-{end}"),
            label = glue("{geneSymbol} ({coords})")
        )

    observe({
        shinyWidgets::updateVirtualSelect(
            "genes",
            choices = shinyWidgets::prepare_choices(promoters, label, coords)
        )
    })

    observeEvent(input$update_gene_coord, {
        shinyWidgets::updateSearchInput(
            session,
            inputId = "coords",
            value = input$genes,
            trigger = TRUE
        )
    })

    update_intervals <- function(new_intervals) {        
        if (new_intervals$end - new_intervals$start > 1.5e6) {
            showNotification("Region is too large")
            req(FALSE)
        }
        
        globals$history <- c(globals$history, list(new_intervals))
        globals$history_iterator <- length(globals$history)
        intervals(new_intervals)
    }

    observeEvent(input$coords_search, {
        req(input$coords)
        new_intervals <- parse_coordinate_text(input$coords)
        req(new_intervals)
        update_intervals(new_intervals)
    })

    observeEvent(input$back, {
        req(globals$history)
        req(globals$history_iterator)
        if (globals$history_iterator > 1) {
            globals$history_iterator <- globals$history_iterator - 1
            intervals(globals$history[[globals$history_iterator]])
        }
    })

    observeEvent(input$forward, {
        req(globals$history)
        req(globals$history_iterator)
        if (globals$history_iterator < length(globals$history)) {
            globals$history_iterator <- globals$history_iterator + 1
            intervals(globals$history[[globals$history_iterator]])
        }
    })

    output$region_plot <- renderPlot({
        req(intervals())
        req(!is.null(input$detect_dca))
        req(input$dca_peak_lf_thresh)
        req(input$dca_trough_lf_thresh)
        req(input$dca_sz_frac_for_peak)
        mct_plot_region(
            shiny_mct, intervals(),
            detect_dca = input$detect_dca,
            gene_annot = TRUE,
            hc = hc,
            peak_lf_thresh1 = input$dca_peak_lf_thresh,
            trough_lf_thresh1 = input$dca_trough_lf_thresh,
            sz_frac_for_peak = input$dca_sz_frac_for_peak
        )
    }) %>% bindCache(
        intervals(),
        input$detect_dca,
        input$dca_peak_lf_thresh,
        input$dca_trough_lf_thresh,
        input$dca_sz_frac_for_peak
    )

    output$current_coords <- renderText({
        if (is.null(intervals())) {
            return("Please enter valid genomic coordinates (e.g. \"chr3:34300000-35000020\")")
        }
        paste0(intervals()$chrom, ":", intervals()$start, "-", intervals()$end, " (", scales::comma(intervals()$end - intervals()$start), " bp)")
    })

    purrr::walk(c(1.5, 3, 10), ~ {
        observeEvent(input[[glue("zoom_in_{.x}")]], {
            req(intervals())
            update_intervals(gintervals.zoom_in(intervals(), .x))
        })
        observeEvent(input[[glue("zoom_out_{.x}")]], {
            req(intervals())
            update_intervals(gintervals.zoom_out(intervals(), .x))
        })
    })

    purrr::walk(c(10, 47.5, 95), ~ {
        observeEvent(input[[glue("shift_left_{.x}")]], {
            req(intervals())
            shift <- (intervals()$end - intervals()$start) * .x / 100
            update_intervals(gintervals.shift_left(intervals(), shift))
        })

        observeEvent(input[[glue("shift_right_{.x}")]], {
            req(intervals())
            shift <- (intervals()$end - intervals()$start) * .x / 100
            update_intervals(gintervals.shift_right(intervals(), shift))
        })
    })

    observeEvent(input$region_brush, {
        req(intervals())
        xmin <- max(0, input$region_brush$xmin)
        xmax <- min(1, input$region_brush$xmax)
        zoom_intervals <- intervals() %>%
            mutate(
                len = end - start,
                end = round(start + len * xmax),
                start = round(start + len * xmin)
            ) %>%
            select(chrom, start, end)
        update_intervals(zoom_intervals)
    })

    observe({
        shinyjs::toggle(id = "dca_peak_lf_thresh", condition = input$detect_dca)
        shinyjs::toggle(id = "dca_trough_lf_thresh", condition = input$detect_dca)
        shinyjs::toggle(id = "dca_sz_frac_for_peak", condition = input$detect_dca)
    })
}

parse_coordinate_text <- function(text) {
    text <- gsub(",", "", text) # remove commas
    text <- gsub(" +", " ", text) # remove multiple spaces
    text <- gsub(": +", ":", text) # remove spaces after colon
    text <- gsub("- +", "-", text) # remove spaces after hyphen
    coords <- stringr::str_split(text, "[:_ -]")[[1]]
    chrom <- coords[1]
    chrom <- gsub("^chrom", "", chrom)
    chrom <- gsub("^chr", "", chrom)
    chrom_str <- paste0("chr", chrom)
    start <- as.numeric(coords[2])
    end <- as.numeric(coords[3])

    validate(
        need(chrom_str %in% gintervals.all()$chrom, "Chromosome not found or not in the genome"),
        need(!is.na(start), "Start coordinate is not specified or not in a valid format"),
        need(!is.na(end), "End coordinate is not specified or not in a valid format"),
        need(start <= end, "Start coordinate is greater or equal than end coordinate"),
        need(end - start < 1.5e6, "Region is too large")
    )

    intervals <- data.frame(chrom = chrom_str, start = start, end = end)
    fixed_intervals <- gintervals.force_range(intervals)

    return(fixed_intervals)
}



#' Run a shiny app for viewing the MCT object
#'
#' @param mct MCT object
#'
#' @examples
#' \dontrun{
#' run_app(mct)
#' }
#'
#' @export
run_app <- function(mct,
                    port = NULL,
                    host = NULL,
                    launch.browser = FALSE) {
    library(misha)
    library(shiny)
    opt <- options(gmultitasking = FALSE, shiny.usecairo = TRUE)
    withr::defer(options(opt))
    shiny_mct <<- mct
    shiny::shinyApp(
        ui = app_ui,
        server = app_server,
        options = list(port = port, host = host, launch.browser = launch.browser)
    )
}
