app_ui <- function(request) {
    css <- "
            .checkbox_inline .shiny-input-container {
            display: inline-block;
            width: auto;
            }
            "
    if (!is.null(shiny_mct@id)) {
        title <- glue("mcATAC viewer: {shiny_mct@id}")
    } else {
        title <- "mcATAC viewer"
    }
    tagList(
        fluidPage(
            tags$style(css),
            titlePanel(title),
            fluidRow(
                column(
                    12,
                    wellPanel(
                        fluidRow(
                            align = "left",
                            column(
                                10,
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
                                        span(class = "checkbox_inline", checkboxInput("detect_dca", "Detect DCA"))
                                    ),
                                    fluidRow(
                                        column(
                                            8,
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

    observeEvent(input$coords_search, {
        req(input$coords)
        new_intervals <- parse_coordinate_text(input$coords)
        req(new_intervals)
        intervals(new_intervals)
    })

    output$region_plot <- renderPlot({
        req(intervals())
        req(!is.null(input$detect_dca))
        mct_plot_region(shiny_mct, intervals(), detect_dca = input$detect_dca, gene_annot = TRUE)
    }) %>% bindCache(intervals(), input$detect_dca)

    output$current_coords <- renderText({
        if (is.null(intervals())) {
            return("Please enter valid genomic coordinates (e.g. \"chr3:34300000-35000020\")")
        }
        paste0(intervals()$chrom, ":", intervals()$start, "-", intervals()$end, " (", scales::comma(intervals()$end - intervals()$start), " bp)")
    })

    purrr::walk(c(1.5, 3, 10), ~ {
        observeEvent(input[[glue("zoom_in_{.x}")]], {
            req(intervals())
            intervals(gintervals.zoom_in(intervals(), .x))
        })
        observeEvent(input[[glue("zoom_out_{.x}")]], {
            req(intervals())
            intervals(gintervals.zoom_out(intervals(), .x))
        })
    })

    purrr::walk(c(10, 47.5, 95), ~ {
        observeEvent(input[[glue("shift_left_{.x}")]], {
            req(intervals())
            shift <- (intervals()$end - intervals()$start) * .x / 100
            intervals(gintervals.shift_left(intervals(), shift))
        })

        observeEvent(input[[glue("shift_right_{.x}")]], {
            req(intervals())
            shift <- (intervals()$end - intervals()$start) * .x / 100
            intervals(gintervals.shift_right(intervals(), shift))
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
        intervals(zoom_intervals)
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
