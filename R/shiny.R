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

    promoters <- misha.ext::get_promoters(upstream = 5e4, downstream = 5e4) %>%
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

        if (is.null(globals$history_iterator)) {
            globals$history <- list(new_intervals)
            globals$history_iterator <- 1
        } else {
            globals$history <- c(globals$history[1:globals$history_iterator], list(new_intervals))
            globals$history_iterator <- length(globals$history)
        }

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
    }) %>%
        bindCache(
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
    shiny_mct <<- mct
    gset_genome(mct@genome)
    shiny::shinyApp(
        ui = app_ui,
        server = app_server,
        options = list(port = port, host = host, launch.browser = launch.browser)
    )
}

#' Create a bundle for the shiny app in order to run it with shiny server
#'
#' Generate a 'deployment ready' bundle of the shiny app
#'
#' Create a minimal shiny app in \code{path} directory which would contain:
#' \itemize{
#' \item{}{app.R file. }
#' \item{}{mct.rds with the mct object. }
#' }
#'
#' The bundle can then be deployed in shiny-server, shinyapps.io or any other environment that supports serving shiny apps.
#'
#' Note: when deploying to these services - make sure you have the mcATAC package installed with the 'Suggests' dependencies.
#'
#' @param mct a MCTracks object
#' @param path Path to the bundle directory
#' @param overwrite overwrite bundle if already exists
#' @param self_contained include the source code of \code{mcATAC} in the bundle
#' and use it to run the app. Use this in order to ensure that the package would always
#' run the same way, regardless of mcATAC changes. When this option is FALSE,
#' the installed version of \code{mcATAC} would be loaded, which can be occasionally
#' updated for all the \code{mcATAC} apps running from a server. By default, the code
#' of the latest \code{mcATAC} release would be used, see \code{branch} for
#' other options.
#' @param branch name of the \code{mcATAC} branch to include when \code{self_contained=TRUE}. By default, the master branch would be used.
#' @param restart add a file named 'restart.txt' to the bundle. This would force shiny-server to restart the app when updated.
#' @param permissions change the file permissions of the bundle after creation, e.g. "777". When NULL -
#' permissions would not be changed.
#'
#' @inheritDotParams gert::git_clone
#'
#' @examples
#' \dontrun{
#' create_bundle(mct, "/path/to/the/bundle/directory")
#' }
#'
#' @export
create_bundle <- function(mct, path, overwrite = FALSE, self_contained = FALSE, branch = "master", restart = overwrite, permissions = NULL, ...) {
    bundle_dir <- path
    if (fs::dir_exists(bundle_dir)) {
        if (overwrite) {
            fs::dir_delete(bundle_dir)
            fs::dir_create(bundle_dir)
            cli::cli_li("Removing previous bundle ({.field overwrite = TRUE})")
        } else {
            cli::cli_abort("{.path {bundle_dir}} already exists. Run with {.code overwrite=TRUE} to force overwriting it.")
        }
    } else {
        fs::dir_create(bundle_dir)
    }

    readr::write_rds(mct, fs::path(bundle_dir, "mct.rds"))
    fs::file_copy(system.file("app.R", package = "mcATAC"), fs::path(bundle_dir, "app.R"))

    if (self_contained) {
        cli::cli_alert("Creating a self-contained bundle")
        code_dir <- fs::path(bundle_dir, "code")
        if (!is.null(branch) && branch == "latest_release") {
            gert::git_clone("git@github.com:tanaylab/mcATAC", path = code_dir, ...)
            tag_list <- gert::git_tag_list(repo = code_dir)
            latest_tag <- tail(tag_list, n = 1)
            gert::git_branch_create(
                branch = latest_tag$name,
                ref = latest_tag$commit,
                repo = code_dir,
                checkout = TRUE
            )
            cli::cli_alert_info("Using latest release: {.file {latest_tag$name}}")
        } else {
            gert::git_clone("git@github.com:tanaylab/mcATAC", path = code_dir, branch = branch, ...)
        }
    }


    if (restart) {
        fs::file_touch(fs::path(bundle_dir, "restart.txt"))
        cli::cli_li("Adding a file called {.field restart.txt}")
    }

    if (!is.null(permissions)) {
        fs::file_chmod(c(bundle_dir, fs::dir_ls(bundle_dir, recurse = TRUE)), mode = permissions)
        cli::cli_li("Changing permissions to {.field {permissions}}")
    }

    cli::cli_li("Bundle files:")
    fs::dir_tree(bundle_dir)
    cli::cat_line("")
    cli::cli_alert_success("created a bundle at {bundle_dir}")
    cli::cli_li("To deploy to shinyapps.io, run: {.field rsconnect::deployApp(appDir = \"{as.character(bundle_dir)}\")}")
    cli::cli_li("To deploy to another shiny-server service, upload {.path {bundle_dir}} to the service.")
}
