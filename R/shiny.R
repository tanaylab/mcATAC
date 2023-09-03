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
                    8,
                    wellPanel(
                        fluidRow(
                            align = "left",
                            column(
                                12,
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
                                        shinyjs::disabled(
                                            span(class = "checkbox_inline", checkboxInput("detect_dca", "Detect DCA", value = has_rna(shiny_mct) && has_cell_type(shiny_mct) && has_cell_type_colors(shiny_mct)))
                                        ),
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
                                            3,
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
                                        )),
                                        column(
                                            1,
                                            shinyWidgets::materialSwitch(
                                                inputId = "show_controls",
                                                label = "Advanced",
                                                value = FALSE
                                            )
                                        )
                                    ),
                                    fluidRow(
                                        column(
                                            5,
                                            em(textOutput("current_coords", inline = TRUE)),
                                        ),
                                        column(
                                            3,
                                            tags$div(numericInput("n_smooth", "Smooth (bp):", value = 200, min = shiny_mct@resolution, step = shiny_mct@resolution), id = "inline", style = "display:inline-block;vertical-align: middle;", width = "100%"),
                                        ),
                                        column(
                                            2,
                                            tags$div(numericInput("min_color", "Min color:", value = 6, min = 0, step = 1), id = "inline", style = "display:inline-block;vertical-align: middle;", width = "15%"),
                                        ),
                                        column(
                                            2,
                                            tags$div(numericInput("max_color", "Max color:", value = 24, min = 0, step = 1), id = "inline", style = "display:inline-block;vertical-align: middle;", width = "15%")
                                        ),
                                        column(
                                            2,
                                            tags$div(numericInput("min_color2", "Min color:", value = 6, min = 0, step = 1), id = "inline", style = "display:inline-block;vertical-align: middle;", width = "15%"),
                                        ),
                                        column(
                                            2,
                                            tags$div(numericInput("max_color2", "Max color:", value = 24, min = 0, step = 1), id = "inline", style = "display:inline-block;vertical-align: middle;", width = "15%")
                                        )
                                    )
                                )
                            )
                        )
                    )
                ),
                column(
                    2,
                    plotOutput(
                        "pheatmap_plot",
                        height = "300px",
                        width = "300px",
                        fill = FALSE,
                    )
                ),
                column(
                    2,
                    plotOutput(
                        "scatter_plot",
                        height = "300px",
                        width = "300px",
                        fill = FALSE,
                    )
                )

            ),
            fluidRow(
                shinycssloaders::withSpinner(
                    plotOutput(
                        "region_plot",
                        height = "400px",
                        width = "2060px",
                        fill = FALSE,
                        hover = hoverOpts(
                            id = "region_hover",
                            delayType = "throttle",
                            delay = 250
                        ),
                        brush = brushOpts(
                            id = "region_brush",
                            direction = "x",
                            resetOnNew = TRUE,
                            delay = 1000
                        )
                    )
                ),
                shinycssloaders::withSpinner(
                    plotOutput(
                        "region_plot_cor",
                        height = "20px",
                        width = "2067px",
                        fill = FALSE,
                        click = "region_plot_cor_click"
                    )
                ),
                shinycssloaders::withSpinner(
                    plotOutput(
                        "comparison_plot",
                        width = "2070px",
                        height = "200px",
                        fill = FALSE
                    )
                ),
                shinycssloaders::withSpinner(
                    plotOutput(
                        "region_plot_cor2",
                        width = "2067px",
                        height = "20px",
                        fill = FALSE,
                        click = "region_plot_cor2_click"
                    )
                ),
                shinycssloaders::withSpinner(
                    plotOutput(
                        "region_plot2",
                        width = "2060px",
                        height = "400px",
                        fill = FALSE
                    )
                )
            ),
            uiOutput("hover_info", style = "pointer-events: none")
        )
    )
}

app_server <- function(input, output, session) {
    intervals <- reactiveVal()
    intervals2 <- reactiveVal()
    intervals_comparison <- reactiveVal()
    globals <- reactiveValues()
    mat <- reactiveVal()
    mat2 <- reactiveVal()
    scatter_data <- reactiveVal()
    comb_mat <- reactiveVal()
    comb_mat_txt <- reactiveVal()

    if (!is.null(shiny_hc)) {
        hc <- shiny_hc
        shinyjs::enable("detect_dca") # Enable DCA detection
    } else if (has_rna(shiny_mct) && has_cell_type(shiny_mct) && has_cell_type_colors(shiny_mct)) {
        hc <- mc_hclust_rna(shiny_mct, force_cell_type = TRUE)
        shinyjs::enable("detect_dca") # Enable DCA detection
    } else {
        hc <- NULL
    }

    if (!is.null(shiny_mct2)) {
        if (!is.null(shiny_hc2)) {
            hc2 <- shiny_hc2
            shinyjs::enable("detect_dca") # Enable DCA detection
        } else if (has_rna(shiny_mct2) && has_cell_type(shiny_mct2) && has_cell_type_colors(shiny_mct2)) {
            hc2 <- mc_hclust_rna(shiny_mct, force_cell_type = TRUE)
            shinyjs::enable("detect_dca") # Enable DCA detection
        } else {
            hc2 <- NULL
        }
        # hc <- NULL
    }

    gset_genome(shiny_mct@genome)
    min_val <- round(gsummary(shiny_mct@tracks[1], intervals = gintervals.all()[1, ])[[6]] * 10 / 2)
    gset_genome(shiny_mct2@genome)
    min_val2 <- round(gsummary(shiny_mct2@tracks[1], intervals = gintervals.all()[1, ])[[6]] * 10 / 2)
    gset_genome(shiny_mct@genome)

    observe({
        updateNumericInput(session = session, inputId = "min_color", value = min_val)
        updateNumericInput(session = session, inputId = "max_color", value = min_val * 4)
        updateNumericInput(session = session, inputId = "min_color2", value = min_val2)
        updateNumericInput(session = session, inputId = "max_color2", value = min_val2 * 4)
        updateNumericInput(session = session, inputId = "n_smooth", value = shiny_mct@resolution * 10)
    })

    observe({
        shinyjs::toggle(id = "min_color", condition = input$show_controls)
        shinyjs::toggle(id = "max_color", condition = input$show_controls)
        shinyjs::toggle(id = "min_color2", condition = input$show_controls)
        shinyjs::toggle(id = "max_color2", condition = input$show_controls)
        shinyjs::toggle(id = "n_smooth", condition = input$show_controls)
    })

    promoters <- misha.ext::get_promoters(upstream = 5e4, downstream = 5e4) %>%
        mutate(
            coords = glue("{chrom}:{start}-{end}"),
            label = glue("{geneSymbol} ({coords})")
        )

    observe({
        shinyWidgets::updateVirtualSelect(
            "genes",
            choices = shinyWidgets::prepare_choices(promoters, label, coords),
            selected = "chrX:60843429-60943430"
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

    observeEvent(input$region_plot_cor_click, {
        mat_x <- input$region_plot_cor_click$x
        req(mat_x)
        mat_y <- input$region_plot_cor_click$y
        req(mat_y)
        m <- mat()
        req(m)
        mat_idxs <- get_heatmap_idx(mat_x, mat_y, m)

        gene <- promoters[promoters$coords == input$genes, ]$geneSymbol
        req(gene %in% rownames(shiny_mct@rna_egc))

        mcs <- intersect(colnames(shiny_mct@rna_egc), colnames(m))

        data <- data.frame(atac = m[mat_idxs[1], mcs], rna = shiny_mct@rna_egc[gene, mcs], metacell = mcs) %>%
            mutate(rna = log2(rna + 1e-5)) %>%
            left_join(shiny_mct@metadata, by = "metacell") %>%
            mutate(coords = rownames(m)[mat_idxs[1]])

        scatter_data(data)
    })
    
    observeEvent(input$region_plot_cor2_click, {
        mat_x <- input$region_plot_cor2_click$x
        req(mat_x)
        mat_y <- input$region_plot_cor2_click$y
        req(mat_y)
        m <- mat2()
        req(m)
        req(intervals_comparison())
        flipped=is_comparison_flipped(intervals_comparison())
        if(flipped){
            m = m[dim(m)[1]:1,]
        }
        mat_idxs <- get_heatmap_idx(mat_x, mat_y, m)

        gene <- toupper(promoters[promoters$coords == input$genes, ]$geneSymbol)
        req(gene %in% rownames(shiny_mct2@rna_egc))

        mcs <- intersect(colnames(shiny_mct2@rna_egc), colnames(m))

        data <- data.frame(atac = m[mat_idxs[1], mcs], rna = shiny_mct2@rna_egc[gene, mcs], metacell = mcs) %>%
            mutate(rna = log2(rna + 1e-5)) %>%
            left_join(shiny_mct2@metadata, by = "metacell") %>%
            mutate(coords = rownames(m)[mat_idxs[1]])

        scatter_data(data)
    })

    output$pheatmap_plot <- renderPlot({
        req(intervals())
        req(intervals2())
        req(intervals_comparison())
        req(dim(intervals_comparison()[!is.na(intervals_comparison()$start1),])[1]>1)
        req(shiny_mct)
        req(shiny_mct2)
        req(has_rna(shiny_mct2))
        req(comb_mat())
        req(comb_mat_txt())
        paletteLength=1000

        myColor <- c(
        #         colorRampPalette(c("black","grey","white","white"))(10),
                colorRampPalette(c("blue","white", "red"))(paletteLength)
                    )
        myBreaks <- c(
        #     seq(-2.0, -1.001, length.out=10),
            seq(-1.0, 0.0, length.out=ceiling(paletteLength/2) + 1),
            seq(0.01, 1.0, length.out=floor(paletteLength/2)))

        pheatmap::pheatmap(
            comb_mat(),
            color = myColor, breaks = myBreaks,
            display_numbers = comb_mat_txt(),
            cluster_rows = F, cluster_cols = F,
            labels_row='', labels_col=''
        )
    })

    output$scatter_plot <- renderPlot({
        data <- scatter_data()
        req(data)

        data %>%
            ggplot(aes(x = atac, y = rna, color = color)) +
            geom_point(size = 2) +
            scale_color_identity() +
            theme_classic() +
            xlab("ATAC") +
            ylab("RNA") +
            ggtitle(glue("Correlation: {round(cor(data$atac, data$rna), digits = 2)} Coverage:{sum(data$atac)}"), subtitle = data$coords[1])
    })

    observe({
        req(intervals())
        req(input$n_smooth)
        req(input$n_smooth >= 1)
        n_smooth <- max(round(input$n_smooth / shiny_mct@resolution), 1)
        raw_mat <- get_raw_mat(
            shiny_mct, intervals(),
            downsample = TRUE,
            downsample_n = 20000000,
            detect_dca = input$detect_dca %||% FALSE,
            gene_annot = TRUE,
            hc = shiny_hc,
            peak_lf_thresh1 = input$dca_peak_lf_thresh,
            trough_lf_thresh1 = input$dca_trough_lf_thresh,
            sz_frac_for_peak = input$dca_sz_frac_for_peak,
            color_breaks = color_breaks,
            n_smooth = n_smooth,
            plot_x_axis_ticks = FALSE
        )
        mat(raw_mat)
    })

    observe({
        req(intervals2())
        req(input$n_smooth)
        req(input$n_smooth >= 1)
        n_smooth <- max(round(input$n_smooth / shiny_mct2@resolution), 1)
        raw_mat2 <- get_raw_mat(
            shiny_mct2, intervals2(),
            downsample = TRUE,
            downsample_n = 2500000,
            detect_dca = input$detect_dca %||% FALSE,
            gene_annot = TRUE,
            hc = shiny_hc2,
            peak_lf_thresh1 = input$dca_peak_lf_thresh,
            trough_lf_thresh1 = input$dca_trough_lf_thresh,
            sz_frac_for_peak = input$dca_sz_frac_for_peak,
            color_breaks = color_breaks2,
            n_smooth = n_smooth,
            plot_x_axis_ticks = FALSE
        )
        mat2(raw_mat2)
    })

    observe({
        req(intervals())
        req(intervals2())
        req(intervals_comparison())
        req(shiny_mct)
        req(shiny_mct2)
        req(shiny_all_peaks)
        req(shiny_all_peaks2)
        gset_genome(shiny_mct@genome)
        scope_peaks = gintervals.neighbors(shiny_all_peaks, intervals(), mindist = 0, maxdist = 0)[,1:length(colnames(shiny_all_peaks))]
        gset_genome(shiny_mct2@genome)
        scope_peaks2 = gintervals.neighbors(gintervals.force_range(as.data.frame(shiny_all_peaks2)), intervals2(), mindist = 0, maxdist = 0)[,1:length(colnames(shiny_all_peaks2))]
        gset_genome(shiny_mct@genome)
        req(!is.null(scope_peaks))
        req(!is.null(scope_peaks2))
        comb_mat_tmp = get_comb_mat(
            promoters[promoters$coords == input$genes, ]$geneSymbol,
            shiny_mct, 
            shiny_mct2, 
            shiny_all_peaks %>% misha.ext::intervs_to_mat(remove_intervalID = TRUE),
            shiny_all_peaks2 %>% misha.ext::intervs_to_mat(remove_intervalID = TRUE), 
            scope_peaks, 
            scope_peaks2,
            flipped=is_comparison_flipped(intervals_comparison()))
        comb_mat(comb_mat_tmp)
    })


    observe({
        req(intervals())
        req(intervals2())
        req(intervals_comparison())
        req(shiny_mct)
        req(shiny_mct2)
        req(shiny_all_peaks)
        req(shiny_all_peaks2)
        req(comb_mat())
        gset_genome(shiny_mct@genome)
        scope_peaks = gintervals.neighbors(shiny_all_peaks, intervals(), mindist = 0, maxdist = 0)[,1:length(colnames(shiny_all_peaks))]
        gset_genome(shiny_mct2@genome)
        scope_peaks2 = gintervals.neighbors(gintervals.force_range(as.data.frame(shiny_all_peaks2)), intervals2(), mindist = 0, maxdist = 0)[,1:length(colnames(shiny_all_peaks2))]
        gset_genome(shiny_mct@genome)
        dot_mat = get_dot_mat(shiny_mct, shiny_mct2, scope_peaks, scope_peaks2)
        if(is_comparison_flipped(intervals_comparison())){
            dot_mat = dot_mat[dim(dot_mat)[1]:1,]
        }
        txt = comb_mat()
        txt = txt-txt
        txt[txt==0] = ""
        txt[which(t(dot_mat) < -1.01, arr.ind = TRUE)+1] = "X"
        txt[cbind(
                (which(dot_mat < -1.01, arr.ind = TRUE)+1)[,2],
                (which(dot_mat < -1.01, arr.ind = TRUE)+1)[,2] + dim(scope_peaks)[1]
                )] = "O"

        txt[cbind(
                (which(dot_mat < -1.01, arr.ind = TRUE)+1)[,1] + dim(scope_peaks2)[1],
                (which(dot_mat < -1.01, arr.ind = TRUE)+1)[,1]
                )] = "M"
        comb_mat_txt(txt)
    })

    output$region_plot <- renderPlot(
        {
            req(intervals())
            req(input$dca_peak_lf_thresh)
            req(input$dca_trough_lf_thresh)
            req(input$dca_sz_frac_for_peak)
            req(input$min_color)
            req(input$max_color)
            req(input$min_color < input$max_color)
            req(input$n_smooth)
            req(input$n_smooth >= 1)
            req(mat())
            color_breaks <- c(0, seq(
                input$min_color,
                input$max_color,
                length.out = 4
            ))

            n_smooth <- max(round(input$n_smooth / shiny_mct@resolution), 1)

            mct_plot_region(
                shiny_mct, intervals(),
                downsample = TRUE,
                downsample_n = 20000000,
                detect_dca = input$detect_dca %||% FALSE,
                gene_annot = TRUE,
                hc = hc,
                peak_lf_thresh1 = input$dca_peak_lf_thresh,
                trough_lf_thresh1 = input$dca_trough_lf_thresh,
                sz_frac_for_peak = input$dca_sz_frac_for_peak,
                color_breaks = color_breaks,
                n_smooth = n_smooth,
                plot_x_axis_ticks = FALSE
            )
        },
        res = 96
    ) %>%
        bindCache(
            intervals(),
            input$detect_dca,
            input$dca_peak_lf_thresh,
            input$dca_trough_lf_thresh,
            input$dca_sz_frac_for_peak,
            input$min_color,
            input$max_color,
            input$n_smooth
        )
    output$region_plot_cor <- renderPlot(
        {
            req(intervals())
            req(input$dca_peak_lf_thresh)
            req(input$dca_trough_lf_thresh)
            req(input$dca_sz_frac_for_peak)
            req(has_rna(shiny_mct))
            req(input$min_color)
            req(input$max_color)
            req(input$min_color < input$max_color)
            req(input$n_smooth)
            req(input$n_smooth >= 1)
            color_breaks <- c(0, seq(
                input$min_color,
                input$max_color,
                length.out = 4
            ))

            n_smooth <- max(round(input$n_smooth / shiny_mct@resolution), 1)

            mct_plot_region(
                shiny_mct, intervals(),
                downsample = TRUE,
                downsample_n = 20000000,
                detect_dca = input$detect_dca %||% FALSE,
                gene_annot = FALSE,
                hc = NULL,
                peak_lf_thresh1 = input$dca_peak_lf_thresh,
                trough_lf_thresh1 = input$dca_trough_lf_thresh,
                sz_frac_for_peak = input$dca_sz_frac_for_peak,
                color_breaks = color_breaks,
                n_smooth = n_smooth,
                plot_x_axis_ticks = FALSE,
                roll_mean = TRUE,
                genes_correlations = promoters[promoters$coords == input$genes, ]$geneSymbol,
                cor_color_breaks = c(-1.0, -0.5, 0, 0.5, 1.0),
                cor_colors = c("blue", "white", "white", "white", "red")
            )
        },
        res = 96
    ) %>%
        bindCache(
            intervals(),
            input$detect_dca,
            input$dca_peak_lf_thresh,
            input$dca_trough_lf_thresh,
            input$dca_sz_frac_for_peak,
            input$min_color,
            input$max_color,
            input$n_smooth
        )

    observe({
        req(intervals())
        if (!is.null(chain_1_to_2)) {
            # intervs2 <- translate_and_center(intervals(), chain_1_to_2, chainscore = 78939364) %>%
            intervs2 <- translate_and_center(intervals(), chain_1_to_2, chainscore = NULL)
            req(intervs2)
            intervs2 <- intervs2 %>%
                select(chrom, start, end) %>%
                as.data.frame()
            intervals2(intervs2)
        } else {
            intervals2(intervals())
        }
    })

    observe({
        req(intervals())
        req(intervals2())
        req(chain_1_to_2)
        req(chain_2_to_1)
        i12 <- compute_intervals_comparison(intervals(), intervals2(), chain_1_to_2, chainscore = NULL, grid_resolution = 100)
        intervals_comparison(i12)
    })

    output$comparison_plot <- renderPlot(
        {
            req(intervals())
            req(intervals2())
            req(intervals_comparison())
            req(dim(intervals_comparison()[!is.na(intervals_comparison()$start1),])[1]>1)
            req(shiny_mct2)
            req(chain_1_to_2)
            req(chain_2_to_1)
            gset_genome(mct@genome)
            layout(matrix(2:1, nrow = 1), w = c(1, 20))
            par(mar = c(0, 0, 0, 2))
            par(xaxs = "i")
            plot_intervals_comparison(
                intervals_comparison(),
                annotations = compute_annotation_comparison(shiny_annotations1, shiny_annotations2, intervals(), intervals2(), chain_1_to_2, chain_2_to_1),
                grid_resolution = 100
            )
        },
        res = 96
    )

    output$region_plot_cor2 <- renderPlot(
        {
            req(intervals())
            req(intervals2())
            req(intervals_comparison())
            req(dim(intervals_comparison()[!is.na(intervals_comparison()$start1),])[1]>1)
            req(shiny_mct2)
            req(has_rna(shiny_mct2))
            req(input$dca_peak_lf_thresh)
            req(input$dca_trough_lf_thresh)
            req(input$dca_sz_frac_for_peak)
            req(input$min_color2)
            req(input$max_color2)
            req(input$min_color2 < input$max_color2)
            req(length(intersect(toupper(promoters[promoters$coords == input$genes, ]$geneSymbol), rownames(shiny_mct2@rna_egc))) > 0)
            req(input$n_smooth)
            req(input$n_smooth >= 1)
            color_breaks2 <- c(0, seq(
                input$min_color2,
                input$max_color2,
                length.out = 4
            ))

            n_smooth <- max(round(input$n_smooth / shiny_mct2@resolution), 1)

            mct_plot_region(
                shiny_mct2, intervals2(),
                downsample = TRUE,
                downsample_n = 2500000,
                detect_dca = input$detect_dca %||% FALSE,
                gene_annot = FALSE,
                hc = NULL,
                peak_lf_thresh1 = input$dca_peak_lf_thresh,
                trough_lf_thresh1 = input$dca_trough_lf_thresh,
                sz_frac_for_peak = input$dca_sz_frac_for_peak,
                color_breaks = color_breaks2,
                n_smooth = n_smooth,
                plot_x_axis_ticks = FALSE,
                genes_correlations = toupper(promoters[promoters$coords == input$genes, ]$geneSymbol),
                roll_mean = TRUE,
                flip = is_comparison_flipped(intervals_comparison()),
                cor_color_breaks = c(-1.0, -0.5, 0, 0.5, 1.0),
                cor_colors = c("blue", "white", "white", "white", "red")
            )

            gset_genome(mct@genome)
        },
        res = 96
    ) %>%
        bindCache(
            intervals2(),
            input$detect_dca,
            input$dca_peak_lf_thresh,
            input$dca_trough_lf_thresh,
            input$dca_sz_frac_for_peak,
            input$min_color2,
            input$max_color2,
            input$n_smooth
        )

    output$region_plot2 <- renderPlot(
        {
            req(intervals())
            req(intervals2())
            req(intervals_comparison())
            req(dim(intervals_comparison()[!is.na(intervals_comparison()$start1),])[1]>1)
            req(shiny_mct2)
            req(input$dca_peak_lf_thresh)
            req(input$dca_trough_lf_thresh)
            req(input$dca_sz_frac_for_peak)
            req(input$min_color2)
            req(input$max_color2)
            req(input$min_color2 < input$max_color2)
            req(input$n_smooth)
            req(input$n_smooth >= 1)
            color_breaks2 <- c(0, seq(
                input$min_color2,
                input$max_color2,
                length.out = 4
            ))

            n_smooth <- max(round(input$n_smooth / shiny_mct2@resolution), 1)

            mct_plot_region(
                shiny_mct2, intervals2(),
                downsample = TRUE,
                downsample_n = 2500000,
                detect_dca = input$detect_dca %||% FALSE,
                gene_annot = TRUE,
                hc = shiny_hc2,
                peak_lf_thresh1 = input$dca_peak_lf_thresh,
                trough_lf_thresh1 = input$dca_trough_lf_thresh,
                sz_frac_for_peak = input$dca_sz_frac_for_peak,
                color_breaks = color_breaks2,
                n_smooth = n_smooth,
                gene_annot_pos = "bottom",
                flip = is_comparison_flipped(intervals_comparison()),
                plot_x_axis_ticks = FALSE
            )
            gset_genome(mct@genome)
        },
        res = 96
    ) %>%
        bindCache(
            intervals2(),
            input$detect_dca,
            input$dca_peak_lf_thresh,
            input$dca_trough_lf_thresh,
            input$dca_sz_frac_for_peak,
            input$min_color2,
            input$max_color2,
            input$n_smooth
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

    output$hover_info <- renderUI({
        hover <- input$region_hover
        req(hover)
        raw_mat <- mat()
        req(raw_mat)
        mat_idxs <- get_heatmap_idx(hover$x, hover$y, raw_mat)

        mc <- colnames(raw_mat)[mat_idxs[2]]
        req(mc)
        mc_metadata <- mct@metadata %>%
            filter(metacell == mc) %>%
            slice(1)
        mc_color <- mc_metadata$color
        mc_cell_type <- mc_metadata$cell_type
        coords <- rownames(raw_mat)[mat_idxs[1]]

        # taken from https://gitlab.com/-/snippets/16220
        left_px <- hover$coords_css$x
        top_px <- hover$coords_css$y
        style <- glue(
            "background-color: {grDevices::rgb(t(grDevices::col2rgb(mc_color) / 255))}; position:absolute; z-index:100; left: {left_px + 2}px; top: {top_px + 200}px;"
        )
        tooltip <- paste(
            glue("Metacell: {mc}"),
            glue("Cell type: {mc_cell_type}"),
            glue("Coordinates: {coords}"),
            sep = "<br/>"
        )

        wellPanel(
            style = style,
            p(HTML(tooltip))
        )
    })

    observe({
        shinyjs::toggle(id = "dca_peak_lf_thresh", condition = input$detect_dca && input$show_controls)
        shinyjs::toggle(id = "dca_trough_lf_thresh", condition = input$detect_dca && input$show_controls)
        shinyjs::toggle(id = "dca_sz_frac_for_peak", condition = input$detect_dca && input$show_controls)
    })
}

parse_coordinate_text <- function(text) {
    text <- gsub(",", "", text) # remove commas
    text <- gsub("\t", " ", text) # transform tabs to spaces
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
#' @param mct2 MCT object from a different genome (optional)
#' @param hc an hclust object with clustering of the metacells (optional, see \code{mct_plot_region}
#' @param hc2 an hclust object with clustering of the metacells for the second genome (optional, see \code{mct_plot_region}
#' @param chain filename of a chain translating from the genome of mct to the genome of mct2, or an rtracklayer chain object (optional)
#' @param chain2 filename of a chain translating from the genome of mct2 to the genome of mct, or an rtracklayer chain object (optional)
#' @param annotations1 an intervals object with annotations for the first genome (optional)
#' @param annotations2 an intervals object with annotations for the second genome (optional)
#'
#' @examples
#' \dontrun{
#' run_app(mct)
#' }
#'
#' @export
run_app <- function(mct,
                    hc = NULL,
                    mct2 = NULL,
                    hc2 = NULL,
                    port = NULL,
                    host = NULL,
                    launch.browser = FALSE,
                    chain = NULL,
                    chain2 = NULL,
                    annotations1 = NULL,
                    annotations2 = NULL,
                    all_peaks = NULL,
                    all_peaks2 = NULL) {
    library(misha)
    library(shiny)
    opt <- options(gmultitasking = FALSE, shiny.usecairo = TRUE)
    shiny_mct <<- mct
    shiny_hc <<- hc
    gset_genome(mct@genome)

    shiny_mct2 <<- mct2
    shiny_hc2 <<- hc2

    chain_1_to_2 <<- NULL
    if (!is.null(chain)) {
        if (is.null(mct2)) {
            cli_abort("Please provide mct2 when providing chain")
        }
        chain_1_to_2 <<- load_chain(chain, genome1 = mct@genome, genome2 = mct2@genome)
    }

    chain_2_to_1 <<- NULL
    if (!is.null(chain2)) {
        if (is.null(mct2)) {
            cli_abort("Please provide mct2 when providing chain")
        }
        chain_2_to_1 <<- load_chain(chain2, genome1 = mct2@genome, genome2 = mct@genome)
    }

    shiny_annotations1 <<- annotations1
    shiny_annotations2 <<- annotations2

    shiny_all_peaks <<- all_peaks
    shiny_all_peaks2 <<- all_peaks2

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
#' @param hc an hclust object with clustering of the metacells (optional, see \code{mct_plot_region})
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
create_bundle <- function(mct, path, hc = NULL, overwrite = FALSE, self_contained = FALSE, branch = "master", restart = overwrite, permissions = NULL, ...) {
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
    if (!is.null(hc)) {
        readr::write_rds(hc, fs::path(bundle_dir, "hc.rds"))
    }
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
