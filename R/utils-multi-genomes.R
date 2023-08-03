flipStrandTricky <- function(strand, flip) {
    strandCodes <- c("+" = 1L, "-" = -1L, "*" = 0L)
    strandInt <- strandCodes[as.vector(strand)]
    flipped <- ifelse(flip, strandInt * -1L, strandInt) + 2L
    strandRevCodes <- factor(c("-", "*", "+"), levels(rtracklayer::strand()))
    strandRevCodes[as.vector(flipped)]
}

custom_lift <- function(x, chain) {
    liftOverSpace <- function(gr, chain, subind) {
        r <- rtracklayer::ranges(gr)
        ol <- IRanges::findOverlaps(r, rtracklayer::ranges(chain))
        shits <- S4Vectors::subjectHits(ol)
        r <- IRanges::overlapsRanges(r, rtracklayer::ranges(chain), ol)
        rev <- as.vector(rtracklayer::reversed(chain)[shits])
        starts <- ifelse(rev,
            start(IRanges::reflect(r, rtracklayer::ranges(chain)[shits])),
            start(r)
        )
        strand <- flipStrandTricky(rtracklayer::strand(gr)[S4Vectors::queryHits(ol)], rev)
        r <- IRanges::IRanges(starts, width = rtracklayer::width(r))
        offsets <- rtracklayer::offset(chain)[shits]
        spaces <- rtracklayer::space(chain)[shits]
        chainscore <- rtracklayer::score(chain)[shits]
        ind[[as.character(GenomeInfoDb::seqnames(gr)[1])]] <<- subind[S4Vectors::queryHits(ol)]
        GenomicRanges::GRanges(spaces,
            IRanges::IRanges(start(r) - offsets, end(r) - offsets),
            strand = strand,
            rtracklayer::values(gr)[S4Vectors::queryHits(ol), , drop = FALSE],
            chainscore
        )
    }
    rl <- split(x, GenomeInfoDb::seqnames(x), drop = TRUE)
    unchainedNames <- setdiff(names(rl), names(chain))
    if (length(unchainedNames)) {
        message(
            "Discarding unchained sequences: ",
            paste(unchainedNames, collapse = ", ")
        )
    }
    sharedNames <- intersect(names(rl), names(chain))
    ind <- split(
        seq_len(length(x)),
        as.vector(GenomeInfoDb::seqnames(x))
    )[sharedNames]
    liftedList <- mapply(liftOverSpace, rl[sharedNames],
        chain[sharedNames], ind,
        SIMPLIFY = FALSE
    )
    lifted <- unlist(GenomicRanges::GRangesList(liftedList), use.names = FALSE)

    f <- structure(as.integer(unlist(ind, use.names = FALSE)),
        levels = seq_len(length(x)), class = "factor"
    )

    GenomicRanges::GRangesList(setNames(split(lifted, f), names(x)))
}


custom_lift2 <- function(x, chain, maxgap = 1000) {
    test2 <- custom_lift(x, chain)
    test3 <- as.data.frame(test2) %>%
        group_by(across(c(-start, -end, -width))) %>%
        summarise(start = min(start), end = max(end), widthsum = sum(width)) %>%
        ungroup()
    test3$alt_start <- test3$start - 1
    test3$match <- test3$widthsum / (test3$end - test3$alt_start)
    atest2 <- as.data.frame(test2) %>%
        left_join(
            ungroup(test3) %>%
                select(seqnames, strand, chainscore, lstart = alt_start, lend = end, widthsum, match) %>%
                mutate(lrow_id = row_number()),
            by = join_by(seqnames, strand, chainscore)
        )
    btest2 <- atest2 %>%
        group_by(
            group, group_name, seqnames, strand, chainscore, lrow_id
        ) %>%
        arrange(
            start,
            .by_group = TRUE
        ) %>%
        mutate(
            gap = start - lag(end)
        ) %>%
        mutate(
            nbreaks = cumsum(ifelse(is.na(gap > maxgap), 0, gap > maxgap))
        ) %>%
        group_by(
            nbreaks,
            .add = TRUE
        ) %>%
        summarise(
            start = min(start), end = max(end), widthsum = sum(width)
        ) %>%
        ungroup()
    btest2$alt_start <- btest2$start - 1
    btest2$match <- btest2$widthsum / (btest2$end - btest2$alt_start)
    return(list(full = test2, gapped = btest2, cont = test3))
}

translate_intervals <- function(intervals, chain, chainscore = NULL) {
    r <- regioneR::toGRanges(intervals %>% rowid_to_column("row_ID") %>% select(chrom, start, end, row_ID))
    lifted_list <- custom_lift2(r, chain)

    if (is.null(chainscore)) {
        chainscore <- max(lifted_list$cont$chainscore)
    }

    lifted_list$cont %>%
        filter(chainscore == !!chainscore) %>%
        select(chrom = seqnames, start, end, row_ID) %>%
        as.data.frame() %>%
        select(chrom, start, end, everything())
}

translate_and_center <- function(intervals, chain, chainscore = NULL) {
    r <- regioneR::toGRanges(intervals)
    lifted_list <- custom_lift2(r, chain)

    if (!is.null(chainscore)) {
        intervs2 <- lifted_list$cont %>%
            filter(chainscore == !!chainscore) %>%
            select(chrom = seqnames, start, end) %>%
            as.data.frame()
    } else {
        intervs2 <- lifted_list$cont %>%
            arrange(desc(chainscore)) %>%
            ungroup() %>%
            slice(1) %>%
            select(chrom = seqnames, start, end) %>%
            as.data.frame()
    }

    expand <- round(abs(intervals$end - intervals$start) / 2)
    intervs2 <- misha.ext::gintervals.centers(intervs2) %>%
        mutate(start = start - expand, end = end + expand)

    return(intervs2)
}


get_comparison <- function(elems, chain, selected_chain_chainscore = NULL, color = "grey") {
    elems_gr <- regioneR::toGRanges(elems %>% rowid_to_column("row_ID") %>% select(chrom, start, end, row_ID))
    lifted_elems <- custom_lift2(elems_gr, chain)
    if (is.null(selected_chain_chainscore)) {
        selected_chain_chainscore <- max(lifted_elems$cont$chainscore)
    }
    proj_elems <- as.comparison(elems %>% rowid_to_column("row_ID") %>% merge(
        lifted_elems$cont %>% filter(chainscore == selected_chain_chainscore),
        by = "row_ID",
    ) %>% mutate(direction = 1) %>% select(start1 = start.x, end1 = end.x, start2 = start.y, end2 = end.y, direction))
    proj_elems$color <- color
    proj_elems
}

orth_one_per_element <- function(comparison_obj, width = 10) {
    comparison_obj$mid1 <- comparison_obj$start1 + as.integer((comparison_obj$end1 - comparison_obj$start1) / 2)
    comparison_obj$mid2 <- comparison_obj$start2 + as.integer((comparison_obj$end2 - comparison_obj$start2) / 2)
    comparison_obj$start1 <- comparison_obj$mid1
    comparison_obj$end1 <- comparison_obj$mid1 + width
    comparison_obj$start2 <- comparison_obj$mid2
    comparison_obj$end2 <- comparison_obj$mid2 + width
    comparison_obj[, 1:6]
}

mct_plot_comparison <- function(mct, intervals, intervals2, chain, chain2, selected_chain_chainscore = NULL, annot1 = NULL, annot2 = NULL) {
    mm_elems <- intervals

    if (!("strand" %in% colnames(mm_elems))) {
        mm_elems$strand <- "1"
    }

    lifted_elems <- custom_lift2(
        regioneR::toGRanges(mm_elems %>% rowid_to_column("row_ID") %>% select(chrom, start, end, row_ID)),
        chain
    )

    if (is.null(selected_chain_chainscore)) {
        selected_chain_chainscore <- max(lifted_elems$cont$chainscore)
    }

    oc_elems <- ungroup(lifted_elems$cont %>% filter(chainscore == selected_chain_chainscore)) %>% select(chrom = seqnames, start = start, end = end, row_ID = row_ID)
    oc_elems$strand <- "1"

    # mm_oc_comp <- get_comparison(mm_elems, chain, selected_chain_chainscore)

    mm_seg <- dna_seg(mm_elems %>% rowid_to_column("row_ID") %>% mutate(name = as.character(row_ID)) %>% select(name, start, end, strand))
    oc_seg <- dna_seg(oc_elems %>% mutate(name = as.character(row_ID)) %>% select(name, start, end, strand))

    dna_segs <- list(
        mm = mm_seg,
        oc = oc_seg
    )

    xlims <- list(c(intervals$start, intervals$end), c(intervals2$start, intervals2$end))

    # xlims <- list()
    # for (i in 1:length(dna_segs)) {
    #     rng <- range(dna_segs[[i]])
    #     xlims[[i]] <- c(rng[1] - 0.02 * diff(rng), rng[2] + 0.02 * diff(rng))
    # }


    grid_resolution <- (intervals$end - intervals$start) / 50
    ruler <- get_comparison(giterator.intervals(intervals = intervals, iterator = grid_resolution), chain, selected_chain_chainscore = selected_chain_chainscore)

    flip <- cor(
        orth_one_per_element(ruler)$start1,
        orth_one_per_element(ruler)$start2
    ) < 0

    if (flip) {
        xlims[[2]] <- c(xlims[[2]][2], xlims[[2]][1])
    }

    seg_plots <- list()
    if (!is.null(annot1)) {
        annot1 <- annot1 %>%
            gintervals.neighbors1(intervals) %>%
            filter(dist == 0)
        annot1_colors <- chameleon::distinct_colors(nrow(annot1))$name

        lifted_annot1 <- translate_intervals(as.data.frame(annot1), chain = chain, chainscore = selected_chain_chainscore) %>%
            arrange(row_ID)
    }

    if (!is.null(annot2)) {
        annot2 <- annot2 %>%
            filter(
                chrom == intervals2$chrom,
                start >= intervals2$start,
                end <= intervals2$end
            )
        annot2_colors <- chameleon::distinct_colors(nrow(annot1) + nrow(annot2))$name[(nrow(annot1) + 1):(nrow(annot1) + nrow(annot2))]
        lifted_annot2 <- translate_intervals(as.data.frame(annot2), chain = chain2, chainscore = selected_chain_chainscore) %>%
            arrange(row_ID)
    }


    if (!is.null(annot1)) {
        if (!is.null(annot2)) {
            seg_plots[[1]] <- seg_plot(
                func = pointsGrob,
                args = list(
                    x = c(annot1$start, lifted_annot2$start),
                    y = c(rnorm(nrow(annot1), sd = 0.1) + 1, rnorm(nrow(lifted_annot2), sd = 0.1) + 5),
                    default.units = "native", pch = 3,
                    gp = gpar(col = c(annot1_colors, annot2_colors), cex = 0.5)
                )
            )
        } else {
            seg_plots <- list(
                seg_plots[[1]],
                NULL
            )
        }
    }

    if (!is.null(annot2)) {
        if (!is.null(annot1)) {
            seg_plots[[2]] <- seg_plot(
                func = pointsGrob,
                args = list(
                    x = c(annot2$start, lifted_annot1$start),
                    y = c(rnorm(nrow(annot2), sd = 0.1) + 1, rnorm(nrow(lifted_annot1), sd = 0.1) + 5),
                    default.units = "native", pch = 3,
                    gp = gpar(col = c(annot2_colors, annot1_colors), cex = 0.5)
                )
            )
        } else {
            seg_plots <- list(
                NULL,
                seg_plots[[2]]
            )
        }
    }


    vp <- grid::plotViewport(margins = c(0, 2, 0, 0))


    grid::pushViewport(vp)
    genoPlotR::plot_gene_map(
        offsets = c(0, 0),
        dna_segs = dna_segs,
        # annotations=list(
        #     annot,
        #     NULL
        # ),
        seg_plots = seg_plots,
        seg_plot_height = 2,
        seg_plot_y_axis = NULL,
        comparisons = list(
            # orth_one_per_element(mm_oc_comp)
            orth_one_per_element(ruler)
        ),
        dna_seg_scale = c(
            TRUE,
            TRUE
        ),
        xlims = xlims,
        plot_new = FALSE
    )
}

load_chain <- function(chain) {
    if (is.character(chain)) {
        cli::cli_alert_info("Loading chain {.file {chain}}")
        chain <<- rtracklayer::import.chain(chain)
        cli::cli_alert_success("Loaded chain {.file {chain}} successfully")
    }

    return(chain)
}

plot_intervals_comparison <- function(intervals, intervals2, chain, chain2, chainscore = NULL, annot1 = NULL, annot2 = NULL) {
    grid_resolution <- round((intervals$end - intervals$start) / 50)

    # create a grid iterator for the first intervals set
    i1 <- tibble(
        chrom = intervals$chrom,
        end = seq(intervals$start, intervals$end, grid_resolution),
        start = end - grid_resolution
    ) %>%
        mutate(end = pmin(end, intervals$end)) %>%
        as.data.frame()


    i12 <- translate_intervals(i1, chain, chainscore = chainscore)

    i12 <- i1 %>%
        rowid_to_column("row_ID") %>%
        left_join(i12 %>% rename(chrom1 = chrom, start1 = start, end1 = end), by = join_by(row_ID)) %>%
        as_tibble()

    # translate start1 and start2 to relative coordinates (from 0 to 1) according to intervals1 and intervals2
    i12 <- i12 %>%
        mutate(
            x1 = (start - intervals$start) / (intervals$end - intervals$start),
            x2 = (start1 - intervals2$start) / (intervals2$end - intervals2$start),
            x2 = ifelse(is.na(chrom1) | as.character(chrom1) != as.character(intervals2$chrom), NA, x2)
        )

    # Set the plot parameters
    plot(1, 1, xlim = c(0, 1), ylim = c(0, 1), type = "n", ann = FALSE, axes = FALSE)

    segments(x0 = i12$x1, y0 = 1, x1 = i12$x1, y1 = 0.99, lwd = 1)
    text(x = i12$x1, y = 0.96, labels = round(i12$start / 1e+6, 3), srt = 90, adj = c(1, 0.5), xpd = TRUE, cex = 0.5)
    segments(x0 = i12$x1, y0 = 0.75, x1 = i12$x1, y1 = 0.7, lwd = 1)
    segments(x0 = i12$x1, y0 = 0.7, x1 = i12$x2, y1 = 0.27, lwd = 1)
    segments(x0 = i12$x2, y0 = 0.27, x1 = i12$x2, y1 = 0.25, lwd = 1)
    text(x = i12$x2[i12$x2 <= 1 & i12$x2 >= 0], y = 0.2, labels = round(i12$start1[i12$x2 <= 1 & i12$x2 >= 0] / 1e+6, 3), srt = 90, adj = c(1, 0.5), xpd = TRUE, cex = 0.5)
    segments(x0 = i12$x2, y0 = 0, x1 = i12$x2, y1 = 0.01, lwd = 1)
}
