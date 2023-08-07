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


orth_one_per_element <- function(comparison_obj, width = 10) {
    comparison_obj$mid1 <- comparison_obj$start1 + as.integer((comparison_obj$end1 - comparison_obj$start1) / 2)
    comparison_obj$mid2 <- comparison_obj$start2 + as.integer((comparison_obj$end2 - comparison_obj$start2) / 2)
    comparison_obj$start1 <- comparison_obj$mid1
    comparison_obj$end1 <- comparison_obj$mid1 + width
    comparison_obj$start2 <- comparison_obj$mid2
    comparison_obj$end2 <- comparison_obj$mid2 + width
    comparison_obj[, 1:6]
}

load_chain <- function(chain) {
    if (is.character(chain)) {
        cli::cli_alert_info("Loading chain {.file {chain}}")
        chain <<- rtracklayer::import.chain(chain)
        cli::cli_alert_success("Loaded chain {.file {chain}} successfully")
    }

    return(chain)
}

compute_annotation_comparison <- function(annot1, annot2, intervals, intervals2, chain_1_to_2, chain_2_to_1, selected_chain_chainscore=NULL){
    annot1 <- annot1 %>%
        gintervals.neighbors1(intervals) %>%
        filter(dist == 0)
    annot1_colors <- chameleon::distinct_colors(nrow(annot1))$name

    lifted_annot1 <- translate_intervals(as.data.frame(annot1), chain = chain_1_to_2, chainscore = selected_chain_chainscore) %>%
        arrange(row_ID)

    annot2 <- annot2 %>%
        filter(
            chrom == intervals2$chrom,
            start >= intervals2$start,
            end <= intervals2$end
        )
    annot2_colors <- chameleon::distinct_colors(nrow(annot1) + nrow(annot2))$name[(nrow(annot1) + 1):(nrow(annot1) + nrow(annot2))]
    lifted_annot2 <- translate_intervals(as.data.frame(annot2), chain = chain_2_to_1, chainscore = selected_chain_chainscore) %>%
        arrange(row_ID)
    list(
         annot1=annot1 %>% mutate(x1 = (start - intervals$start) / (intervals$end - intervals$start)), 
         lifted_annot1=lifted_annot1 %>% mutate(x2 = (start - intervals2$start) / (intervals2$end - intervals2$start)), 
         annot1_colors=annot1_colors, 
         annot2=annot2 %>% mutate(x2 = (start - intervals2$start) / (intervals2$end - intervals2$start)), 
         lifted_annot2=lifted_annot2 %>% mutate(x1 = (start - intervals$start) / (intervals$end - intervals$start)), 
         annot2_colors=annot2_colors)
}


compute_intervals_comparison <- function(intervals, intervals2, chain, chainscore = NULL, grid_resolution=100) {
    grid_resolution <- round((intervals$end - intervals$start) / grid_resolution)

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

    return(i12)
}

is_comparison_flipped <- function(intervals_comparison) {
    i12 <- intervals_comparison

    cr <- i12 %>%
        mutate(
            mid1 = start + (end - start) / 2,
            mid2 = start1 + (end1 - start1) / 2
        ) %>%
        summarise(cr = cor(mid1, mid2, use = "pairwise.complete.obs"))

    return(cr < 0)
}

plot_intervals_comparison <- function(intervals_comparison, annotations = NULL, grid_resolution=100) {
    i12 <- intervals_comparison

    if (is_comparison_flipped(i12)) {
        i12$x2 <- 1 - i12$x2
    }
    # Set the plot parameters
    plot(1, 1, xlim = c(0, 1), ylim = c(-0.02, 1.02), type = "n", ann = FALSE, axes = FALSE)

    segments(x0 = i12$x1, y0 = 1, x1 = i12$x1, y1 = 0.99, lwd = 1)
    text(x = i12$x1, y = 0.96, labels = round(i12$start / 1e+6, 3), srt = 90, adj = c(1, 0.5), xpd = TRUE, cex = 0.5)
    segments(x0 = i12$x1, y0 = 0.75, x1 = i12$x1, y1 = 0.7, lwd = 1, col=ifelse(1:length(i12$x1)%%round(grid_resolution/10)==0,"black","grey"))
    segments(x0 = i12$x1, y0 = 0.7, x1 = i12$x2, y1 = 0.27, lwd = 1, col=ifelse(1:length(i12$x1)%%round(grid_resolution/10)==0,"black","grey"))
    segments(x0 = i12$x2, y0 = 0.27, x1 = i12$x2, y1 = 0.25, lwd = 1, col=ifelse(1:length(i12$x1)%%round(grid_resolution/10)==0,"black","grey"))
    text(x = i12$x2[i12$x2 <= 1 & i12$x2 >= 0], y = 0.2, labels = round(i12$start1[i12$x2 <= 1 & i12$x2 >= 0] / 1e+6, 3), srt = 90, adj = c(1, 0.5), xpd = TRUE, cex = 0.5)
    segments(x0 = i12$x2, y0 = 0, x1 = i12$x2, y1 = 0.01, lwd = 1)    
    if (!is.null(annotations)) {
        annot1 = annotations[["annot1"]]
        lifted_annot1 = annotations[["lifted_annot1"]]
        annot1_colors = annotations[["annot1_colors"]]
        annot2 = annotations[["annot2"]]
        lifted_annot2 = annotations[["lifted_annot2"]]
        annot2_colors = annotations[["annot2_colors"]]
        if (is_comparison_flipped(i12)) {
            annot2$x2 <- 1 - annot2$x2
            lifted_annot1$x2 <- 1 - lifted_annot1$x2
        }
        points(
            x = c(annot1$x1, lifted_annot2$x1),
            y = c(rnorm(nrow(annot1), sd = 0.001) + 1.019, rnorm(nrow(lifted_annot2), sd = 0.001) + 0.981),
            col = c(annot1_colors, annot2_colors), 
            pch = 3, cex = 0.5
        )
        points(
            x = c(annot2$x2, lifted_annot1$x2),
            y = c(rnorm(nrow(annot2), sd = 0.001) - 0.019, rnorm(nrow(lifted_annot1), sd = 0.001) + 0.019),
            col = c(annot2_colors, annot1_colors),
            pch = 3, cex = 0.5
        )

    }
}
