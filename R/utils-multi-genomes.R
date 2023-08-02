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
        summarise(start = min(start), end = max(end), widthsum = sum(width))
    test3$alt_start <- test3$start - 1
    test3$match <- test3$widthsum / (test3$end - test3$alt_start)
    atest2 <- as.data.frame(test2) %>% left_join(ungroup(test3) %>% select(seqnames, strand, chainscore, lstart = alt_start, lend = end, widthsum, match) %>% mutate(lrow_id = row_number()), by = join_by(seqnames, strand, chainscore))
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
        )
    btest2$alt_start <- btest2$start - 1
    btest2$match <- btest2$widthsum / (btest2$end - btest2$alt_start)
    return(list(full = test2, gapped = btest2, cont = test3))
}

translate_intervals <- function(intervals, chain) {
    # regioneR::toGRanges(mm_tad)
    # lifted_list_mm_tad_to_oc = custom_lift2(gr_mm10_tad, mm2oc)
    r <- regioneR::toGRanges(intervals)
    lifted_list <- custom_lift2(r, chain)
    lifted_list$cont %>%
        arrange(desc(chainscore)) %>%
        ungroup() %>%
        slice(1) %>%
        select(chrom = seqnames, start, end) %>%
        as.data.frame()
}
