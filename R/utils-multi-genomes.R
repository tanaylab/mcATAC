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
            by = c('seqnames', 'strand', 'chainscore')
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
    gset_genome(chain@metadata[["genome2"]])
    ff_lifted_list <- gintervals.force_range(
        lifted_list$cont %>%
            filter(chainscore == !!chainscore) %>%
            select(chrom = seqnames, start, end, row_ID) %>%
            as.data.frame() %>%
            select(chrom, start, end, everything())
    )
    gset_genome(chain@metadata[["genome1"]])
    ff_lifted_list
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
    gset_genome(chain@metadata[["genome2"]])
    intervs2 <- gintervals.force_range(intervs2)
    gset_genome(chain@metadata[["genome1"]])
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

load_chain <- function(chain, genome1 = NULL, genome2 = NULL) {
    if (is.character(chain)) {
        cli::cli_alert_info("Loading {genome1}->{genome2} chain {.file {chain}}")
        chain <- rtracklayer::import.chain(chain)
        cli::cli_alert_success("Loaded chain {.file {chain}} successfully")
    }

    chain@metadata <- list(genome1 = genome1, genome2 = genome2)

    return(chain)
}

compute_annotation_comparison <- function(annot1, annot2, intervals, intervals2, chain_1_to_2, chain_2_to_1, selected_chain_1_to_2_chainscore = NULL, selected_chain_2_to_1_chainscore = NULL) {
    annot1 <- annot1 %>%
        gintervals.neighbors1(intervals) %>%
        filter(dist == 0)
    if(nrow(annot1)>0){
        annot1_colors <- chameleon::distinct_colors(nrow(annot1))$name
    }else{
        annot1_colors <- NULL
    }

    lifted_annot1 <- translate_intervals(as.data.frame(annot1), chain = chain_1_to_2, chainscore = selected_chain_1_to_2_chainscore) %>%
        arrange(row_ID)

    annot2 <- annot2 %>%
        filter(
            chrom == intervals2$chrom,
            start >= intervals2$start,
            end <= intervals2$end
        )
    annot2_colors <- chameleon::distinct_colors(nrow(annot1) + nrow(annot2))$name[(nrow(annot1) + 1):(nrow(annot1) + nrow(annot2))]
    lifted_annot2 <- translate_intervals(as.data.frame(annot2), chain = chain_2_to_1, chainscore = selected_chain_2_to_1_chainscore) %>%
        arrange(row_ID)
    list(
        annot1 = annot1 %>% mutate(x1 = (start - intervals$start) / (intervals$end - intervals$start)),
        lifted_annot1 = lifted_annot1 %>% mutate(x2 = (start - intervals2$start) / (intervals2$end - intervals2$start)),
        annot1_colors = annot1_colors,
        annot2 = annot2 %>% mutate(x2 = (start - intervals2$start) / (intervals2$end - intervals2$start)),
        lifted_annot2 = lifted_annot2 %>% mutate(x1 = (start - intervals$start) / (intervals$end - intervals$start)),
        annot2_colors = annot2_colors
    )
}


get_dot_mat = function(mct, mct_oc, mm_scope_peaks, oc_scope_peaks, mm2oc_chainscore=NULL, oc2mm_chainscore=NULL){
    lifted_mm <- translate_intervals(mm_scope_peaks %>% arrange(chrom, start), chain = mm2oc, chainscore = mm2oc_chainscore) %>% arrange(row_ID)
    lifted_oc <- translate_intervals(oc_scope_peaks %>% arrange(chrom, start), chain = oc2mm, chainscore = oc2mm_chainscore) %>% arrange(row_ID)
    dot_mat = matrix(nrow = dim(mm_scope_peaks)[1], ncol = dim(oc_scope_peaks)[1],0)
    gset_genome(mct_oc@genome)
    mm2oc_dot = gintervals.neighbors(oc_scope_peaks %>% arrange(start), lifted_mm, mindist = 0, maxdist = 0, na.if.notfound = T)$row_ID
    names(mm2oc_dot) = 1:dim(oc_scope_peaks)[1]
    for(i in 1:dim(oc_scope_peaks)[1]){
        if(!is.na(mm2oc_dot[i])){
            dot_mat[mm2oc_dot[i], i] = dot_mat[mm2oc_dot[as.character(i)], i] + 0.7
        }
    }
    gset_genome(mct@genome)
    oc2mm_dot = gintervals.neighbors(mm_scope_peaks %>% arrange(start), lifted_oc, mindist = 0, maxdist = 0, na.if.notfound = T)$row_ID
    names(oc2mm_dot) = 1:dim(mm_scope_peaks)[1]
    for(i in 1:dim(mm_scope_peaks)[1]){
        if(!is.na(oc2mm_dot[i])){
            dot_mat[i, oc2mm_dot[i]] = dot_mat[i, oc2mm_dot[as.character(i)]] + 0.3
        }
    }
    -1.01-dot_mat
}


get_comb_mat = function(gene_name, mct, mct_oc, mm_mat_peaks, oc_mat_peaks, mm_scope_peaks, oc_scope_peaks, flipped=FALSE, mm2oc_chainscore=NULL, oc2mm_chainscore=NULL){
    mm_mat_peaks_n = t(t(mm_mat_peaks)/mct@total_cov)
    mm_mat_peaks_nl = log2(1e-5+mm_mat_peaks_n)
    oc_mat_peaks_n = t(t(oc_mat_peaks)/mct_oc@total_cov)
    oc_mat_peaks_nl = log2(1e-5+oc_mat_peaks_n)

    gset_genome(mct@genome)
    mm_gene = t(log2(1e-5 + mct@rna_egc)[gene_name,,drop = FALSE])
    mm_atac = t(mm_mat_peaks_nl[rownames(mm_scope_peaks %>% misha.ext::intervs_to_mat(remove_intervalID = TRUE)),])
    gset_genome(mct_oc@genome)
    oc_gene = t(log2(1e-5 + mct_oc@rna_egc)[toupper(gene_name),,drop = FALSE])
    oc_atac = t(oc_mat_peaks_nl[rownames(oc_scope_peaks %>% misha.ext::intervs_to_mat(remove_intervalID = TRUE)),])
    common_mcs = intersect(rownames(mm_atac), rownames(oc_atac))

    mm_gene__mm_atac = tgstat::tgs_cor(
            mm_gene,
            mm_atac
        )
    mm_atac__mm_atac = tgstat::tgs_cor(
        mm_atac
    )
    diag(mm_atac__mm_atac) = 0
    oc_gene__oc_atac = tgstat::tgs_cor(
            oc_gene,
            oc_atac
        )
    oc_atac__oc_atac = tgstat::tgs_cor(
        oc_atac
    )
    diag(oc_atac__oc_atac) = 0

    oc_atac__mm_atac = tgstat::tgs_cor(
        oc_atac[common_mcs,],
        mm_atac[common_mcs,]
    )

    gset_genome(mct@genome)

    dot_mat = get_dot_mat(mct, mct_oc, mm_scope_peaks, oc_scope_peaks, mm2oc_chainscore=mm2oc_chainscore, oc2mm_chainscore=oc2mm_chainscore)
    
    colnames(oc_atac__mm_atac) = mm_scope_peaks$peak_name
    rownames(mm_atac__mm_atac) = mm_scope_peaks$peak_name
    colnames(mm_atac__mm_atac) = mm_scope_peaks$peak_name
    dot_mat = get_dot_mat(mct, mct_oc, mm_scope_peaks, oc_scope_peaks, mm2oc_chainscore=mm2oc_chainscore, oc2mm_chainscore=oc2mm_chainscore)

    colnames(mm_gene__mm_atac) = mm_scope_peaks$peak_name
    colnames(oc_atac__mm_atac) = mm_scope_peaks$peak_name
    rownames(mm_atac__mm_atac) = mm_scope_peaks$peak_name
    colnames(mm_atac__mm_atac) = mm_scope_peaks$peak_name
    rownames(dot_mat) = mm_scope_peaks$peak_name
    
    if(flipped){
        oc_atac__oc_atac = oc_atac__oc_atac[dim(oc_atac__oc_atac)[1]:1, dim(oc_atac__oc_atac)[2]:1]
        oc_gene__oc_atac = oc_gene__oc_atac[, dim(oc_gene__oc_atac)[2]:1, drop = FALSE]
        oc_atac__mm_atac = oc_atac__mm_atac[dim(oc_atac__mm_atac)[1]:1, ]
        dot_mat = dot_mat[dim(dot_mat)[1]:1,]
    }
    
    oc_gene_oc_atac_padded = t(cbind(0, oc_gene__oc_atac))
    rownames(oc_gene_oc_atac_padded)[1] = gene_name

    comb_mat = cbind(
            rbind(
                cbind(oc_gene_oc_atac_padded, rbind(mm_gene__mm_atac, oc_atac__mm_atac)), 
                cbind(0,mm_atac__mm_atac)
            ),
            rbind(0,
                oc_atac__oc_atac, 
                matrix(nrow = dim(mm_scope_peaks)[1], ncol = dim(oc_scope_peaks)[1],0))
        )
    comb_mat
}


compute_intervals_comparison <- function(intervals, intervals2, chain, chainscore = NULL, grid_resolution = 100) {
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
        left_join(i12 %>% rename(chrom1 = chrom, start1 = start, end1 = end), by = "row_ID") %>%
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

plot_intervals_comparison <- function(intervals_comparison, annotations = NULL, correlations = NULL, grid_resolution = 100) {
    i12 <- intervals_comparison

    if (is_comparison_flipped(i12)) {
        i12$x2 <- 1 - i12$x2
    }
    # Set the plot parameters
    plot(1, 1, xlim = c(0, 1), ylim = c(-0.02, 1.02), type = "n", ann = FALSE, axes = FALSE)

    segments(x0 = i12$x1, y0 = 1, x1 = i12$x1, y1 = 0.99, lwd = 1)
    #text(x = i12$x1, y = 0.96, labels = round(i12$start / 1e+6, 3), srt = 90, adj = c(1, 0.5), xpd = TRUE, cex = 0.5)
    #segments(x0 = i12$x1, y0 = 0.75, x1 = i12$x1, y1 = 0.7, lwd = 1, col = ifelse(1:length(i12$x1) %% round(grid_resolution / 10) == 0, "black", "grey"))
    segments(x0 = i12$x1, y0 = 0.99, x1 = i12$x2, y1 = 0.01, lwd = 1, col = ifelse(1:length(i12$x1) %% round(grid_resolution / 10) == 0, "black", "grey"))
    #segments(x0 = i12$x2, y0 = 0.27, x1 = i12$x2, y1 = 0.25, lwd = 1, col = ifelse(1:length(i12$x1) %% round(grid_resolution / 10) == 0, "black", "grey"))
    #text(x = i12$x2[i12$x2 <= 1 & i12$x2 >= 0], y = 0.2, labels = round(i12$start1[i12$x2 <= 1 & i12$x2 >= 0] / 1e+6, 3), srt = 90, adj = c(1, 0.5), xpd = TRUE, cex = 0.5)
    segments(x0 = i12$x2, y0 = 0, x1 = i12$x2, y1 = 0.01, lwd = 1)
    if (!is.null(correlations)) {
    }
    if (!is.null(annotations)) {
        annot1 <- annotations[["annot1"]]
        lifted_annot1 <- annotations[["lifted_annot1"]]
        annot1_colors <- annotations[["annot1_colors"]]
        annot2 <- annotations[["annot2"]]
        lifted_annot2 <- annotations[["lifted_annot2"]]
        annot2_colors <- annotations[["annot2_colors"]]
        if (is_comparison_flipped(i12)) {
            annot2$x2 <- 1 - annot2$x2
            lifted_annot1$x2 <- 1 - lifted_annot1$x2
        }
        points(
            x = c(annot1$x1, lifted_annot2$x1),
            y = c(rnorm(nrow(annot1), sd = 0.001) + 1.019, rnorm(nrow(lifted_annot2), sd = 0.001) + 0.981),
            col = c(annot1_colors, annot2_colors),
            pch = 3, cex = 1
        )
        points(
            x = c(annot2$x2, lifted_annot1$x2),
            y = c(rnorm(nrow(annot2), sd = 0.001) - 0.019, rnorm(nrow(lifted_annot1), sd = 0.001) + 0.019),
            col = c(annot2_colors, annot1_colors),
            pch = 3, cex = 1
        )
    }
}

#' Generate an hclust object from a data.frame with cell_type, group and order columns
#'
#' @param df data.frame with cell_type, group and order
#'
#' @return hclust object with cell_type as labels and order as order of the leaves
#'
#' @examples
#' df <- data.frame(
#'     cell_type = c(
#'         "PGC", "Surface ectoderm", "Neural crest",
#'         "Neural tube/Floor plate", "Forebrain/Midbrain/Hindbrain", "Caudal neural plate",
#'         "Rostral neural plate", "Neural plate boundary", "Definitive ectoderm",
#'         "Epiblast", "Primitive streak", "Caudal epiblast", "Tail bud - neural",
#'         "Tail bud - mesoderm", "Early nascent mesoderm"
#'     ), order = 1:15,
#'     group = c(
#'         "ecto", "ecto", "ecto", "ecto", "ecto", "ecto",
#'         "ecto", "ecto", "ecto", "early", "early", "early", "early",
#'         "Meso", "Meso"
#'     )
#' )
#'
#' group_to_hclust(df)
#'
#' @export
group_to_hclust <- function(df) {
    # Create an empty matrix of 1s
    distance_matrix <- matrix(1, nrow = nrow(df), ncol = nrow(df))
    rownames(distance_matrix) <- df$cell_type
    colnames(distance_matrix) <- df$cell_type

    # set the distance between elements the same group to 0 (using df$group)
    for (grp in unique(df$group)) {
        distance_matrix[df$group == grp, df$group == grp] <- 0
    }

    # Hierarchical clustering
    hc <- hclust(as.dist(distance_matrix))

    # Label leaves with cell_type
    hc$labels <- df$cell_type

    hc <- as.hclust(reorder(as.dendrogram(hc), df$order, agglo.FUN = mean))

    return(hc)
}

get_heatmap_idx <- function(x, y, m) {
    n_x <- nrow(m)
    n_y <- ncol(m)
    x_idx <- floor(x * n_x) + 1
    y_idx <- n_y - max(floor(y * n_y), 0)
    return(c(x_idx, y_idx))
}
