#' Plot a scatter of gene expression vs an atac profile
#'
#' @param mc_atac a McATAC object
#' @param gene name of the gene to plot
#' @param atac_promoter name of the promoter to plot. By default this would be the rna gene \code{gene}.
#' @param mc_rna a \code{metacell1} 'mc' object or a \code{metacells} metacell UMI matrix (a matrix where each row is a gene and each column is a metacell). Can be NULL if \code{mc_atac} already contains the gene expression data (added by \code{add_mc_rna}).
#' @param peak name of the peak to plot. If NULL - the promoter of \code{atac_promoter} would be shown. You can get the peak names from \code{atac@peaks$peak_name}
#' @param normalize_atac whether to use normalized atac profiles (default: TRUE)
#' @param eps_rna added regularization when calculating log expression (Default: 1e-5). Promoter ATAC signal has been observed empirically to often be linearly correlated with log of gene expression.
#' @param plot_object_id plot the object id of the \code{mc_atac} object in the bottom left corner of the plot (default: TRUE)
#'
#' @return a ggplot object with the scatter plot
#'
#' @examples
#' \dontrun{
#' p1 <- plot_atac_rna(mc_atac = pbmc_atac_mc, mc_rna = pbmc_rna_mc, gene = "CD4")
#'
#' p2 <- plot_atac_rna(mc_atac = pbmc_atac_mc, mc_rna = pbmc_rna_mc, gene = "CD4", atac_promoter = "CD8")
#'
#' # Plot gene vs some peak of interest
#' peak <- mc_atac@peaks[1]$peak_name
#' p3 <- plot_atac_rna(mc_atac = pbmc_atac_mc, mc_rna = pbmc_rna_mc, gene = "CD4", peak = peak)
#' }
#'
#' @inheritParams get_promoter_peaks
#' @export
plot_atac_rna <- function(mc_atac, gene, atac_promoter = gene, mc_rna = NULL, peak = NULL, max_dist_to_promoter_peak = 5e+2, normalize_atac = TRUE, eps_rna = 1e-5, tss_intervals = "intervs.global.tss", plot_object_id = TRUE) {
    assert_atac_object(mc_atac, class = "McATAC")

    if (!has_rna(mc_atac) && is.null(mc_rna)) {
        cli_abort("No gene expression data in the McATAC object. Use {.code add_mc_rna()} to add it or provide the {.field mc_rna} parameter.")
    }

    if (!is.null(mc_rna)) {
        mc_atac <- add_mc_rna(mc_atac, mc_rna)
    }

    rv <- log2(mc_atac@rna_egc[gene, ] + eps_rna)

    if (is.null(peak)) {
        peaks <- get_promoter_peaks(mc_atac@peaks, atac_promoter, max_dist_to_promoter_peak = max_dist_to_promoter_peak, tss_intervals = tss_intervals)

        if (is.null(peak) && length(peaks) == 0) {
            cli_abort("{.field peak} was not supplied and no relevant promoter peak was found for gene {.val {atac_promoter}}. Check if your ATAC matrix should have a peak for this gene's promoter.")
        } else if (length(peaks) > 1) {
            cli_alert("The gene {.val {atac_promoter}} has multiple ({.val {length(peaks)}}) peaks within {.val {max_dist_to_promoter_peak}} bp of its TSS. Summing the ATAC signal from all of them.")
        }

        if (length(peaks) > 1) {
            peak_range <- misha.ext::convert_10x_peak_names_to_misha_intervals(peaks) %>%
                summarise(chrom = chrom[1], start = min(start), end = max(end)) %>%
                peak_names()
            peak_str <- glue("{peak_range} ({length(peaks)} peaks)")
        } else {
            peak_str <- peaks
        }
    } else {
        if (peak %!in% mc_atac@peaks$peak_name) {
            cli_abort("{.val {peak}} is not a peak in the McATAC object.")
        }
        peaks <- peak
        peak_str <- peak
    }

    if (normalize_atac) {
        av <- colSums(mc_atac@egc[peaks, , drop = FALSE])
        ylab <- "ATAC (normalized)"
    } else {
        av <- colSums(mc_atac@mat[peaks, , drop = FALSE])
        ylab <- "ATAC (not normalized)"
    }

    # take only metacells which exist in both atac and rna
    both_mcs <- intersect(names(av), names(rv))
    av <- av[both_mcs]
    rv <- rv[both_mcs]

    cor_r2 <- cor.test(rv, av)$estimate

    df <- tibble(
        metacell = names(av),
        atac = av,
        rna = rv
    )

    if (all(has_name(mc_atac@metadata, c("cell_type", "color")))) {
        df <- df %>%
            left_join(mc_atac@metadata %>%
                mutate(metacell = as.character(metacell)) %>%
                select(metacell, cell_type, color), by = "metacell")
        gg <- ggplot(df, aes(x = rna, y = atac, fill = cell_type)) +
            geom_point(shape = 21) +
            scale_fill_manual(name = "Cell type", values = get_cell_type_colors(mc_atac@metadata))
    } else {
        gg <- ggplot(df, aes(x = rna, y = atac)) +
            geom_point()
    }

    if (!is.null(peak)) {
        title <- gene
        subtitle <- glue("ATAC of {peak_str} vs. RNA, R^2 = {round(cor_r2, digits=2)}")
        caption <- ggplot2::waiver()
    } else if (atac_promoter == gene) {
        title <- gene
        subtitle <- glue("ATAC of promoter vs. RNA, R^2 = {round(cor_r2, digits=2)}")
        caption <- glue("Promoter: {peak_str}")
    } else { # atac_promoter != gene and no peak was supplied
        title <- glue("ATAC of {atac_promoter} promoter vs. {gene} RNA")
        subtitle <- glue("R^2 = {round(cor_r2, digits=2)}")
        caption <- glue("Promoter: {peak_str}")
    }

    if (plot_object_id) {
        caption <- paste(caption, glue("object id: {mc_atac@id}"), sep = "\n")
    }

    gg <- gg +
        labs(
            title = title,
            subtitle = subtitle,
            x = "log2(gene expression)",
            y = ylab,
            caption = caption
        ) +
        theme(
            plot.subtitle = ggtext::element_markdown(),
            plot.caption = element_text(hjust = 0),
            aspect.ratio = 1
        )

    return(gg)
}

#' Plot a correlation matrix of ATAC metacells
#'
#' @param mc_atac McATAC object
#' @param sp_f whether to use Spearman correlation (default) or Pearson
#'
#' @return p a pheatmap of ATAC metacell correlations
#'
#' @examples
#' \dontrun{
#' p1 <- plot_atac_atac_cor(my_atac_mc)
#' p2 <- plot_atac_atac_cor(my_atac_mc, sp_f = F)
#' }
#' @export
plot_atac_atac_cor <- function(mc_atac, sp_f = TRUE) {
    if (all(!grepl("cell_type", colnames(mc_atac@metadata)))) {
        if (all(!grepl("cluster_k_", colnames(mc_atac@peaks)))) {
            k <- round(ncol(mc_atac@mat) / 10)
            message(glue::glue('There is no "cell_type" or "cluster" field in metadata, clustering with k == {k}'))
            mc_atac <- gen_atac_mc_clust(mc_atac, k = k, use_prior_annot = F)
            clust_vec <- unlist(mc_atac@metadata[, paste0("cluster_k_", k)])
        } else {
            clust_vec <- unlist(mc_atac@metadata[, grep("cluster_k_", colnames(mc_atac@metadata))[[1]]])
        }
    } else {
        clust_vec <- unlist(mc_atac@metadata[, "cell_type"])
    }
    col_annot <- as.data.frame(list("cluster" = clust_vec))
    rownames(col_annot) <- 1:nrow(mc_atac@metadata)
    if (any(grepl("color", colnames(mc_atac@metadata)))) {
        col_key <- unique(mc_atac@metadata[, c("cell_type", "color")])
        ann_colors <- list("cluster" = setNames(col_key$color, col_key$cell_type))
    } else {
        ann_colors <- list("cluster" = setNames(sample(
            grep("white|gray|grey", colors(), v = T, inv = T),
            length(unique(clust_vec))
        ), unique(clust_vec)))
    }
    cor_mat <- tgs_cor(mc_atac@mat, spearman = sp_f, pairwise.complete.obs = TRUE)
    p <- pheatmap::pheatmap(cor_mat[order(clust_vec), order(clust_vec)],
        cluster_rows = FALSE, cluster_cols = FALSE,
        show_rownames = FALSE, show_colnames = FALSE,
        annotation_col = col_annot, annotation_colors = ann_colors
    )
    return(p)
}

#' Plot a cross-correlation matrix between RNA metacells and ATAC metacell scores
#'
#' @description use peak gene annotations to match between RNA metacells and ATAC metacells and plot
#' a cross-correlation matrix.
#'
#' @param mc_atac McATAC object
#' @param rna_mat an RNA metacell count matrix, where metacells are in columns and genes are in rows
#' @param gene_field (optional) either \code{closest_tss} or \code{closest_exon_gene} -- field names in \code{mc_atac@peaks} which contain gene names. If NULL - the peaks would be
#' transformed to promoter peaks and the gene names would be taken from the promoter gene names.
#' @param tss_dist (optional) maximal absolute distance to a TSS to be considered a promoter peak
#'
#' @return a pheatmap object
#'
#' @examples
#' \dontrun{
#' # Plot correlation of ATAC promoter peaks vs. log2 gene expression fraction (regularized)
#' arc_prom <- plot_atac_rna_cor(mc_atac = mc_atac, rna_mat = log2(atac_mc@rna_egc + 1e-05))
#'
#' # Plot correlation of all available ATAC peaks (whose nearest TSS is of an expressed gene) vs. log2 gene expression fraction (regularized) of that gene
#' arc_tss <- plot_atac_rna_cor(mc_atac = mc_atac, rna_mat = log2(atac_mc@rna_egc + 1e-05), gene_field = "closest_tss")
#'
#' # Plot correlation of all available ATAC peaks (whose nearest exon is of an expressed gene) vs. log2 gene expression fraction (regularized) of that gene
#' arc_exon <- plot_atac_rna_cor(mc_atac = mc_atac, rna_mat = log2(atac_mc@rna_egc + 1e-05), gene_field = "closest_exon_gene")
#' }
#' @export
plot_atac_rna_cor <- function(mc_atac, rna_mat, gene_field = NULL, tss_dist = 5e+2) {
    if (!is.null(gene_field) && has_name(mc_atac@peaks, gene_field)) {
        gb <- intersect(unique(unlist(mc_atac@peaks[, gene_field])), rownames(rna_mat))
        ulgf <- unlist(mc_atac@peaks[, gene_field])
        non_na_inds <- which(!is.na(ulgf) & ulgf %in% gb)
        genes_of_peaks <- unlist(mc_atac@peaks[non_na_inds, gene_field])
        atac_mat <- mc_atac@egc[non_na_inds, ]
    } else {
        tss <- gintervals.load("intervs.global.tss")
        nei_peak_prom <- gintervals.neighbors(as.data.frame(mc_atac@peaks), tss, mindist = -tss_dist, maxdist = tss_dist)
        prom_peaks_genes <- nei_peak_prom[, c("peak_name", "geneSymbol", "dist")]
        min_inds <- sapply(unique(prom_peaks_genes$geneSymbol), function(u) {
            inds <- which(prom_peaks_genes$geneSymbol == u)
            res <- inds[which.min(prom_peaks_genes$dist[inds])]
            return(res)
        })
        atac_mat <- mc_atac@egc[prom_peaks_genes$peak_name[min_inds], ]
        genes_of_peaks <- prom_peaks_genes$geneSymbol[min_inds]
        gb <- intersect(genes_of_peaks, rownames(rna_mat))
        gp_in <- genes_of_peaks %in% gb
        atac_mat <- atac_mat[gp_in, ]
        genes_of_peaks <- genes_of_peaks[gp_in]
    }
    atac_mat_ord <- atac_mat

    rownames(atac_mat_ord) <- genes_of_peaks
    rna_match <- match(genes_of_peaks, rownames(rna_mat))
    rna_mat_ord <- rna_mat[rna_match, ]

    atac_rna_cor <- tgs_cor(atac_mat_ord, rna_mat_ord, spearman = T)
    if (all(has_name(mc_atac@metadata, c("metacell", "cell_type")))) {
        col_annot <- tibble::column_to_rownames(mc_atac@metadata[, c("metacell", "cell_type")], "metacell")
        ann_colors <- list("cell_type" = setNames(unlist(mc_atac@metadata[, "color"]), unlist(mc_atac@metadata[, "cell_type"])))
    } else {
        cli_alert_info("No metacell annotation detected. Clustering metacells.")
        mc_annot <- generate_mc_annotation(mc_atac = mc_atac)
        col_annot <- mc_annot[[1]]
        ann_colors <- mc_annot[[2]]
    }
    p <- pheatmap::pheatmap(atac_rna_cor,
        cluster_rows = TRUE, cluster_cols = TRUE,
        show_rownames = TRUE, show_colnames = TRUE,
        annotation_col = col_annot,
        annotation_row = col_annot,
        annotation_colors = ann_colors, silent = FALSE
    )
    return(p)
}

#' Plot normalized accessibility of peaks over metacells, ordered by clustering
#'
#' @param mc_atac McATAC object
#' @param mc_atac_clust output of \code{gen_atac_mc_clust} (for meaningful visuals, make sure this is ordered by cluster)
#' @param peak_clust output of \code{gen_atac_peak_clust} or other clustering of peaks (for meaningful visuals, make sure this is ordered by cluster)
#' @param peak_annotation (optional) a list of a named vector and a dataframe conforming to the pheatmap \code{annotation_colors} and \code{annotation_row} conventions
#' @param filename (optional) path and filename of where to save the figure; if unspecified, figure isn't saved
#' @param dev (optional; default - png) graphical device to save figure with
#' @param colors (optional) colorRampPalette vector of colors for scaling colors in heatmap
#'
#' @inheritDotParams save_pheatmap
#'
#' @return a pheatmap figure.
#' @examples
#' \dontrun{
#' peak_clust <- gen_atac_peak_clust(my_mcatac, 16)
#' plot_atac_peak_map(my_mcatac, order(my_mcatac@metadata$cell_type), order(peak_clust), "./my_figures/my_mcatac_heatmap.png")
#'
#' ## Peak annotation example
#' is_dyn <- setNames(as.numeric(rownames(my_mcatac@mat) %in% my_mcatac_dynamic_peaks_only@peaks$peak_name), rownames(my_mcatac@mat))
#' pa1 <- list("is_dyn" = setNames(c("black", "red"), c(0, 1)))
#' pa2 <- tibble::column_to_rownames(enframe(is_dyn, name = "peak_name", value = "is_dyn"), "peak_name")
#' pa <- list(pa1, pa2)
#' plot_atac_peak_map(my_mcatac, mc_atac_clust = order(my_mcatac@metadata$cell_type), peak_annotation = pa)
#' }
#' @export
plot_atac_peak_map <- function(mc_atac, mc_atac_clust = NULL, peak_clust = NULL,
                               peak_annotation = NULL, filename = NULL,
                               dev = png, main = mc_atac@id,
                               colors = colorRampPalette(c("blue4", "white", "red4"))(100),
                               ...) {
    if (is.null(mc_atac_clust)) {
        if (all(has_name(mc_atac@metadata, c("metacell", "cell_type")))) {
            mc_atac_clust <- deframe(mc_atac@metadata[, c("metacell", "cell_type")])
        } else {
            hc <- hclust(tgs_dist(t(mc_atac@fp)))
            mc_atac_clust <- match(as.numeric(hc$labels), hc$order)
        }
    }
    lmcoefs <- setNames(c(15.9252182670872, 0.00089588822113778), c("(Intercept)", "x"))
    row_annot <- NULL
    if (all(has_name(mc_atac@metadata, c("metacell", "cell_type")))) {
        col_annot <- tibble::column_to_rownames(mc_atac@metadata[, c("metacell", "cell_type")], "metacell")
        ann_colors <- list("cell_type" = setNames(unlist(mc_atac@metadata[, "color"]), unlist(mc_atac@metadata[, "cell_type"])))
    } else {
        mc_annot <- generate_pheatmap_annotation(mc_atac_clust, feature_type = "metacell", feature_annotation = "cluster")
        col_annot <- mc_annot[[1]]
        ann_colors <- mc_annot[[2]]
    }
    mca_lfc <- mc_atac@fp
    brks <- c(
        seq(min(mca_lfc), 0, l = 50),
        seq(0.01 * (max(mca_lfc) - min(mca_lfc)), max(mca_lfc), l = 51)
    )
    colnames(mca_lfc) <- 1:ncol(mca_lfc)
    if (!is.null(peak_annotation)) {
        if (length(peak_annotation) != 2 || sapply(peak_annotation, class) != c("list", "data.frame")) {
            cli_abort("{.field peak_annotation} must be a list containing: 1. a list containing a named vector; 2. a dataframe of one column. See examples and compare with pheatmap docs/examples.")
        }
        if (!class(peak_annotation[[2]][, 1]) %in% c("numeric", "character")) {
            cli_abort("Peak annotation column in peak annotation dataframe must be of a numeric or character class")
        }
        row_annot <- peak_annotation[[2]]
        ann_colors[names(peak_annotation[[1]])] <- peak_annotation[[1]]
    } else {
        if (is.null(peak_clust)) {
            cli_alert_info("No peak clustering specified. Generating peak clusters.")
            peak_clust <- gen_atac_peak_clust(mc_atac, clustering_algoritm = "louvain")
        }
        peak_annot <- generate_pheatmap_annotation(peak_clust, feature_type = "peak", feature_annotation = "cluster")
        row_annot <- peak_annot[[1]]
        ann_colors[names(peak_annot[[1]])] <- peak_annot[[2]]
    }
    mc_atac_clust <- order(mc_atac_clust)
    peak_clust <- order(peak_clust)
    cli_alert_info("Expected time to plot is roughly {.val {round(lmcoefs[[1]] + length(peak_clust)*lmcoefs[[2]], 0)}}s")
    pp <- pheatmap::pheatmap(mca_lfc[peak_clust, mc_atac_clust],
        annotation_col = subset(col_annot, select = cell_type),
        annotation_legend = FALSE,
        annotation_colors = ann_colors,
        annotation_row = row_annot, main = main,
        color = colors, breaks = brks, cluster_cols = FALSE, cluster_rows = FALSE, show_colnames = FALSE, show_rownames = FALSE
    )
    if (!is.null(filename)) {
        save_pheatmap(pp, filename = filename, dev = dev, ...)
    }
    return(pp)
}


#' Plot metacell tracks around locus
#'
#' @param tracks (optional) all tracks to plot
#' @param gene (optional) which gene to plot around
#' @param intervals (optional) what genomic interval to plot
#' @param mc_rna (optional) RNA metacell object to extract gene expression data from
#' @param atac (optional) ScATAC, McATAC or PeakIntervals object from which to extract peaks in locus
#' @param track_regex (optional) regular expression for matching tracks to plot
#' @param rna_legc_eps (optional) what regularization value to add to mc_rna@e_gc when calculating log
#' @param gene_annot (optional) whether to add gene annotations
#' @param silent (optional) whether to print generated plot
#' @param annotation_row pheatmap-format annotation
#' @param annotation_col pheatmap-format annotation
#' @param annotation_colors pheatmap-format annotation
#' @param iterator (optional) misha iterator
#' @param extend (optional) how much to extend \code{intervals} by on each side
#' @param colors (optional) colorRampPalette vector of colors for scaling colors in heatmap
#'
#' @inheritParams ComplexHeatmap
#' @inheritDotParams save_pheatmap
#'
#' @return a ComplexHeatmap figure
#' @examples
#' \dontrun{
#'     library(metacell)
#'     scdb_init("scdb")
#'     mc_rna = scdb_mc('rna_w_color_key')
#'     mcmd = readr::read_csv('./data/mcmd.csv')
#'     color_key = unique(mcmd[,c('st', 'color')])
#'     colnames(color_key) = c('cell_type', 'color')
#'     col_annot = data.frame(cell_type = mcmd$st)
#'     rownames(col_annot) = mcmd$mc
#'     ann_colors = list(cell_type = setNames(color_key$color, color_key$cell_type))
#'     plot_tracks_at_locus(tracks = pbmc_tracks, extend = 5e+4,
#'                     gene = "ATF3", 
#'                     mc_rna = mc_rna, 
#'                     gene_annot = T, 
#'                     order_rows = T,
#'                     annotation_row = col_annot, 
#'                     annotation_colors = ann_colors)
#' }
#' @export
plot_tracks_at_locus <- function(tracks = NULL, 
                                gene = NULL, 
                                intervals = NULL, 
                                mc_rna = NULL,
                                gene_feature = "exon",
                                track_regex = NULL, 
                                iterator = 100,
                                extend = 0, 
                                order_rows = FALSE,
                                row_order = NULL,
                                gene_annot = TRUE,
                                annotation_row = NULL,
                                annotation_col = NULL,
                                annotation_colors = NULL,
                                silent = FALSE,
                                rna_legc_eps = 1e-5,
                                colors = colorRampPalette(c("white", "darkblue","red"))(100),
                               ...) {
    if (is.null(tracks)) {
        if (is.null(track_regex)) {
            cli_abort("Must specify either {.var tracks} or {.var track_regex")
        }
        else {
            tracks <- gtrack.ls(track_regex)
            if (length(tracks) == 0) {
                cli_abort("No tracks matching {.var track_regex} were found.")
            }
        }
    }
    if (!is.null(gene)) {
        if (gene_feature %!in% c('tss', 'exon')) {
            cli_abort("{.var gene_feature} should be either 'tss' or 'exon'")
        }
        feature_df <- gintervals.load(glue::glue("intervs.global.{gene_feature}"))
        feature_df <- dplyr::filter(feature_df, !grepl("_", feature_df$chrom))
        gene_features <- dplyr::filter(feature_df, geneSymbol == gene)
        if (nrow(gene_features) == 0) {
            cli_abort("No {.val gene_feature}s matching {.val gene} were found. Maybe check \\
                    that this feature exists in the field {.field geneSymbol} of gintervals.load('intervs.global.{gene_feature}')")
        }
        intervals <- gintervals(unique(gene_features$chrom), min(gene_features$start) - extend, max(gene_features$end) + extend)
    } else {
        intervals <- dplyr::mutate(intervals, start = start - extend, end = end + extend)
    }
    if (gene_annot) {
        gene_annots <- make_gene_annot(intervals, iterator)
    } else {
        gene_annots <- NULL
    }
    if (order_rows) {
        if (!is.null(annotation_row) && has_name(annotation_row, "cell_type")) {
            if (!is.null(row_order)) {
                cli_alert_info("Both {.var order_rows} == TRUE and {.var row_order} specified. Ordering rows and ignoring {.var row_order}.")    
            }
            row_order <- order(annotation_row[,"cell_type"])
        }
        else {
            cli_alert_info("No appropriate metacell annotation provided for ordering tracks. Tracks will not be ordered.")
        }
    } else {
        row_order <- 1:length(tracks)
    }
    if (!is.null(mc_rna)) {
        if (length(tracks) != ncol(mc_rna@e_gc)) {
            cli_abort("Number of tracks and number of RNA metacells do not match.")
        }
        rna_ha <- make_rna_annot(mc_rna, gene, rna_legc_eps, row_order)
    }
    mc_gene_vals <- gextract(tracks, intervals = intervals, iterator = iterator)
    mat <- t(subset(mc_gene_vals, select = -c(chrom, start, end, intervalID)))
    rownames(mat) <- 1:nrow(mat)
    mat_n <- apply(mat, 2, function(x) {x[is.na(x)] <- 0; return(x)})
    if (!is.null(annotation_row)) {
        if (is.null(annotation_colors)) {
            cli_abort("Must specify {.var annotation_colors} if {.var annotation_row} is specified.")
        }
        ct_ha <- HeatmapAnnotation(cell_type = unlist(annotation_row[row_order,]), which = 'row', col = annotation_colors)
    } else {
        ct_ha <- HeatmapAnnotation(cell_type = row_order, which = 'row')
    }
    ch <- Heatmap(mat_n[row_order,], 
                    use_raster = F,
                    top_annotation = gene_annots, 
                    left_annotation = ct_ha, 
                    right_annotation = rna_ha, 
                    col = colors,
                    cluster_rows = F, 
                    cluster_columns = F, 
                    show_row_names = F, 
                    show_column_names = F)
    return(ch)
}


#' @param intervals what genomic interval to generate gene annotations for
#' @param iterator misha iterator of the associated matrix
#' @return gene annotations: a binary vector for exon locations and associated text labels for starts and ends of transcripts
#' @noRd
make_gene_annot <- function(intervals, iterator) {
    file_path <- file.path(dirname(GROOT), "annots", "refGene.txt")
    refgene <- tgutil::fread(file_path, col.names = c("bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames"))
    rg_here <- dplyr::filter(refgene, chrom == as.character(intervals$chrom), txStart >= intervals$start, txEnd <= intervals$end)
    make_binary_exon_vec <- function(exon_intervals, genomic_intervals, iterator) {
        vec = matrix(0, nrow = 1, ncol = round((genomic_intervals[['end']] - genomic_intervals[['start']])/iterator + 1))
        exon_bins = apply(exon_intervals[,2:3], 2, function(x) {
            round((x - genomic_intervals[['start']])/iterator)
        })
        if (nrow(exon_intervals) > 1) {
            for (i in 1:nrow(exon_intervals)) {
                vec[exon_bins[i,1]:exon_bins[i,2]] = 1
            }
        } else {
            vec[exon_bins[[1]]:exon_bins[[2]]] = 1
        }
        return(as.numeric(vec))
    }
    if (nrow(rg_here) > 0) {
        vecs <- t(apply(rg_here, 1, function(x) {
                        starts = as.numeric(stringi::stri_remove_empty(unlist(stringr::str_split(x[['exonStarts']], ','))))
                        ends = as.numeric(stringi::stri_remove_empty(unlist(stringr::str_split(x[['exonEnds']], ','))))
                        df <- data.frame('chrom' = rep(x[['chrom']], length(starts)), 'start' = starts, 'end' = ends)
                        vec <- make_binary_exon_vec(df, intervals, iterator)
                        return(vec)
                    }))
        vecsums = as.numeric(pmin(colSums(vecs), 1))
        print(length(vecsums))
        genes_ha <- HeatmapAnnotation(gene_names = anno_mark(labels = rep(rg_here$name2, 2),
                                                        at = c(round((rg_here$txEnd - intervals$start)/iterator), 
                                                                round((rg_here$txStart - intervals$start)/iterator)),
                                                        ),
                                        genes = vecsums, col = list(genes = setNames(c('white', 'black'), c(0,1))),
                                        which = 'column')
    } else {
        genes_ha <- NULL
    }
    return(genes_ha)
}

#' @param mc_rna RNA metacell object
#' @param gene gene of interest
#' @param rna_legc_eps regularization value when taking log of mc_rna@e_gc
#' @param row_order order of rows
#' @noRd
make_rna_annot <- function(mc_rna, gene, rna_legc_eps, row_order) {
    rna_vals <- log2(mc_rna@e_gc[gene,] + rna_legc_eps)
    rna_vals <- rna_vals - median(rna_vals)
    rna_vals <- rna_vals[row_order]
    rna_ha <- HeatmapAnnotation('rna\nlegc\nminus\nmedian' = anno_barplot(rna_vals), which = 'row')
    return(rna_ha)
}