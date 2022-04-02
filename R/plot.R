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
        caption <- ""
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
#' @param gene_field either \code{closest_tss} or \code{closest_exon_gene} -- field names in \code{mc_atac@peaks} which contain gene names. If NULL - the peaks would be
#' @param tss_dist (optional) maximal absolute distance to a TSS to be considered a promoter peak
#' transformed to promoter peaks and the gene names would be taken from the promoter gene names.
#'
#' @export
plot_atac_rna_cor <- function(mc_atac, rna_mat, gene_field, tss_dist = 5e+2) {
    if (has_name(mc_atac@peaks, gene_field)) {
        gb <- intersect(unique(unlist(mc_atac@peaks[,gene_field])), rownames(rna_mat))
        non_na_inds <- which(!is.na(unlist(mc_atac@peaks[,gene_field])))
        genes_of_peaks <- unlist(mc_atac@peaks[non_na_inds,gene_field])
        atac_mat <- mc_atac@mat[non_na_inds,]
    }
    else {
        nei_peak_prom <- gintervals.neighbors(intervals, tss, mindist = -tss_dist, maxdist = tss_dist)
        prom_peaks_genes <- nei_peak_prom[,c('peak_name', 'geneSymbol', 'dist')]
        min_inds <- sapply(unique(prom_peaks_genes$geneSymbol), function(u) {
            inds <- which(prom_peaks_genes$geneSymbol == u)
            return(inds[which.min(prom_peaks_genes$dist)])
        })
        atac_mat <- atac_mat[prom_peaks_genes$peak_name[min_inds],]
        genes_of_peaks <- prom_peaks_genes$geneSymbol[min_inds]
        gb <- intersect(genes_of_peaks, rownames(rna_mat))
    }
    gb <- gb[order(match(gb, rownames(rna_mat)))]
    atac_mat_ord <- atac_mat[order(match(genes_of_peaks, gb)),]
    rna_mat_ord <- rna_mat[match(genes_of_peaks, gb),]
    atac_rna_cor <- tgs_cor(atac_mat_ord, rna_mat_ord, spearman = T)
    if (all(has_name(mc_atac@metadata, c("metacell", "cell_type")))) {
        col_annot <- tibble::column_to_rownames(mc_atac@metadata[, c("metacell", "cell_type")], "metacell")
        ann_colors <- list("cell_type" = setNames(unlist(mc_atac@metadata[, "color"]), unlist(mc_atac@metadata[, "cell_type"])))
    }
    else {
        cli_alert_info('No metacell annotation detected. Clustering metacells.')

    }
    return(atac_rna_cor)
}

#' Plot normalized accessibility of peaks over metacells, ordered by clustering
#'
#' @param mc_atac McATAC object
#' @param mc_atac_clust output of \code{gen_atac_mc_clust} (for meaningful visuals, make sure this is ordered by cluster)
#' @param peak_clust output of \code{gen_atac_peak_clust} or other clustering of peaks (for meaningful visuals, make sure this is ordered by cluster)
#' @param peak_annotation (optional) a list of a named vector and a dataframe conforming to the pheatmap \code{annotation_colors} and \code{annotation_row} conventions
#' @param filename (optional) path and filename of where to save the figure; if unspecified, figure isn't saved
#' @param dev (optional; default - png) graphical device to save figure with
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
plot_atac_peak_map <- function(mc_atac, mc_atac_clust, peak_clust,
                               peak_annotation = NULL,
                               filename = NULL,
                               dev = png,
                               main = mc_atac@id,
                               clrs = colorRampPalette(c("blue4", "white", "red4"))(100),
                               ...) {
    if (is.null(mc_atac_clust)) {
        cli_abort("Must specify clustering of metacells (e.g. using {.code gen_atac_mc_clust})")
    }
    if (is.null(peak_clust)) {
        cli_abort("Must specify clustering of peaks (e.g. using {.code gen_atac_peak_clust})")
    }
    annotation_row <- NULL
    mc_annot <- generate_mc_annotation(mc_atac)
    if (all(has_name(mc_atac@metadata, c("metacell", "cell_type")))) {
        col_annot <- tibble::column_to_rownames(mc_atac@metadata[, c("metacell", "cell_type")], "metacell")
        ann_colors <- list("cell_type" = setNames(unlist(mc_atac@metadata[, "color"]), unlist(mc_atac@metadata[, "cell_type"])))
    } else {
        cts <- unique(mc_atac_clust)
        color_key <- enframe(setNames(chameleon::distinct_colors(length(cts)), cts), name = "cell_type", value = "color")
        col_annot <- enframe(setNames(
            as.numeric(names(mc_atac_clust)),
            color_key$color[match(mc_atac_clust, color_key$cell_type)]
        ),
        name = "metacell", value = "cell_type"
        )
        col_annot <- tibble::column_to_rownames(col_annot, "metacell")
        ann_colors <- list("cell_type" = deframe(color_key))
    }
    col_annot <- mc_annot[[1]]
    ann_colors <- mc_annot[[2]]
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
        annotation_row <- peak_annotation[[2]]
        ann_colors[names(peak_annotation[[1]])] <- peak_annotation[[1]]
    }
    pp <- pheatmap::pheatmap(mca_lfc[peak_clust, mc_atac_clust],
        annotation_col = subset(col_annot, select = cell_type),
        annotation_legend = FALSE,
        annotation_colors = ann_colors,
        annotation_row = annotation_row, main = main,
        color = clrs, breaks = brks, cluster_cols = FALSE, cluster_rows = FALSE, show_colnames = FALSE, show_rownames = FALSE
    )

    if (!is.null(filename)) {
        save_pheatmap(pp, filename = filename, dev = dev, ...)
    }
    return(pp)
}
