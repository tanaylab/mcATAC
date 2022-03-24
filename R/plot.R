#' Plot a scatter of gene expression vs an atac profile
#'
#' @param mc_atac a McATAC object
#' @param mc_rna a \code{metacell1} 'mc' object or a \code{metacells} metacell UMI matrix
#' @param gene name of the gene to plot
#' @param peak name of the peak to plot. If NULL - the promoter of \code{gene} would be shown
#' @param max_dist_to_promoter_peak how far from \code{gene}'s TSS to search for a promoter-proximal peak
#' @param eps_q quantile of mean expression distribution to add as regularization when calculating log expression. Promoter ATAC signal has been observed empirically to often be linearly correlated with log of gene expression.
#' @return a ggplot object with the scatter plot
#'
#' @examples
#' \dontrun{
#' p1 <- plot_atac_rna(mc_atac = pbmc_atac_mc, mc_rna = pbmc_rna_mc, gene = "CD4")
#'
#' ## Plot gene vs some peak of interest; here - expression of CD4 vs. promoter ATAC signal of CD8
#' tss <- gintervals.load("intervs.global.tss")
#' cd8_tss <- tss[which(tss$geneSymbol == "CD8")[[1]], ]
#' peaks_near_cd8_tss <- misha::gintervals.neighbors(pbmc_atac_mc@peaks, cd8_tss, mindist = -2e+3, maxdist = 2e+3)
#' nearest_peak <- peaks_near_cd8_tss[which.min(peaks_near_cd8_tss$dist), ]
#' p2 <- plot_atac_rna(mc_atac = pbmc_atac_mc, mc_rna = pbmc_rna_mc, gene = "CD4", peak = nearest_peak[, "peak_name"])
#' }
#' @export
plot_atac_rna <- function(mc_atac, mc_rna, gene, peak = NULL, max_dist_to_promoter_peak = 5e+2, eps_q = 0.05) {

    ### Placeholder for checking class of mc_rna
    if (class(mc) == "tgMCCov") {
        eps_rna <- quantile(rowMeans(mc_rna@e_gc), 0.05)
        rv <- log2(mc_rna@e_gc[gene, ] + eps_rna)
    }
    # else if (mc_rna is a metacells metacell UMI matrix) {
    #   eps_rna <- quantile(rowMeans(mc2_mc@UMI_mat), 0.05)
    #   rv <- log2(mc2_mc@UMI_mat[gene, ] + eps_rna)
    # }
    # else {cli_abort('mc_rna is not a valid metacell1 or metacells metacell matrix')}
    ### end placeholder

    if (is.null(peak)) {
        tss <- gintervals.load("tss")
        tss_gene <- tss[tss$geneSymbol == gene, c("chrom", "start", "end", "geneSymbol")]
        if (nrow(tss_gene) == 0) {
            cli_abort("Peak not supplied and gene TSS location not found")
        }
        nei_peaks_tss <- gintervals.neighbors(tss_gene[, 1:3], as.data.frame(mc_atac@peaks),
            mindist = -max_dist_to_promoter_peak,
            maxdist = max_dist_to_promoter_peak
        )
        if (nrow(nei_peaks_tss) == 0) {
            cli_abort("Peak not supplied and no relevant promoter peak found for gene. Check if your ATAC matrix should have a peak for this gene's promoter.")
        }
        prom_peak <- nei_peaks_tss[which.min(nei_peaks_tss$dist), 4:6]
        if (length(unlist(stringr::str_split(rownames(mc_atac@mat_[[1]]), "-"))) == 3) {
            peak <- paste0(unlist(prom_peak), collapse = "-")
        } else {
            peak <- paste0(prom_peak[[1]], ":", prom_peak[[2]], "-", prom_peak[[3]])
        }
        av <- mc_atac@mat[peak, ]
    } else {
        av <- mc_atac[peak, ]
    }

    if (any(grepl("color", colnames(mc_atac@metadata)))) {
        clrs <- unique(mc_atac@metadata[, c("cell_type", "color")])
    } else {
        clrs <- as.data.frame(list("cell_type" = "all_mc", "color" = "black"))
    }
    df <- as.data.frame(list("atac" = av, "rna" = rv, "cell_type" = unlist(mc_atac@metadata[, "cell_type"])))
    gg <- ggplot2::ggplot(df, aes(x = rna, y = atac, color = cell_type)) +
        ggplot2::geom_point() +
        labs(title = glue::glue("{gene} - ATAC of {peak} vs. RNA"), x = "log2 e_gc", y = "ATAC (not normalized)") +
        scale_color_manual(values = setNames(unlist(clrs[, "color"]), unlist(clrs[, "cell_type"])))
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
#' @param gene_field name of a field in \code{mc_atac@peaks} which contains the gene names. If NULL - the peaks would be
#' transformed to promoter peaks and the gene names would be taken from the promoter gene names.
#'
#' @export
plot_atac_rna_cor <- function(mc_atac, rna_mat) {

}

#' Plot normalized accessibility of peaks over metacells, ordered by clustering
#'
#' @inheritDotParams save_pheatmap_png
#' @param mc_atac McATAC object
#' @param mc_atac_clust a clustering of metacells (e.g. output of {.code gen_atac_mc_clust})
#' @param filename file name of saved image
#' @param peak_clust (optional) output of \code{gen_atac_peak_clust}
#' @examples
#' \dontrun{
#'
#' }
#' @export
plot_atac_peak_map <- function(mc_atac, mc_atac_clust, filename, peak_clust, eps_q = 0.1, ...) {
    # the central heat map showing normalized accessibility of peaks over metacells, ordered by clustering
    clrs = colorRampPalette(c('blue4', 'white', 'red4'))(100)
    col_annot <- tibble::column_to_rownames(mc_atac@metadata[,c('metacell', 'cell_type')], 'metacell')
    ann_colors <- list('cell_type' = setNames(unlist(mc_atac@metadata[,'color']), unlist(mc_atac@metadata[,'cell_type'])))
    eps = quantile(rowMeans(mc_atac@mat), eps_q)
    if (is.null(mc_atac_clust)) {
        cli_abort('Must specify clustering of metacells (e.g. using {.code gen_atac_mc_clust})')
    }
    if (is.null(peak_clust)) {
        cli_abort('Must specify clustering of peaks (e.g. using {.code gen_atac_peak_clust})')
    }
    mca_lfc = t(apply(mc_atac@mat, 1, function(x) log2((x + eps)/median(x + eps))))
    brks = c(seq(min(mca_lfc), 0, l=50),
        seq(0.01*(max(mca_lfc) - min(mca_lfc)), max(mca_lfc), l=51))
    colnames(mca_lfc) = 1:ncol(mca_lfc)
    pp = pheatmap::pheatmap(mca_lfc[peak_clust,mc_atac_clust], annotation_col = subset(col_annot, select = cell_type), 
                        annotation_legend = FALSE,
                        annotation_colors = ann_colors['cell_type'], 
                        color = clrs, breaks = brks, cluster_cols = F, cluster_rows = F, show_colnames = F, show_rownames = F)
    if (!dir.exists('./figs')) {dir.create('./figs')}
    save_pheatmap_png(pp, glue::glue('./figs/{filename}'), ...)
    return(pp)
}