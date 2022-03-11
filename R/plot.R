#' Plot a scatter of gene expression vs an atac profile
#'
#' @param mc_atac a McATAC object
#' @param mc_rna a metacell 'mc' object
#' @param gene name of the gene to plot
#' @param peak name of the peak to plot. If NULL - the promoter of \code{gene} would be shown
#' @param max_dist_to_promoter_peak how far from \code{gene}'s TSS to search for a promoter-proximal peak
#' @param eps_q quantile of mean expression distribution to add as regularization when calculating log expression. Promoter ATAC signal has been observed empirically to often be linearly correlated with log of gene expression.
#' @return a ggplot object with the scatter plot
#'
#' @export
plot_atac_rna <- function(mc_atac, mc_rna, gene, peak = NULL, max_dist_to_promoter_peak = 5e+2, eps_q = 0.05) {
    if (is.null(peak)) {
        tss = gintervals.load('tss')
        tss_gene = tss[tss$geneSymbol == gene,c('chrom', 'start', 'end', 'geneSymbol')]
        if (nrow(tss_gene) == 0) {
            cli_abort('Peak not supplied and gene TSS location not found')
        }
        nei_peaks_tss = gintervals.neighbors(tss_gene[,1:3], as.data.frame(mc_atac$peaks), 
            mindist = -max_dist_to_promoter_peak, 
            maxdist = max_dist_to_promoter_peak)
        if (nrow(nei_peaks_tss) == 0) {
            cli_abort("Peak not supplied and no relevant promoter peak found for gene. Check if your ATAC matrix should have a peak for this gene's promoter.")
        }
        prom_peak = nei_peaks_tss[which.min(nei_peaks_tss$dist), 4:6]
        if (length(unlist(stringr::str_split(rownames(mc_atac$mat_[[1]]), '-'))) == 3) {
            peak = paste0(unlist(prom_peak), collapse = '-')
        }
        else {peak = paste0(prom_peak[[1]], ':', prom_peak[[2]], '-', prom_peak[[3]])}
        av = mc_atac$mat[peak,]
    }
    else {av = mc_atac[peak,]}
    eps_rna = quantile(apply(mc_rna@e_gc, 1, mean), eps_q)
    rv = log2(mc_rna@e_gc[gene,] + eps_rna)
    if (any(grepl('color', colnames(mc_atac$metadata)))) {clrs = unique(mc_atac$metadata[,c('cell_type', 'color')])}
    else {clrs = as.data.frame(list('cell_type' = 'all_mc', 'color' = 'black'))}
    df = as.data.frame(list('atac' = av, 'rna' = rv, 'cell_type' = unlist(mc_atac$metadata[,'cell_type'])))
    gg = ggplot2::ggplot(df, aes(x = rna, y = atac, color = cell_type)) + ggplot2::geom_point() + 
            labs(title=glue::glue('{gene} - ATAC of {peak} vs. RNA'), x = 'log2 e_gc', y = 'ATAC (not normalized)') + 
            scale_color_manual(values = tibble::deframe(clrs[,c('color', 'cell_type')]))
    return(gg)
}

#' Plot a correlation matrix of ATAC metacells
#'
#' @param mc_atac McATAC object
#' @param sp_f whether to use Spearman correlation (default) or Pearson
#'
#' @return p a pheatmap of ATAC metacell correlations
#'
#' @export
plot_atac_atac_cor <- function(mc_atac, sp_f = TRUE) {
    if (all(!grepl('cell_type', colnames(mc_atac$metadata)))) {
        if (all(!grepl('cluster_k_', colnames(mc_atac$peaks)))) {
            k = round(ncol(mc_atac$mat)/10)
            message(glue::glue('There is no "cell_type" or "cluster" field in metadata, clustering with k == {k}'))
            mc_atac = gen_atac_mc_clust(mc_atac,k=k,use_prior_annot=F)
            clust_vec = unlist(mc_atac$metadata[,paste0('cluster_k_',k)])
        }
        else  {
            clust_vec = unlist(mc_atac$metadata[,grep('cluster_k_',colnames(mc_atac$metadata))[[1]]])
        }
    }
    else {clust_vec = unlist(mc_atac$metadata[,'cell_type'])}
    col_annot = as.data.frame(list('cluster' = clust_vec))
    rownames(col_annot) = 1:nrow(mc_atac$metadata)
    if (any(grepl('color', colnames(mc_atac$metadata)))) {
        col_key = unique(mc_atac$metadata[,c('cell_type', 'color')])
        ann_colors = list('cluster' = setNames(col_key$color, col_key$cell_type))
    }
    else {
        ann_colors = list('cluster' = setNames(sample(grep('white|gray|grey', colors(), v=T, inv=T), 
                    length(unique(clust_vec))), unique(clust_vec)))
    }
    cor_mat = tgs_cor(mc_atac$mat, spearman = sp_f, pairwise.complete.obs = TRUE)
    p = pheatmap::pheatmap(cor_mat[order(clust_vec), order(clust_vec)], 
                            # cluster_rows = F, cluster_cols = F, 
                            show_rownames = F, show_colnames = F, 
                            annotation_col = col_annot, annotation_colors = ann_colors)
    return(p)
}

#' Plot a cross-correlation matrix between RNA metacells and ATAC metacell scores
#'
#' @description use peak gene annotations to match between RNA metacells and ATAC metacells and plot
#' a cross-correlation matrix.
#'
#' @param mc_atac McATAC object
#' @param rna_mat an RNA metacell count matrix, where metacells are in columns and genes are in rows
#' @param gene_field name of a field in \code{mc_atac$peaks} which contains the gene names. If NULL - the peaks would be
#' transformed to promoter peaks and the gene names would be taken from the promoter gene names.
#'
#' @export
plot_atac_rna_cor <- function(mc_atac, rna_mat) {

}

#' Plot normalized accessibility of peaks over metacells, ordered by clustering
#'
#' @param mc_atac McATAC object
#' @param mc_atac_clust output of \code{gen_atac_mc_clust}
#' @param peak_clust output of \code{gen_atac_peak_clust}
#'
#' @export
plot_atac_peak_map <- function(mc_atac, mc_atac_clust, peak_clust) {
    # the central heat map showing normalized accessibility of peaks over metacells, ordered by clustering
}