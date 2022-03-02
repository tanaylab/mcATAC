#' Plot a scatter of gene expression vs an atac profile
#'
#' @param mc_atac a McATAC object
#' @param gene name of the gene top plot
#' @param peak name of the peak to plot. If NULL - the promoter of \code{gene} would be shown
#'
#' @return a ggplot object with the scatter plot
#'
#' @export
plot_atac_rna <- function(mc_atac, gene, peak = NULL) {

}

#' Plot a correlation matrix of ATAC metacells
#'
#' @param mc_atac McATAC object
#'
#' @export
plot_atac_atac_cor <- function(mc_atac) {

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
