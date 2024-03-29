% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rna.R
\name{get_genes_atac_fp}
\alias{get_genes_atac_fp}
\title{Return a relative fold change matrix of ATAC peaks for a list of genes}
\usage{
get_genes_atac_fp(
  atac_mc,
  genes = NULL,
  metacells = NULL,
  rm_zeros = TRUE,
  spearman = TRUE,
  pairwise.complete.obs = TRUE
)
}
\arguments{
\item{atac_mc}{a McPeaks/McTracks object}

\item{genes}{a list of genes.}

\item{rm_zeros}{remove genes with no RNA expression in any metacell. Default: TRUE}

\item{spearman}{if 'TRUE' Spearman correlation is computed, otherwise
Pearson}

\item{pairwise.complete.obs}{see below}

\item{metacell}{select only a subset of the metacells.}
}
\description{
This function returns a relative fold change matrix of ATAC peaks for a list of genes matched using \code{rna_atac_cor_knn}.
}
\examples{
\dontrun{
marker_genes <- get_rna_markers(atac_mc)
atac_fp <- get_genes_atac_fp(atac_mc, genes = marker_genes)
rna_fp <- get_rna_marker_matrix(atac_mc, genes = marker_genes)
}

}
