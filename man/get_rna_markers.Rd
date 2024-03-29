% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rna.R
\name{get_rna_markers}
\alias{get_rna_markers}
\title{Calculate gene markers of metacells gene expression}
\usage{
get_rna_markers(
  atac_mc,
  n_genes = 100,
  minimal_max_log_fraction = -13,
  minimal_relative_log_fraction = 2,
  fold_change_reg = 0.1,
  genes = NULL,
  rm_zeros = TRUE,
  epsilon = 0.00001
)
}
\arguments{
\item{atac_mc}{a McPeaks/McTracks object with RNA expression (using \code{add_mc_rna})}

\item{n_genes}{maximal number of genes to return. Default: 100}

\item{minimal_max_log_fraction}{take only genes with at least one value
(in log fraction units - normalized egc) above this threshold}

\item{minimal_relative_log_fraction}{take only genes with at least one value with relative
log fraction (mc_fp) above this this value}

\item{fold_change_reg}{regularization factor for the fold change calculation (fold_change would be changed to
\code{fold_change = fold_change + fold_change_reg})}

\item{genes}{list of genes to match. Default (NULL): all genes}

\item{rm_zeros}{remove genes with no RNA expression in any metacell. Default: TRUE}

\item{epsilon}{regularization factor added to the log normalized expression}
}
\description{
The function first removes genes without sufficient expression in any metacell \code{minimal_max_log_fraction}, or without sufficient fold
change over the median (\code{minimal_relative_expression}), and then ranks the fold change of genes within each metacell. The markers
are then the genes with the highest rank, up to \code{n_genes} genes.
}
\examples{
\dontrun{
get_rna_markers(atac_mc)
}

}
