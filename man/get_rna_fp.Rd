% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rna.R
\name{get_rna_fp}
\alias{get_rna_fp}
\title{Get enrichment of normalized counts of gene expression over the median}
\usage{
get_rna_fp(atac_mc, genes = NULL, rm_zeros = TRUE, epsilon = 0.00001)
}
\arguments{
\item{atac_mc}{a McPeaks/McTracks object with RNA expression (using \code{add_mc_rna})}

\item{genes}{list of genes to match. Default (NULL): all genes}

\item{rm_zeros}{remove genes with no RNA expression in any metacell. Default: TRUE}

\item{epsilon}{regularization factor added to the log normalized expression}
}
\value{
a matrix with normalized counts of gene expression for each gene (rows) and metacell (columns)
}
\description{
The function first transforms the UMI matrix to fractions per metacell, and then calculates the enrichment of each gene
over the median (dividing the gene expression by the median).
}
\examples{
\dontrun{
get_rna_fp(atac_mc)
get_rna_fp(atac_mc, "GNLY")
get_rna_fp(atac_mc, "GNLY", epsilon = NULL)
}

}
