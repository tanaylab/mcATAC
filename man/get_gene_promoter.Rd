% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-intervals.R
\name{get_gene_promoter}
\alias{get_gene_promoter}
\title{Get an intervals set of a gene promoter.}
\usage{
get_gene_promoter(gene, upstream = 500, downstream = 50, unify = FALSE)
}
\arguments{
\item{upstream}{bp upstream to TSS}

\item{downstream}{bp downstread from tss}

\item{unify}{unify the intervals if there are multiple promoters. The start coordinate of the intervals will be the minimum of the starts and the end coordinate of the intervals will be the maximum of the ends.}

\item{name}{The name of the gene ("name2" in UCSC files).}
}
\value{
an intervals set with the promoter of the gene. If \code{unify} is TRUE - the intervals will be unified to a single row.
}
\description{
Get an intervals set of a gene promoter.
}
\examples{
\dontrun{
get_gene_promoter("GZMK")
}

}
