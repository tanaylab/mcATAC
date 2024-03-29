% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/McCounts.R
\name{scc_to_mcc}
\alias{scc_to_mcc}
\alias{scc_project_on_mc}
\title{Given metacells (usually from RNA data), project ATAC counts to get a McCounts object}
\usage{
scc_to_mcc(sc_counts, cell_to_metacell, ignore_metacells = c(-1, -2))

scc_project_on_mc(sc_counts, cell_to_metacell, ignore_metacells = c(-1, -2))
}
\arguments{
\item{sc_counts}{A ScCounts object}

\item{cell_to_metacell}{a data.frame with columns \code{cell_id} and \code{metacell} containing the mapping from single cell names to metacell names, or the name of an 'h5ad' file containing this information at the 'obs' slot. In such a case, the 'obs' slot should contain
a column named \code{metacell} and the rownames should be the cell names.}

\item{ignore_metacells}{a vector of metacells to ignore. Default: \link{-1, -2} (the "outliers" metacell in the metacell2 python package).}
}
\value{
A McCounts object
}
\description{
Given cell to metacell association, summarise atac read counts to generate a McCounts object. This can
take a while - around 5 minutes using 24 cores on the PBMC dataset.
}
\examples{
\dontrun{
data(cell_to_metacell_pbmc_example)
scc_to_mcc(sc_counts, cell_to_metacell_pbmc_example)
}

}
