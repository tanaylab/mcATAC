% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/McCounts.R
\name{mcc_marginal}
\alias{mcc_marginal}
\alias{scc_marginal}
\title{Return the total coverage of each non-zero coordinate in an ScCounts/McCounts object}
\usage{
mcc_marginal(mc_counts, metacells = NULL)

scc_marginal(sc_counts, cells = NULL)
}
\arguments{
\item{mc_counts}{a McCounts object}

\item{metacells}{names of metacells to include. Default: all metacells.}

\item{cells}{names of cells to include. Default: all cells.}
}
\value{
an intervals set with an additional column called "cov" with the total coverage for each coordinate
}
\description{
Return the total coverage of each non-zero coordinate in an ScCounts/McCounts object
}
\examples{
\dontrun{
# cell counts
sc_counts <- scc_read("pbmc_reads")
sc_mars <- scc_marginal(sc_counts)

# metacell counts
mc_counts <- mcc_read("pbmc_reads_mc")
mc_mars <- mcc_marginal(mc_counts)
}

}
