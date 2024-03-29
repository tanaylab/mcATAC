% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/McCounts.R
\name{scc_extract}
\alias{scc_extract}
\title{Extract summary statistics per cell for a set of intervals}
\usage{
scc_extract(scc, intervals, cells = NULL)
}
\arguments{
\item{scc}{a ScCounts object}

\item{intervals}{an intervals set. Can have a column called "peak_name" with the peak name.}

\item{cells}{a vector of cell names to include. Default: all cells.}
}
\value{
a sparse matrix where rows are the intervals, columns are the cells, and values are the counts summed. If the intervals have a column called "peak_name", the rows will be the peak names, otherwise the rows will be of the form "{chr}:{start}_{end}"
}
\description{
Extract summary statistics per cell for a set of intervals
}
\examples{
\dontrun{
scc_extract(scc, gintervals(1, 0, 100))
}

}
