% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ATAC.R
\name{mc_order}
\alias{mc_order}
\title{Order the metacells in an MC object}
\usage{
mc_order(obj, order)
}
\arguments{
\item{obj}{an McPeaks or McTracks object}

\item{order}{a vector of integers with the order of the metacells}
}
\description{
Order the metacells in an MC object
}
\examples{
\dontrun{
atac_mc <- mc_order(atac_mc, rev(1:length(atac_mc@metacells)))
}

}
