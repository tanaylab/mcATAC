% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/McCounts.R
\name{mcc_read}
\alias{mcc_read}
\title{Read an McCounts object from a directory}
\usage{
mcc_read(path, id = NULL, description = NULL, verbose = TRUE)
}
\arguments{
\item{path}{path to the directory containing the object (which was created by \code{scc_from_bam})}

\item{id}{an identifier for the object (optional)}

\item{description}{description of the object (optional)}

\item{verbose}{print success message (default: TRUE)}
}
\value{
A McCounts object
}
\description{
Read an McCounts object from a directory
}
\examples{
\dontrun{
mc_counts <- mcc_read("pbmc_reads_mc")
}

}
