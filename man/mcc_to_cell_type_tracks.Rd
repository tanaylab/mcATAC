% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/McCounts.R
\name{mcc_to_cell_type_tracks}
\alias{mcc_to_cell_type_tracks}
\title{Convert MC counts to cell type tracks}
\usage{
mcc_to_cell_type_tracks(
  mc_counts,
  track_prefix,
  overwrite = FALSE,
  resolution = 20,
  window_size = NULL
)
}
\arguments{
\item{mc_counts}{A data frame containing MC counts.}

\item{track_prefix}{A character string specifying the prefix for the output track names.}

\item{overwrite}{A logical value indicating whether to overwrite existing tracks with the same name. Default is FALSE.}

\item{resolution}{An integer specifying the resolution for the tracks. Default is 20.}

\item{window_size}{An optional integer specifying the window size for smoothing the tracks. Default is NULL.}
}
\value{
A list of cell type tracks.
}
\description{
This function takes in MC counts and converts them into cell type tracks.
}
\examples{
\dontrun{
mcc_to_cell_type_tracks(mc_counts, "pbmc_cell_types")
}

}
