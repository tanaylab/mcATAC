% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_selection.R
\name{find_blacklist_overlaps}
\alias{find_blacklist_overlaps}
\title{Find overlaps with ENCODE blacklists}
\usage{
find_blacklist_overlaps(
  atac = NULL,
  peaks = NULL,
  genome = NULL,
  max_dist_to_blacklist_region = 0,
  blacklist_name = "ENCODE.blacklist"
)
}
\arguments{
\item{atac}{(optional) an ScPeaks or McPeaks object}

\item{peaks}{(optional) the intervals set to check}

\item{genome}{(optional, required if checking peaks directly) the genome of the peaks}

\item{max_dist_to_blacklist_region}{(optional) distance to nearest blacklist region which still qualifies for being blacklisted}

\item{blacklist_name}{name of the blacklist intervals to use (default: "ENCODE.blacklist")}
}
\value{
blacklist_overlaps - a PeakIntervals object with peaks identified as overlapping blacklisted regions
}
\description{
Identify peaks in the data which overlap (or are adjacent to?) regions blacklisted by ENCODE as having universally high DNAse HS or ChIP signal (basically mapping artifacts)
See https://doi.org/10.1038/s41598-019-45839-z for more details
}
\examples{
\dontrun{
blacklist_overlaps <- find_blacklist_overlaps(atac = my_atac_mc, max_dist_to_blacklist_region = 100)
blacklist_overlaps <- find_blacklist_overlaps(peaks = my_peak_set, genome = "mm10")
}
}
