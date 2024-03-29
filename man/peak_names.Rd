% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PeakIntervals.R
\name{peak_names}
\alias{peak_names}
\title{Return the names of PeakIntervals}
\usage{
peak_names(peaks, tad_based = TRUE, promoters = FALSE)
}
\arguments{
\item{peaks}{a PeakIntervals object or misha intervals}

\item{tad_based}{whether to name peaks based on TADs. If an intervals set named "intervs.global.tad_names" does not
exists - the function will return coordinate based names.}

\item{promoters}{are the peaks promoters? When true, the peak names (which are gene names) are returned.}
}
\value{
the field called 'peak_name' if exists in \code{peaks}, or a tad based peak name if \code{tad_based=TRUE} and the coordinates in
the form of "{chrom}:{start}-{end}" otherwise.
}
\description{
Return the names of PeakIntervals
}
\examples{
\dontrun{
misha.ext::gset_genome("hg38")
head(peak_names(atac_sc@peaks)) # return the "peak_name" field
head(peak_names(atac_sc@peaks[, -4])) # return tad based peak names
head(peak_names(atac_sc@peaks[, -4], tad_based = FALSE)) # return coordinate based peak names
}

}
