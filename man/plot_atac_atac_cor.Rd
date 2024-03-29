% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_atac_atac_cor}
\alias{plot_atac_atac_cor}
\title{Plot a correlation matrix of ATAC metacells}
\usage{
plot_atac_atac_cor(atac_mc, sp_f = TRUE)
}
\arguments{
\item{atac_mc}{McPeaks object}

\item{sp_f}{whether to use Spearman correlation (default) or Pearson}
}
\value{
p a pheatmap of ATAC metacell correlations
}
\description{
Plot a correlation matrix of ATAC metacells
}
\examples{
\dontrun{
p1 <- plot_atac_atac_cor(my_atac_mc)
p2 <- plot_atac_atac_cor(my_atac_mc, sp_f = F)
}
}
