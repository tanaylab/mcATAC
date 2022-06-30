#' ATACTracks objects
#'
#'
#' @description ATACTracks is a shallow object holding ATAC data over cells/metacells. Minimally it should include a vector of track names which hold the names of the tracks with the ATAC data.
#' McTracks extend the ATACTracks object by adding metadata and additional slots.
#'
#' @slot id an identifier for the object, e.g. "pbmc".
#' @slot description description of the object, e.g. "PBMC from a healthy donor - granulocytes removed through cell sorting (10k)"
#' @slot genome genome assembly of the peaks. e.g. "hg38", "hg19", "mm9", "mm10"
#' @slot metadata data frame with a column called 'metacell' and additional metacell annotations. The constructor can also include or the name of a delimited file which contains such annotations.
#'
#' @exportClass ATACPeaks
ATACTracks <- setClass(
    "ATACTracks",
    slots = c(
        tracks = "vector"
    ),
    contains = c("ATAC", "VIRTUAL")
)
