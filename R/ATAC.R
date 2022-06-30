setClassUnion("any_matrix", c("sparseMatrix", "matrix"))
setClassUnion("data.frame_or_null", c("data.frame", "NULL"))

setOldClass("PeakIntervals")

#' ATAC objects
#'
#'
#' @description ATAC is a shallow object holding ATAC data over cells/metacells. Minimally it needs to have an id, a description and a genome id. \cr
#' ATACPeaks classes hold the data in the format of peaks and are used mostly for "trans" analysis. \cr
#' ATACCounts classes hold the raw data in the format of sparse matrices - it can then be transformed into misha tracks. \cr
#' ATACTracks classes hold pointers to tracks that hold the data at higher resolution (e.g. 20 bp), and can be used for "cis" analysis.
#'
#' @slot id an identifier for the object, e.g. "pbmc".
#' @slot description description of the object, e.g. "PBMC from a healthy donor - granulocytes removed through cell sorting (10k)"
#' @slot genome genome assembly of the peaks. e.g. "hg38", "hg19", "mm9", "mm10"
#' @slot metadata data frame with a column called 'metacell' and additional metacell annotations for McATAC, or 'cell_id' and per-cell annotations for ScATAC. The constructor can also include or the name of a delimited file which contains such annotations.
#' @slot path original path from which the object was loaded (optional)
#'
#' @exportClass ATAC
ATAC <- setClass(
    "ATAC",
    slots = c(
        id = "character",
        description = "character",
        genome = "character",
        metadata = "data.frame_or_null",
        path = "character"
    ),
    contains = "VIRTUAL"
)
