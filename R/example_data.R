
#' Download example data from 10x
#'
#' @description
#' Download the following dataset from 10x site:
#' PBMC from a healthy donor - granulocytes removed through cell
#' sorting (10k)
#' When \code{fragments=TRUE} the fragment file and its index are downloaded as well.
#' Link: https://support.10xgenomics.com/single-cell-multiome-atac-gex/datasets/1.0.0/pbmc_granulocyte_sorted_10k \cr
#' Processed data: Filtered feature barcode matrix MEX (DIR) \cr
#' Fragments: ATAC Per fragment information file (TSV.GZ) \cr
#'
#' @param dir directory to download the data to
#' @param fragments download fragments file (and its index)
#'
#' @examples
#' \dontrun{
#' download_pbmc_example_data()
#' download_pbmc_example_data(fragments = TRUE)
#' }
#'
#' @export
download_pbmc_example_data <- function(dir = "pbmc_data", fragments = FALSE) {
    url <- "https://support.10xgenomics.com/single-cell-multiome-atac-gex/datasets/1.0.0/pbmc_granulocyte_sorted_10k"
    file_url <- "https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.tar.gz"

    temp_file <- tempfile()
    download.file(file_url, temp_file)
    untar(temp_file)

    if (!dir.exists("filtered_feature_bc_matrix")) {
        cli_abort("Download failed. Please try to download manually from {.url {url}}")
    }
    file.rename("filtered_feature_bc_matrix", dir)
    cli_alert_info("downloaded processed matrix")

    if (fragments) {
        fragments_url <- "https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
        fragment_index_url <- "https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi"
        download.file(fragments_url, file.path(dir, "fragments.tsv.gz"))
        download.file(fragment_index_url, file.path(dir, "fragments.tsv.gz.tbi"))
        cli_alert_info("downloaded fragments")
    }

    cli_alert_success("successfully downloaded data to {.file {dir}}")
}
