
#' Download example data from 10x
#'
#' @description
#' Download the following dataset from 10x site:
#' PBMC from a healthy donor - granulocytes removed through cell
#' sorting (10k)
#' Link: https://support.10xgenomics.com/single-cell-multiome-atac-gex/datasets/1.0.0/pbmc_granulocyte_sorted_10k
#' File: Filtered feature barcode matrix MEX (DIR)
#'
#' @param dir directory to download the data to
#'
#' @examples
#' \dontrun{
#' download_pbmc_example_data()
#' }
#'
#' @export
download_pbmc_example_data <- function(dir = "pbmc_data") {
    url <- "https://support.10xgenomics.com/single-cell-multiome-atac-gex/datasets/1.0.0/pbmc_granulocyte_sorted_10k"
    file_url <- "https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.tar.gz"

    temp_file <- tempfile()
    download.file(file_url, temp_file)
    untar(temp_file)

    if (!dir.exists("filtered_feature_bc_matrix")) {
        cli_abort("Download failed. Please try to download manually from {.url {url}}")
    }
    file.rename("filtered_feature_bc_matrix", dir)

    cli_alert_success("successfully downloaded data to {.dir {dir}}")
}
