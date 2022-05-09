
#' Download example data from 10x
#'
#' @description
#' Download the following dataset from 10x site:
#' PBMC from a healthy donor - granulocytes removed through cell
#' sorting (10k)
#' \code{download_pbmc_example_data} downloads the processed data, while \code{download_pbmc_example_data_raw} downloads the raw data. \cr
#' Link: https://support.10xgenomics.com/single-cell-multiome-atac-gex/datasets/1.0.0/pbmc_granulocyte_sorted_10k \cr
#' Processed data: Filtered feature barcode matrix MEX (DIR) \cr
#' Raw data: ATAC Position-sorted alignments (BAM), ATAC Position-sorted alignments (BAM index) \cr
#'
#' @param dir directory to download the data to
#'
#' @examples
#' \dontrun{
#' download_pbmc_example_data()
#' download_pbmc_example_data_raw()
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

    cli_alert_success("successfully downloaded data to {.file {dir}}")
}

#' @rdname download_pbmc_example_data
#' @export
download_pbmc_example_data_raw <- function(dir = "pbmc_data") {
    url <- "https://support.10xgenomics.com/single-cell-multiome-atac-gex/datasets/1.0.0/pbmc_granulocyte_sorted_10k"
    bam_url <- "https://cg.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_possorted_bam.bam"
    index_url <- "https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_possorted_bam.bam.bai"

    if (!dir.exists(dir)) {
        dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }

    withr::with_options(list(timeout = 1e5), {
        download.file(bam_url, file.path(dir, "possorted_bam.bam"))
        download.file(index_url, file.path(dir, "possorted_bam.bam.bai"))
    })

    cli_alert_success("successfully downloaded raw data to {.file {dir}}")
}
