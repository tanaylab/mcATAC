#' Write a sparse-matrix file from an indexed bam file
#'
#' @description This function reads a bam file and writes a sparse-matrix file, where rows are the genomic coordinates
#' and columns are the cells. It does so by:
#' \enumerate{
#' \item filtering the bam file by the given chromosome (using the "regions" paramter)
#' \item extracting the cell name using the "CB" tag added by cellranger
#' \item use awk to extract the start coordinates and the tag
#' \item use unix "sort" and "uniq" to count the number of times each coordinate appears for each cell
#' \item use an R script to convert the cell names to indices
#' \item write the sparse-matrix file (matrix-market format) and zip it
#' } \cr
#' Note (1): In order to use this function, "samtools" (>= 1.15), "awk", "sed", "sort" and "uniq" must be available at the unix command line.
#' Note (2): The function takes a while to run - around 13 minutes using 24 cores on the PBMC data.
#' 
#'
#'
#' @param bam_file name of the bam file
#' @param out_file name of the output file. A ".gz" extension will be added to the file name.
#' @param cell_names a vector with the cell names or an ScATAC object
#' @param chrom a string with a chromosome to filter. See samtools docs <http://www.htslib.org/doc/samtools-view.html> for details
#' @param genome genome name (e.g. hg19). Will be inferred from the ScATAC object if provided
#' @param min_mapq minimal mapping quality (optional)
#' @param samtools_bin path to samtools executable
#' @param samtools_opts additional options for samtools (e.g. "--subsample 0.1")
#' @param num_reads number of reads to process (optional)
#' @param verbose verbose output (optional)
#'
#' @return None
#'
#' @examples
#' \dontrun{
#' write_sparse_matrix_from_bam("pbmc_data/possorted_bam.bam", "chrom1.mm", chrom = "chr1", cell_names = atac_sc)
#' }
#'
#' @export
write_sparse_matrix_from_bam <- function(bam_file, out_file, cell_names, chrom, genome = NULL, min_mapq = NULL, samtools_bin = "/home/feshap/src/samtools-1.15.1/samtools", samtools_opts = NULL, num_reads = NULL, verbose = TRUE) {
    if (class(cell_names) == "ScATAC") {
        genome <- genome %||% cell_names@genome
        cell_names <- colnames(cell_names@mat)
    }
    gset_genome(genome)

    if (!file.exists(paste0(bam_file, ".bai"))) {
        cli_abort("Index file not found for {.file {bam_file}}. Please run 'samtools index {bam_file}'.")
    }

    cell_names_file <- paste0(out_file, ".colnames")
    write.table(cell_names, cell_names_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
    withr::defer(unlink(cell_names_file))

    intervals_df <- gintervals.all() %>%
        mutate(l = end - start) %>%
        filter(chrom == !!chrom)

    if (nrow(intervals_df) == 0) {
        cli_abort("Chromosome {.val {chrom}} not found in the genome")
    }

    mapq_filt_cmd <- ""
    if (!is.null(min_mapq)) {
        mapq_filt_cmd <- glue("-q {min_mapq}")
    }

    chrom <- chrom %||% ""
    samtools_opts <- samtools_opts %||% ""

    num_reads_cmd <- ""
    if (!is.null(num_reads)) {
        num_reads_cmd <- glue("| head -n {num_reads}")
    }

    cell_name_to_index_bin <- system.file("cell_name_to_index.R", package = "mcATAC")

    dims <- c(sum(intervals_df$l), length(cell_names))

    output <- file(out_file, open = "w")
    writeLines("%%MatrixMarket matrix coordinate real general", output)
    writeLines(paste(dims[1], dims[2], "unknown", sep = " "), output) # we do not know yet what would be the length of the file
    close(output)

    cmd <- glue("{samtools_bin} view --keep-tag CB {mapq_filt_cmd} {samtools_opts} {bam_file} {chrom} | grep CB | {awk_cmd} {num_reads_cmd} | sort | uniq -c | sed 's/^ *//g' | {cell_name_to_index_bin} {cell_names_file} >> {out_file}", awk_cmd = "awk '{print $4,substr($12, 6)}'")
    system(cmd)

    # we now replace the second header line with the actual length of the file
    num_rows <- as.numeric(system(glue("wc -l {out_file} | cut -d ' ' -f 1"), intern = TRUE)) - 2
    system(glue("sed -i '2 s/^.*$/{header}/' {out_file}", header = paste(dims[1], dims[2], as.integer(num_rows), sep = " ")))

    cli_alert("zipping {out_file} to {out_file}.gz")
    system(glue("gzip {out_file} && rm -f {out_file}"))

    if (verbose) {
        cli_alert_success("Saved sparse matrix of chromosome {.val {chrom}} with {.val {num_rows}} rows to {.file {out_file}}. Dimensions: {.file {dims[1]} x {dims[2]}}")
    }
}

#' Write single cell counts from an indexed bam file
#'
#' @param out_dir output directory.
#' @param id an identifier for the object, e.g. "pbmc".
#' @param description description of the object, e.g. "PBMC from a healthy donor - granulocytes removed through cell sorting (10k)"
#' @param num_cores number of cores to use (optional)
#'
#' @return None
#'
#' @inheritParams write_sparse_matrix_from_bam
#'
#' @examples
#' \dontrun{
#' write_sc_counts_from_bam("pbmc_data/possorted_bam.bam", "pbmc_reads", cell_names = atac_sc)
#' }
#'
#' @export
write_sc_counts_from_bam <- function(bam_file, out_dir, cell_names, genome = NULL, id = "", description = "", min_mapq = NULL, samtools_bin = "/home/feshap/src/samtools-1.15.1/samtools", samtools_opts = NULL, num_reads = NULL, verbose = TRUE, num_cores = parallel::detectCores()) {
    data_dir <- file.path(out_dir, "data")
    if (!dir.exists(data_dir)) {
        dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
        cli_alert("Created directory {.file {out_dir}}")
    }

    if (class(cell_names) == "ScATAC") {
        genome <- genome %||% cell_names@genome
        cell_names <- colnames(cell_names@mat)
    }
    gset_genome(genome)

    chroms <- intersect(gintervals.all()$chrom, bam_file_chromosomes(bam_file))

    doMC::registerDoMC(num_cores)
    plyr::l_ply(chroms, function(chrom) {
        cli_alert("Processing chromosome {.val {chrom}}")
        write_sparse_matrix_from_bam(bam_file, paste0(data_dir, "/", chrom, ".mtx"), cell_names, chrom, genome, min_mapq, samtools_bin, samtools_opts, num_reads, verbose)
        cli_alert_success("Finished processing chromosome {.val {chrom}}")
    }, .parallel = TRUE)

    counts_md <- list(
        cell_names = cell_names,
        genome = genome,
        id = id,
        description = description,
        data_dir = "data",
        chromosomes = chroms
    )

    yaml::write_yaml(counts_md, file = file.path(out_dir, "metadata.yaml"))

    cli_alert_info("Created sparse matrices for {.val {length(chroms)}} chromosomes")
    cli_alert_success("Created an ScCounts object at {.file {out_dir}}")
}


#' Get chromosomes from a bam file
#'
#' @param bam_file name of the bam file
#' @param samtools_bin path to samtools executable
#'
#' @return a vector of strings with the chromosome names
#'
#' @examples
#' \dontrun{
#' bam_file_chromosomes("pbmc_data/possorted_bam.bam")
#' }
#'
#' @export
bam_file_chromosomes <- function(bam_file, samtools_bin = "/home/feshap/src/samtools-1.15.1/samtools") {
    cmd <- glue("{samtools_bin} view -H {bam_file} | grep '^@SQ' | cut -f 2 | cut -d ':' -f 2 | cut -d ' ' -f 1 | sort | uniq", samtools_bin = samtools_bin)
    return(system(cmd, intern = TRUE))
}

#' Read an ScCounts object from a directory
#'
#' @param path path to the directory containing the object (which was created by \code{write_sc_counts_from_bam})
#' @param id an identifier for the object (optional)
#' @param description description of the object (optional)
#'
#' @return an ScCounts object
#'
#' @examples
#' \dontrun{
#' counts <- read_sc_counts("pbmc_reads")
#' }
#'
#' @export
read_sc_counts <- function(path, id = NULL, description = NULL) {
    md_file <- file.path(path, "metadata.yaml")
    if (!file.exists(md_file)) {
        cli_abort("Directory {.file {path}} does not contain a valid ScCounts object")
    }
    md <- yaml::read_yaml(md_file)
    data_dir <- file.path(path, md$data_dir)
    required_fields <- c("cell_names", "genome", "id", "description", "data_dir", "chromosomes")
    purrr::walk(required_fields, ~ {
        if (.x %!in% names(md)) {
            cli_abort("Directory {.file {path}} does not contain a valid ScCounts object (the metadata file {.file {md_file}} is missing the {.field {.x}} field.")
        }
    })

    data <- read_sc_counts_data(data_dir, md$chromosomes, md$cell_names, md$genome)

    counts <- new("ScCounts", data = data, cell_names = md$cell_names, genome = md$genome, chromosomes = md$chromosomes, id = id %||% md$id, description = description %||% md$description, path = path)
    return(counts)
}

read_sc_counts_data <- function(data_dir, chromosomes, cell_names, genome, num_cores = parallel::detectCores()) {
    gset_genome(genome)
    intervals <- gintervals.all()

    if (any(setdiff(chromosomes, intervals$chrom))) {
        missing_chroms <- setdiff(chromosomes, intervals$chrom)
        cli_abort("The following chromosomes are not present in the genome: {.val {missing_chroms}}")
    }

    chromosomes <- intersect(chromosomes, intervals$chrom)
    doMC::registerDoMC(num_cores)

    data <- plyr::llply(chromosomes, function(chrom) {
        chrom_file <- file.path(data_dir, paste0(chrom, ".mtx.gz"))
        if (!file.exists(chrom_file)) {
            cli_abort("File {.file {chrom_file}} does not exist")
        }
        m <- tgutil::fread_mm(chrom_file)
        colnames(m) <- cell_names
        return(m)
    }, .parallel = TRUE)

    names(data) <- chromosomes

    return(data)
}
