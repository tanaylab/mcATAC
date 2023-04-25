#' Write a sparse-matrix file from an indexed bam file
#'
#' @description This function reads a bam file and writes a sparse-matrix file, where rows are the genomic coordinates
#' and columns are the cells. It does so by:
#' \enumerate{
#' \item filtering the bam file by the given region (e.g. chromosom, using the "regions" paramter)
#' \item extracting the cell name using the "CB" tag added by cellranger, while excluding reads which were marked as PCR duplicates (1024)
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
#' @param cell_names a vector with the cell names or an ATAC object
#' @param region intervals set to filter. See samtools docs <http://www.htslib.org/doc/samtools-view.html> for details
#' @param genome genome name (e.g. hg19). Will be inferred from the ScPeaks object if provided
#' @param min_mapq minimal mapping quality (optional)
#' @param samtools_bin path to samtools executable
#' @param samtools_opts additional options for samtools (e.g. "--subsample 0.1")
#' @param num_reads number of reads (within the \code{region}) to process (optional).
#' @param verbose verbose output (optional)
#' @param overwrite overwrite existing files (optional)
#'
#' @return None
#'
#' @examples
#' \dontrun{
#' write_sparse_matrix_from_bam("pbmc_data/possorted_bam.bam", "chrom1.mm", region = gintervals.all()[1, ], cell_names = atac_sc)
#' }
#'
#' @export
write_sparse_matrix_from_bam <- function(bam_file, out_file, cell_names, region, genome = NULL, min_mapq = NULL, samtools_bin = "samtools", samtools_opts = NULL, num_reads = NULL, verbose = TRUE, overwrite = FALSE) {
    withr::local_options(list(scipen = 1e5))
    if (class(cell_names) == "ScPeaks") {
        genome <- genome %||% cell_names@genome
        cell_names <- colnames(cell_names@mat)
    }
    gset_genome(genome)

    if (!file.exists(paste0(bam_file, ".bai"))) {
        cli_abort("Index file not found for {.file {bam_file}}. Please run 'samtools index {bam_file}'.")
    }

    overwrite_file(out_file, overwrite)
    overwrite_file(paste0(out_file, "gz"), overwrite)

    cell_names_file <- paste0(out_file, ".colnames")
    write.table(cell_names, cell_names_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
    withr::defer(unlink(cell_names_file))

    if (nrow(region) != 1 || !all(c("chrom", "start", "end") %in% colnames(region))) {
        cli_abort("Invalid region parameter. Must be a single row with columns 'chrom', 'start' and 'end'.")
    }

    if (region$chrom %!in% gintervals.all()$chrom) {
        cli_abort("Chromosome {.val {chrom}} not found in the genome")
    }

    fixed_region <- gintervals.force_range(as.data.frame(region))
    if (fixed_region$start != region$start || fixed_region$end != region$end) {
        cli_warn("Region {.val {region}} was adjusted to {.val {fixed_region}}")
    }

    region_str <- paste0(fixed_region$chrom, ":", fixed_region$start, "-", fixed_region$end)

    mapq_filt_cmd <- ""
    if (!is.null(min_mapq)) {
        mapq_filt_cmd <- glue("-q {min_mapq}")
    }

    samtools_opts <- samtools_opts %||% ""

    num_reads_cmd <- ""
    if (!is.null(num_reads)) {
        num_reads_cmd <- glue("| head -n {num_reads}")
    }

    dims <- c(abs(fixed_region$end - fixed_region$start), length(cell_names))

    output <- file(out_file, open = "w")
    writeLines("%%MatrixMarket matrix coordinate real general", output)
    writeLines(paste(dims[1], dims[2], "unknown", sep = " "), output) # we do not know yet what would be the length of the file
    close(output)

    cmd <- glue("{samtools_bin} view --keep-tag CB --exclude-flags 1024 {mapq_filt_cmd} {samtools_opts} {bam_file} {region_str} | grep CB | {awk_cmd} {num_reads_cmd} | sort | uniq -c | sed 's/^ *//g' | {cell_name_to_index_bin} {cell_names_file} {fixed_region$start} {fixed_region$end} >> {out_file}", awk_cmd = "awk '{print $4,substr($12, 6)}'", cell_name_to_index_bin = system.file("exec", "cell_name_to_index.R", package = "mcATAC"))
    system(cmd)

    # we now replace the second header line with the actual length of the file
    num_rows <- as.numeric(system(glue("wc -l {out_file} | cut -d ' ' -f 1"), intern = TRUE)) - 2
    system(glue("sed -i '2 s/^.*$/{header}/' {out_file}", header = paste(dims[1], dims[2], as.integer(num_rows), sep = " ")))

    if (verbose) {
        cli_alert("zipping {out_file} to {out_file}.gz")
    }
    system(glue("gzip {out_file} && rm -f {out_file}"))

    if (verbose) {
        cli_alert_success("Saved sparse matrix of {.val {region_str}} with {.val {num_rows}} rows to {.file {out_file}}. Dimensions: {.file {dims[1]} x {dims[2]}}")
    }
}

#' Write a sparse-matrix file from an 10x fragments file for a specific region
#'
#' @description This function reads a "fragments.tsv.gz" which is an output of the 10x pipeline ("Barcoded and aligned fragment file"), and writes a sparse-matrix file, where rows are the genomic coordinates and columns are the cells.
#'
#' @param fragments_file path to fragments file
#' @param out_file name of the output file. A ".gz" extension will be added to the file name.
#' @param cell_names a vector with the cell names or an ScPeaks object
#' @param region intervals set to filter.
#' @param overwrite overwrite existing file (optional)
#' @param verbose verbose output (optional)
#' @param tabix_bin path to the tabix binary (optional)
#' @param use_tabix use tabix to speed up region filtering (optional)
#'
#' @return None
#'
#' @examples
#' \dontrun{
#' write_sparse_matrix_from_fragments("pbmc_data/fragments.tsv.gz", "matrix.mm", cell_names = atac_sc, gintervals(1, 0, 1e5))
#' }
#'
#' @noRd
write_sparse_matrix_from_fragments <- function(fragments_file, out_file, cell_names, region, num_reads = NULL, genome = NULL, overwrite = FALSE, verbose = FALSE, tabix_bin = "tabix", use_tabix = FALSE) {
    withr::with_options(list(scipen = 1e5), {
        if (class(cell_names) == "ScPeaks") {
            genome <- genome %||% cell_names@genome
            cell_names <- colnames(cell_names@mat)
        }

        gset_genome(genome)

        overwrite_file(out_file, overwrite)
        overwrite_file(paste0(out_file, "gz"), overwrite)

        cell_names_file <- paste0(out_file, ".colnames")
        write.table(cell_names, cell_names_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
        withr::defer(unlink(cell_names_file))

        if (nrow(region) != 1 || !all(c("chrom", "start", "end") %in% colnames(region))) {
            cli_abort("Invalid region parameter. Must be a single row with columns 'chrom', 'start' and 'end'.")
        }

        if (region$chrom %!in% gintervals.all()$chrom) {
            cli_abort("Chromosome {.val {chrom}} not found in the genome")
        }

        fixed_region <- gintervals.force_range(as.data.frame(region))
        if (fixed_region$start != region$start || fixed_region$end != region$end) {
            cli_warn("Region {.val {region}} was adjusted to {.val {fixed_region}}. This usually means we constructing a matrix for the wrong assembly. The current genome assembly is {.val {genome}}.")
        }

        dims <- c(abs(fixed_region$end - fixed_region$start), length(cell_names))

        output <- file(out_file, open = "w")
        writeLines("%%MatrixMarket matrix coordinate real general", output)
        writeLines(paste(dims[1], dims[2], "unknown", sep = " "), output) # we do not know yet what would be the length of the file
        close(output)

        # select only a small number of reads (for debugging purposes)
        num_reads_cmd <- ""
        if (!is.null(num_reads)) {
            num_reads_cmd <- glue("head -n {num_reads} | ")
        }

        if (grepl("gz$", tools::file_ext(fragments_file))) {
            cat_cmd <- "zcat"
        } else {
            cat_cmd <- "cat"
            use_tabix <- FALSE
        }

        if (use_tabix) {
            filter_cmd <- glue("{tabix_bin} {fragments_file} {fixed_region$chrom}:{fixed_region$start}-{fixed_region$end} | {num_reads_cmd} awk '{{print $2,$4; print $3,$4}}'")
        } else {
            filter_cmd <- glue("{num_reads_cmd} awk '$1 == \"{region$chrom}\" && $2 >= {region$start} && $3 < {region$end} {{print $2,$4; print $3,$4}}'")
        }

        cat_cmd <- glue("{cat_cmd} {fragments_file} | {filter_cmd} ")

        cmd <- glue("{cat_cmd} | sort | uniq -c | sed 's/^ *//g' | {cell_name_to_index_bin} {cell_names_file} {fixed_region$start} {fixed_region$end} >> {out_file}", cell_name_to_index_bin = system.file("exec", "cell_name_to_index.R", package = "mcATAC"))
        system(cmd)

        # we now replace the second header line with the actual length of the file
        num_rows <- as.numeric(system(glue("wc -l {out_file} | cut -d ' ' -f 1"), intern = TRUE)) - 2
        system(glue("sed -i '2 s/^.*$/{header}/' {out_file}", header = paste(dims[1], dims[2], as.integer(num_rows), sep = " ")))

        if (verbose) {
            cli_alert("zipping {out_file} to {out_file}.gz")
        }
        system(glue("gzip {out_file} && rm -f {out_file}"))

        if (verbose) {
            cli_alert_success("Saved sparse matrix of {.val {region_str}} with {.val {num_rows}} rows to {.file {out_file}}. Dimensions: {.file {dims[1]} x {dims[2]}}")
        }
    })
}

#' Write single cell counts from 10x 'outs' directory
#'
#' @description The function reads the "outs/atac_fragments.tsv.gz" file, and "outs/filtered_feature_bc_matrix/barcodes.tsv.gz" file and
#' writes an ScCounts object to "out_dir".
#'
#' @param input_dir path to the 10x 'outs' directory
#'
#' @return None
#'
#' @examples
#' \dontrun{
#' write_sc_counts_from_10x("outs", "pbmc_reads", cell_suffix = "batch1", genome = "hg38")
#' }
#'
#' @inheritParams write_sc_counts_from_fragments
#'
#' @export
write_sc_counts_from_10x <- function(input_dir, out_dir, cell_suffix = NULL, genome = NULL, bin_size = 5e7, overwrite = FALSE, id = "", description = "", verbose = FALSE, tabix_bin = "tabix", chromosomes = NULL) {
    fragments_file <- file.path(input_dir, "atac_fragments.tsv.gz")
    if (!file.exists(fragments_file)) {
        cli_abort("File {.file {fragments_file}} does not exist.")
    }
    barcodes_file <- file.path(input_dir, "filtered_feature_bc_matrix", "barcodes.tsv.gz")
    if (!file.exists(barcodes_file)) {
        cli_abort("File {.file {barcodes_file}} does not exist.")
    }
    barcodes <- tgutil::fread(barcodes_file, header = FALSE, stringsAsFactors = FALSE)[, 1]
    write_sc_counts_from_fragments(fragments_file, out_dir, barcodes, cell_suffix = cell_suffix, genome = genome, bin_size = bin_size, overwrite = overwrite, id = id, description = description, verbose = verbose, tabix_bin = tabix_bin, chromosomes = chromosomes)
}

#' Write single cell counts from a 10x fragments file
#'
#' @description This function writes an ScCounts object from a "fragments.tsv.gz" which is an output of the 10x pipeline ("Barcoded and aligned fragment file").
#'
#' @param fragments_file path to fragments file. Note that in order to use tabix, the file must be compressed with gzip and have a ".gz" extension.
#' @param out_dir output directory.
#' @param cell_names a vector with the cell names or an ScPeaks object
#' @param cell_suffix a suffix to add to the cell names (optional). This is useful when you want to merge multiple batches of cells.
#' The suffix would be added in the following way: `paste0(cell_name, "-", cell_suffix)`.
#' @param id an identifier for the object, e.g. "pbmc".
#' @param description description of the object, e.g. "PBMC from a healthy donor - granulocytes removed through cell sorting (10k)"
#' @param bin_size Size of the genomic bins to use (in bp). Each chromsome will be chunked into bins with size which is
#' smaller than this value. Default is 50Mb.
#' @param overwrite overwrite existing directory (optional)
#' @param verbose verbose output (optional)
#' @param tabix_bin path to the tabix binary (optional)
#' @param chromosomes process only specific chromosomes (optional).
#'
#' @return None
#'
#' @examples
#' \dontrun{
#' write_sc_counts_from_fragments("pbmc_data/fragments.tsv.gz", "pbmc_reads", cell_names = atac_sc)
#' write_sc_counts_from_fragments("pbmc_data/fragments.tsv.gz", "pbmc_reads", cell_names = atac_sc, chromosomes = paste0("chr", c(1:22, "X", "Y")))
#' }
#'
#' @export
write_sc_counts_from_fragments <- function(fragments_file, out_dir, cell_names, cell_suffix = NULL, genome = NULL, bin_size = 5e7, overwrite = FALSE, id = "", description = "", verbose = FALSE, tabix_bin = "tabix", chromosomes = gintervals.all()) {
    withr::local_options(list(scipen = 1e5))
    data_dir <- file.path(out_dir, "data")
    if (dir.exists(out_dir)) {
        if (!overwrite) {
            cli_abort("Output directory {.file {out_dir}} already exists. Use 'overwrite = TRUE' to overwrite.")
        } else {
            cli_warn("Output directory {.file {out_dir}} already exists. Overwriting.")
            unlink(out_dir, recursive = TRUE)
        }
    }

    dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
    cli_alert("Writing to {.file {out_dir}}")

    if (class(cell_names) == "ScPeaks") {
        genome <- genome %||% cell_names@genome
        cell_names <- colnames(cell_names@mat)
    }
    gset_genome(genome)
    use_tabix <- bin_exists(tabix_bin) && file.exists(paste0(fragments_file, ".tbi")) && grepl("gz$", tools::file_ext(fragments_file))

    chromosomes <- chromosomes %||% gintervals.all()$chrom

    genomic_bins <- giterator.intervals(iterator = bin_size, intervals = gintervals.all() %>% filter(chrom %in% chromosomes)) %>%
        mutate(name = glue("{chrom}_{start}_{end}"))

    cli_alert_info("Processing {.val {nrow(genomic_bins)}} genomic bins of maximal size {.val {bin_size}}")

    plyr::a_ply(genomic_bins, 1, function(region) {
        cli_alert("Processing {.val {region$name}}")
        write_sparse_matrix_from_fragments(
            fragments_file,
            glue("{data_dir}/{region$chrom}_{region$start}_{region$end}.mtx"),
            cell_names,
            region = region,
            genome = genome,
            overwrite = overwrite,
            verbose = verbose,
            tabix_bin = tabix_bin,
            use_tabix = use_tabix
        )
        cli_alert_success("Finished processing {.val {region$name}}")
    }, .parallel = getOption("mcatac.parallel"))

    if (!is.null(cell_suffix)) {
        cell_names <- paste0(cell_names, "-", cell_suffix)
    }

    counts_md <- list(
        cell_names = cell_names,
        genome = genome,
        id = id,
        description = description,
        data_dir = "data",
        genomic_bins = genomic_bins$name
    )

    yaml::write_yaml(counts_md, file = file.path(out_dir, "metadata.yaml"))

    cli_alert_info("Created sparse matrices for {.val {nrow(genomic_bins)}} genomic bins")
    cli_alert_success("Created an ScCounts object at {.file {out_dir}}")
}

#' Write single cell counts from an indexed bam file
#'
#' @param out_dir output directory.
#' @param id an identifier for the object, e.g. "pbmc".
#' @param description description of the object, e.g. "PBMC from a healthy donor - granulocytes removed through cell sorting (10k)"
#' @param bin_size Size of the genomic bins to use (in bp). Each chromsome will be chunked into bins with size which is
#' smaller than this value. Default is 50Mb.
#' @param overwrite overwrite existing directory (optional)
#'
#' @return None
#'
#' @inheritParams write_sparse_matrix_from_bam
#' @inheritParams write_sc_counts_from_fragments
#'
#' @examples
#' \dontrun{
#' write_sc_counts_from_bam("pbmc_data/possorted_bam.bam", "pbmc_reads", cell_names = atac_sc)
#' }
#'
#' @export
write_sc_counts_from_bam <- function(bam_file, out_dir, cell_names, cell_suffix = NULL, genome = NULL, bin_size = 5e7, id = "", description = "", min_mapq = NULL, samtools_bin = "samtools", samtools_opts = NULL, num_reads = NULL, verbose = FALSE, overwrite = FALSE) {
    withr::local_options(list(scipen = 1e5))
    data_dir <- file.path(out_dir, "data")
    if (dir.exists(out_dir)) {
        if (!overwrite) {
            cli_abort("Output directory {.file {out_dir}} already exists. Use 'overwrite = TRUE' to overwrite.")
        } else {
            cli_warn("Output directory {.file {out_dir}} already exists. Overwriting.")
            unlink(out_dir, recursive = TRUE)
        }
    }

    dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
    cli_alert("Writing to {.file {out_dir}}")

    if (class(cell_names) == "ScPeaks") {
        genome <- genome %||% cell_names@genome
        cell_names <- colnames(cell_names@mat)
    }
    gset_genome(genome)

    chroms <- intersect(gintervals.all()$chrom, bam_file_chromosomes(bam_file))


    genomic_bins <- giterator.intervals(iterator = bin_size, intervals = gintervals.all() %>% filter(chrom %in% chroms)) %>%
        mutate(name = glue("{chrom}_{start}_{end}"))

    cli_alert_info("Processing {.val {nrow(genomic_bins)}} genomic bins of maximal size {.val {bin_size}}")

    plyr::a_ply(genomic_bins, 1, function(region) {
        cli_alert("Processing {.val {region$name}}")
        write_sparse_matrix_from_bam(bam_file, glue("{data_dir}/{region$chrom}_{region$start}_{region$end}.mtx"), cell_names, region, genome, min_mapq, samtools_bin, samtools_opts, num_reads, verbose, overwrite = overwrite)
        cli_alert_success("Finished processing {.val {region$name}}")
    }, .parallel = getOption("mcatac.parallel"))

    if (!is.null(cell_suffix)) {
        cell_names <- paste0(cell_names, "-", cell_suffix)
    }

    counts_md <- list(
        cell_names = cell_names,
        genome = genome,
        id = id,
        description = description,
        data_dir = "data",
        genomic_bins = genomic_bins$name
    )

    yaml::write_yaml(counts_md, file = file.path(out_dir, "metadata.yaml"))

    cli_alert_info("Created sparse matrices for {.val {nrow(genomic_bins)}} genomic bins")
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
bam_file_chromosomes <- function(bam_file, samtools_bin = "samtools") {
    cmd <- glue("{samtools_bin} view -H {bam_file} | grep '^@SQ' | cut -f 2 | cut -d ':' -f 2 | cut -d ' ' -f 1 | sort | uniq", samtools_bin = samtools_bin)
    return(system(cmd, intern = TRUE))
}
