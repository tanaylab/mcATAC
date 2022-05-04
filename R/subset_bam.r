#' Subset 10X atac_possorted.bam file into per-metacell BAM files
#' 
#' @param bam_path path to the 10X atac_possorted.bam file
#' @param mcatac McATAC object
#' @param out_dir (optional) directory to output per-metacell bam files
#' @param c2mc_path (optional) directory to output cell-to-metacell mappings
#' @inheritParams write_metacell_cell_names
generate_per_metacell_bams <- function(bam_path, mcatac, out_dir = NULL, c2mc_path = NULL) {
    if (is.null(out_dir)) {
        out_dir <- file.path(bam_path, paste0(mcatac@id, "_mc_bams"))
    }
    if (is.null(c2mc_path)) {
        c2mc_path <- file.path(dirname(bam_path), paste0(mcatac@id, "_c2mc"))
    }
    write_metacell_cell_names(mcatac, c2mc_path)
    fp_sb2mc <- system.file("exec", "split_bam_to_metacells.sh", package = "mcATAC")
    fp_rgp <- system.file("exec", "run_gnu_parallel.sh", package = "mcATAC")
    error_log <- withr::with_path(new = "/usr/wisdom/parallel", 
                                    code = system2(command = fp_sb2mc,
                                                    args = c(bam_path, 
                                                             c2mc_path, 
                                                             out_dir,
                                                             fp_rgp))
                                )
    return(error_log)
}

#' Write WIGs from BAMs
#' 
#' @param bam_folder_path path to BAMs
#' @param output_path (optional) where to put WIGs
#' @param parallel (optional) whether to do parallel computations
#' @inheritParams bam_to_wig
generate_wigs_from_bams <- function(bam_folder_path, track_name_prefix, output_path = NULL, parallel = TRUE) {
    bam_folder_path = normalizePath(bam_folder_path)
    if (is.null(output_path)) {
        output_path = file.path(bam_folder_path, "wig_output")
    }
    output_path = normalizePath(output_path)
    if (!dir.exists(output_path)) {dir.create(output_path)}
    bams <- grep('\\.bam$', list.files(bam_folder_path), v=T)
    if (parallel) {
        p = parallel::detectCores()
        error_log <- parallel::mclapply(bams, FUN = function(fl) {
                        output_filename = file.path(output_path, gsub('\\.bam$', '.wig', basename(fl)));
                        bam_to_wig(file.path(bam_folder_path, fl), output_filename)
        }, mc.cores = p)
    }
    else {
        error_log <- sapply(bams, function(fl) bam_to_wig(file.path(bam_folder_path, fl), output_filename))
    }
    return(error_log)
}


#' Utility function to write metacell names to files
#' 
#' @param mcatac an McATAC object (from which to get metacell assignments)
#' @param c2mc_path (optional) path to write out metacell assignments to
write_metacell_cell_names <- function(mcatac, c2mc_path = NULL) {
    if (is.null(c2mc_path)) {
        c2mc_path <- file.path('data', paste0(mcatac@id, "_c2mc"))
    }
    if (!dir.exists(c2mc_path)) {dir.create(c2mc_path)}
    c2mc <- mcatac@cell_to_metacell
    dummy <- sapply(sort(unique(c2mc$metacell)), function(mci) {
        nmc <- c2mc$cell_id[c2mc$metacell == mci]
        write(nmc,file=file.path(c2mc_path, paste0("mc", mci, ".txt")),sep='\n')
    })
}

#' Merge BAMs by metacells
#' 
#' @param bam_path path to the per-metacell BAM tracks
#' @param output_filename what to call the output BAM file
#' @param mcs vector of mcs to merge BAMs of
#' @param parallel (optional) whether to do parallel computations
merge_metacell_bams <- function(bam_path, output_filename, mcs, parallel = TRUE) {
    file_list <- paste(paste0(bam_path,'/mc',mcs, ".bam"))
    if (!grepl(".bam$", output_filename)) {
        output_filename <- paste0(output_filename, ".bam")
    }
    if (all(file.exists(file_list))) {
        if (parallel) {
            p = parallel::detectCores() - 1
        }
        else {p = 0}
        command <- glue::glue("samtools merge -f --threads {p}")
        command <-paste(command, output_filename, paste0(file_list, collapse = ' '))
        system(command)
    }
    else {
        bad_mcs <- mcs[!file.exists(file_list)]
        cli_abort("BAM files for MCs {.val {bad_mcs}} were not found in {.val {bam_path}}")
    }
}

#' Convert BAM file to WIG file
#' 
#' @param bam_path path to the BAM file
#' @param output_filename what to call the output WIG file
#' @param track_name_prefix (optional) name to prepend to track (e.g. for display in UCSC genome browser)
bam_to_wig = function(bam_path, output_filename, track_name_prefix = NULL, header_line = NULL) {
    # base_name <- gsub('.bam$', '', fl)
    base_name = basename(bam_path)
    dir_name = dirname(normalizePath(bam_path))
    if (is.null(output_filename)) {
        output_filename = file.path(dir_name, gsub('\\.bam', '.wig', base_name))
        print(glue::glue("Converting {base_name} to wig. Saving at {output_filename}."))
    }
    # full_path = file.path(output_path, paste0(base_name, '.wig'))
    if (is.null(header_line)) {
        if (is.null(track_name_prefix)) {
            header_line <- "track type=wiggle_0 name='{base_name}'"
        }
        else {
            header_line <- "track type=wiggle_0 name='{track_name_prefix} - {base_name}'"
        }
    }
    hli = glue::glue(header_line)
    create_file_command <- glue::glue('echo "{hli}" > {output_filename}')
    system(command = create_file_command)
    smt_cmd = glue::glue("samtools mpileup -BQ0 {bam_path}")
    perl_cmd = paste0("perl -pe '($c, $start, undef, $depth) = split;if ($c ne $lastC || $start != $lastStart+1) {print ", 
                                '"fixedStep chrom=$c start=$start step=1 span=1\n"',
                                ";}$_ = $depth.\"\n\";($lastC, $lastStart) = ($c, $start);'")
    out_cmd = glue::glue(" >> {output_filename}")
    bam2wig_command = paste0(paste0(c(smt_cmd, perl_cmd), collapse = " | "), out_cmd)
    system(command = bam2wig_command)
}