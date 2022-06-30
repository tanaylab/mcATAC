#' Subset 10X atac_possorted.bam file into per-metacell BAM files
#'
#' @param bam_path path to the 10X atac_possorted.bam file
#' @param out_dir (optional) directory to output per-metacell bam files. Default: /bam_path/mcatac@id_mc_bams
#' @param gparallel_path (optional) path to GNU parallel
#' @param samtools_path (optional) path to samtools (requires a version which supports the -D/--tag-file option in samtools view)
#' @return error log
#' @inheritParams write_metacell_cell_names
#' @export
generate_per_metacell_bams <- function(bam_path,
                                       mcatac,
                                       out_dir = file.path(bam_path, paste0(mcatac@id, "_mc_bams")),
                                       gparallel_path = "/usr/wisdom/parallel",
                                       samtools_path = "/home/feshap/src/samtools-1.15.1/samtools") {
    if (!file.exists(paste0(bam_path, ".bai"))) {
        cli_abort("Index file not found for {.file {bam_path}}. Please run 'samtools index {bam_path}'.")
    }
    c2mc_path <- file.path(dirname(bam_path), paste0(mcatac@id, "_c2mc"))
    write_metacell_cell_names(mcatac, c2mc_path)
    fp_sb2mc <- system.file("exec", "split_bam_to_metacells.sh", package = "mcATAC")
    fp_rgp <- system.file("exec", "run_gnu_parallel.sh", package = "mcATAC")
    error_log <- withr::with_path(
        new = gparallel_path,
        code = system2(
            command = fp_sb2mc,
            args = c(
                bam_path,
                c2mc_path,
                out_dir,
                fp_rgp,
                gparallel_path,
                samtools_path
            )
        )
    )
    system(glue::glue("rm -rf {c2mc_path}"))
    return(error_log)
}

#' Write WIGs from BAMs
#'

#' @param bam_folder_path path to BAMs
#' @param output_path (optional) where to put WIGs. Default - bam_folder_path/wig_output
#' @param parallel (optional) whether to do parallel computations
#' @param nc (optional) number of cores for parallel computations
#' @return error log
#' @export
generate_wigs_from_bams <- function(bam_folder_path,
                                    output_path = file.path(bam_folder_path, "wig_output"),
                                    parallel = getOption("mcatac.parallel"),
                                    nc = getOption("mcatac.parallel.nc")) {
    bam_folder_path <- normalizePath(bam_folder_path)
    output_path <- normalizePath(output_path)
    if (!dir.exists(output_path)) {
        dir.create(output_path, recursive = TRUE)
    }
    bams <- list.files(bam_folder_path, pattern = "\\.bam$")
    if (parallel) {
        error_log <- parallel::mclapply(bams, FUN = function(fl) {
            output_filename <- file.path(output_path, gsub("\\.bam$", ".wig", basename(fl)))
            bam_to_wig(file.path(bam_folder_path, fl), output_filename)
        }, mc.cores = nc)
    } else {
        error_log <- sapply(bams, function(fl) bam_to_wig(file.path(bam_folder_path, fl), output_filename))
    }
    return(error_log)
}


#' Utility function to write metacell names to files
#'
#' @param mcatac an McPeaks object (from which to get metacell assignments)
#' @param c2mc_path (optional) path to write out metacell assignments to
#'
write_metacell_cell_names <- function(mcatac, c2mc_path = NULL) {
    if (is.null(c2mc_path)) {
        c2mc_path <- file.path("data", paste0(mcatac@id, "_c2mc"))
    }
    if (!dir.exists(c2mc_path)) {
        dir.create(c2mc_path, recursive = TRUE)
    }
    c2mc <- mcatac@cell_to_metacell
    dummy <- sapply(sort(unique(c2mc$metacell)), function(mci) {
        nmc <- c2mc$cell_id[c2mc$metacell == mci]
        write(nmc, file = file.path(c2mc_path, paste0("mc", mci, ".txt")), sep = "\n")
        c2mc_path <- file.path("data", paste0(mcatac@id, "_c2mc"))
    })
    if (!dir.exists(c2mc_path)) {
        dir.create(c2mc_path, recursive = TRUE)
    }
    c2mc <- mcatac@cell_to_metacell
    dummy <- sapply(sort(unique(c2mc$metacell)), function(mci) {
        nmc <- c2mc$cell_id[c2mc$metacell == mci]
        write(nmc, file = file.path(c2mc_path, paste0("mc", mci, ".txt")), sep = "\n")
    })
    return(NULL)
}

#' Merge BAMs by metacells
#'
#' @param bam_path path to the per-metacell BAM tracks
#' @param output_filename what to call the output BAM file
#' @param mcs vector of mcs to merge BAMs of
#' @inheritParams generate_wigs_from_bams
#' @export

merge_metacell_bams <- function(bam_path, output_filename, mcs, parallel = getOption("mcatac.parallel"), nc = getOption("mcatac.parallel.nc")) {
    file_list <- paste(paste0(bam_path, "/mc", mcs, ".bam"))
    if (!grepl(".bam$", output_filename)) {
        output_filename <- paste0(output_filename, ".bam")
    }
    if (all(file.exists(file_list))) {
        if (!parallel) {
            nc <- 0
        }
        command <- glue::glue("samtools merge -f --threads {nc}")
        command <- paste(command, output_filename, paste0(file_list, collapse = " "))
        system(command)
    } else {
        bad_mcs <- mcs[!file.exists(file_list)]
        cli_abort("BAM files for MCs {.val {bad_mcs}} were not found in {.val {bam_path}}")
    }
    return(NULL)
}

#' Convert BAM file to WIG file
#'
#' @param bam_path path to the BAM file
#' @param output_filename (optional) what to call the output WIG file. By default it saves the file at the same directory and with the same name as the .bam input (with the appropriate file type suffix).
#' @param track_name_prefix (optional) name to prepend to track (e.g. for display in UCSC genome browser)
#' @export
bam_to_wig <- function(bam_path, output_filename = NULL, track_name_prefix = NULL, header_line = NULL) {
    # base_name <- gsub('.bam$', '', fl)
    base_name <- basename(bam_path)
    dir_name <- dirname(normalizePath(bam_path))
    if (is.null(output_filename)) {
        output_filename <- file.path(dir_name, gsub("\\.bam", ".wig", base_name))
        print(glue::glue("Converting {base_name} to wig. Saving at {output_filename}."))
    }
    # full_path = file.path(output_path, paste0(base_name, '.wig'))
    if (is.null(header_line)) {
        if (is.null(track_name_prefix)) {
            header_line <- "track type=wiggle_0 name='{base_name}'"
        } else {
            header_line <- "track type=wiggle_0 name='{track_name_prefix} - {base_name}'"
        }
    }
    hli <- glue::glue(header_line)
    create_file_command <- glue::glue('echo "{hli}" > {output_filename}')
    system(command = create_file_command)
    smt_cmd <- glue::glue("samtools mpileup -BQ0 {bam_path}")
    perl_cmd <- paste0(
        "perl -pe '($c, $start, undef, $depth) = split;if ($c ne $lastC || $start != $lastStart+1) {print ",
        '"fixedStep chrom=$c start=$start step=1 span=1\n"',
        ";}$_ = $depth.\"\n\";($lastC, $lastStart) = ($c, $start);'"
    )
    out_cmd <- glue::glue(" >> {output_filename}")
    bam2wig_command <- paste0(paste0(c(smt_cmd, perl_cmd), collapse = " | "), out_cmd)
    system(command = bam2wig_command)
    return(NULL)
}

#' Convert WIG files to misha tracks
#'
#' @param wig_dir path to folder containing WIG files
#' @param track_name_prefix name to prepend to track.
#' @param description (optional) description for misha tracks
#' @param parallel (optional) whether to use parallel computation
#' @param force (optional) whether to force rewrite of existing misha tracks
#' @return error log
#' @inheritParams generate_wigs_from_bams
#' @export
convert_wigs_to_tracks <- function(wig_dir,
                                   track_name_prefix = NULL,
                                   description = NULL,
                                   parallel = TRUE,
                                   force = FALSE,
                                   nc = parallel::detectCores()) {
    wig_file_paths <- file.path(wig_dir, grep("\\.wig$", list.files(wig_dir), v = T))
    track_name_suffix <- "unnorm"
    if (is.null(track_name_prefix)) {
        track_name_prefix <- basename(wig_dir)
    }
    if (is.null(description)) {
        description <- "ATAC fragment track for {track_name_prefix} - {bn}"
    }
    tracks_folder <- file.path(GWD, track_name_prefix)
    if (!dir.exists(tracks_folder)) {
        misha.ext::gtrack.create_dirs(tracks_folder)
    }
    if (parallel) {
        error_log <- parallel::mclapply(wig_file_paths,
            FUN = function(x) {
                make_misha_track_from_wig(x,
                    track_name_prefix = track_name_prefix,
                    track_name_suffix = track_name_suffix,
                    description = description,
                    force = force
                )
            },
            mc.cores = nc
        )
    } else {
        error_log <- sapply(wig_file_paths, make_misha_track_from_wig)
    }
    gdb.reload()
    temp_tracks <- gtrack.ls(glue::glue("{track_name_prefix}.*{track_name_suffix}"))
    error_log_2 <- normalize_tracks(tracks = temp_tracks, track_name_suffix = track_name_suffix)
    dummy <- sapply(temp_tracks, gtrack.rm, f = T)
    gdb.reload()
    return(list(error_log, error_log_2))
}


# Convert one WIG file to misha track
#' @param fp path to WIG file
#' @inheritParams convert_wigs_to_tracks
#' @noRd
make_misha_track_from_wig <- function(fp, track_name_prefix = NULL, description = NULL, force = FALSE) {
    bn <- gsub(".wig$", "", basename(fp))
    track_name <- glue::glue("{track_name_prefix}.{bn}")
    track_exists <- gtrack.exists(track_name)
    if (track_exists) {
        if (force) {
            cli_alert_warning("Deleting track {.val track_name}")
            gtrack.rm(track_name, force = force)
            gdb.reload()
        } else {
            cli_alert_info("Track {.val track_name} exists and {.var force} is {.val force}. Skipping import of track.")
        }
    }
    dummy <- gtrack.import(
        track = track_name,
        description = glue::glue(description),
        file = fp,
        binsize = 0
    )
    return(dummy)
}

# Normalize ATAC track to sum of fragments
#' @param tracks path to WIG file
#' @inheritParams make_misha_track_from_wig
#' @inheritParams convert_wigs_to_tracks
#' @noRd
normalize_tracks <- function(tracks, track_name_suffix, parallel = TRUE, nc = parallel::detectCores()) {
    new_track_names <- gsub(glue::glue("_{track_name_suffix}"), "", tracks)
    if (parallel) {
        gsums <- parallel::mclapply(tracks, gsummary, mc.cores = nc)
        error_log <- parallel::mcmapply(
            FUN = function(ntrki, trki, sumi) {
                gtrack.create(
                    track = ntrki,
                    description = glue::glue("Track {trki} normalized to library size"),
                    expr = glue::glue("{trki}/{sumi[['Sum']]}")
                )
            },
            new_track_names,
            tracks,
            gsums,
            mc.cores = nc
        )
    } else {
        gsums <- lapply(tracks, gsummary)
        error_log <- mapply(
            FUN = function(ntrki, trki, sumi) {
                gtrack.create(
                    track = ntrki,
                    description = glue::glue("Track {trki} normalized to library size"),
                    expr = glue::glue("{trki}/{sumi[['Sum']]}")
                )
            },
            new_track_names,
            tracks,
            gsums
        )
    }
    gdb.reload()
    return(error_log)
}
