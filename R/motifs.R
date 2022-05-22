#' Generate peak-motif matrix
#'
#' This function calculates aggregated log-PSSM energies (as derived by the misha package's \code{create_pssm_energy} module) for a set of peaks
#' (default - all peaks in the dataset) and for a set of motifs (default - all available)
#'
#' @param atac - an ScATAC/McATAC or PeakIntervals object
#' @param peak_width (optional) - size of region around peak centers to extract motif energies for
#' @param pssm_path (optional) - path to directory containing misha-formatted pssm files (e.g. \code{motifs.key} and \code{motifs.data})
#' @param datasets_of_interest (optional) - names of pssm datasets ([name].key-[name].data file combinations) to calculate PSSM values for
#' @param motif_tracks (optional) - misha track names for which to extract motif PSSMs
#' @param motif_regex (optional) - a vector of regular expressions for which to match motif track names and extract motif PSSMs
#' @param parallel (optional) - whether to use parallel computations
#' @param nc (optional) - number of cores to use for parallel computations
#' @return a matrix of peaks (rows) vs. aggregated motif energies (columns)
#' @examples
#' \dontrun{
#' peak_motif_mat <- generate_motif_pssm_matrix(
#'     peak_set = head(atac_sc@peaks),
#'     motif_regex = c("Bcl", "Atf"),
#'     datasets_of_interest = c("homer", "jaspar", "jolma"),
#'     parallel = F
#' )
#' }
#' @export
generate_motif_pssm_matrix <- function(atac = NULL,
                                       peak_width = 200,
                                       pssm_path = NULL,
                                       datasets_of_interest = NULL,
                                       motif_tracks = NULL,
                                       motif_regex = NULL,
                                       parallel = getOption("mcatac.parallel"),
                                       nc = getOption("mcatac.parallel.nc")) {
    cli_alert_warning("The runtime of this function call may take several minutes, depending on the numbers of peaks and motifs evaluated, and the number of processors available. Consider running it in a separate terminal and saving the output.")
    withr::local_options(list(gmax.data.size = 1e8))
    suffix <- stringi::stri_rand_strings(1, 5)
    peaks <- get_peaks_for_pssm(atac)
    if (!is.null(motif_tracks)) {
        tracks_exist <- sapply(motif_tracks, gtrack.exists)
        bad_tracks <- motif_tracks[!tracks_exist]
        tracks <- motif_tracks[tracks_exist]
        if (length(bad_tracks) > 0) {
            cli_alert_warning("Tracks {.val {paste(head(bad_tracks))}} and {.val {pmax(0, length(bad_tracks) - 6)}} more were not found in misha gdb. Check that you are querying appropriate motif tracks.")
        }
        if (length(tracks) == 0) {
            cli_abort("No tracks matching {.var motif_tracks} were found")
        } else {
            peak_motif_matrix <- gextract(tracks, intervals = peaks, iterator = peaks, colnames = tracks)
            return(peak_motif_matrix)
        }
    } else if (!is.null(motif_regex)) {
        all_pssms <- get_available_pssms(pssm_path, datasets_of_interest)
        pssm_filt_inds <- lapply(all_pssms[["keys"]], function(kdf) {
            sapply(motif_regex, function(mre) {
                kdf$key[grep(mre, kdf$track, ignore.case = TRUE)]
            }) %>% unlist()
        })
        keys_filt <- lapply(seq_along(all_pssms[["keys"]]), function(i) {
            filter(all_pssms[["keys"]][[i]], key %in% pssm_filt_inds[[i]])
        })
        datasets_filt <- lapply(seq_along(all_pssms[["datasets"]]), function(i) {
            filter(all_pssms[["datasets"]][[i]], key %in% pssm_filt_inds[[i]])
        })
        names(keys_filt) <- names(all_pssms[["keys"]])
        names(datasets_filt) <- names(all_pssms[["datasets"]])
        pssms <- list("keys" = keys_filt, "datasets" = datasets_filt)
    } else {
        cli_alert_info("No motif tracks or regular expressions were specified. Extracting PSSMs for all available motifs.")
        pssms <- get_available_pssms(pssm_path, datasets_of_interest)
    }
    peak_mids <- round((peaks$end + peaks$start) / 2)
    peaks <- mutate(peaks, start = round(peak_mids - peak_width / 2), end = round(peak_mids + peak_width / 2))
    peaks <- fix_missing_chroms_in_peaks(peaks)
    if (!parallel) {
        nc <- 2
    }
    res <- mapply(FUN = function(kdf, n) {
        parallel::mcmapply(FUN = function(.x, .y) {
            if (!dir.exists(file.path(GWD, n))) {
                dir.create(file.path(GWD, n))
            }
            track_name <- paste0(n, ".", gsub("-|\\||\\.", "_", .y), "_", suffix)
            if (!gtrack.exists(track_name)) {
                gtrack.create_pwm_energy(
                    track = track_name,
                    description = glue::glue("Temporary track created for {track_name}"),
                    pssmset = n,
                    pssmid = .x,
                    prior = 0.01,
                    iterator = peaks
                )
            }
            return(track_name)
        }, kdf$key, kdf$track, mc.cores = max(1, nc - 1))
    }, pssms[["keys"]], names(pssms[["keys"]])) %>% unlist()
    gdb.reload()
    trks_motifs <- res
    if (length(trks_motifs) >= 1e+2) {
        divs <- 1:4
        trk_inds <- unlist(do.call("c", lapply(divs, function(i) {
            ifelse(as.numeric(i) != max(divs),
                return(rep.int(x = i, times = ceiling(length(trks_motifs) / length(divs)))),
                return(rep.int(x = i, times = floor(length(trks_motifs) / length(divs))))
            )
        })))
        df_list <- lapply(sort(unique(trk_inds)), function(u) {
            gextract(trks_motifs[trk_inds == u],
                intervals = peaks,
                iterator = peaks,
                colnames = gsub(paste0("_", suffix), "", trks_motifs[trk_inds == u])
            )
        })
        peak_motif_matrix <- plyr::join_all(df_list, by = c("chrom", "start", "end"), type = "left")
    } else {
        peak_motif_matrix <- gextract(res, intervals = peaks, iterator = peaks, colnames = gsub(paste0("_", suffix), "", trks_motifs))
    }
    withr::defer(purrr::walk(trks_motifs, gtrack.rm, force = TRUE))
    peak_motif_matrix <- peak_motif_matrix[peak_motif_matrix$end - peak_motif_matrix$start >= peak_width, ]
    return(peak_motif_matrix)
}

#' Generate random genome motif PSSM matrix
#'
#' @param num_peaks total number of intervals (will be divided proportionately between chromosomes) to use as background
#' @param bp_from_chrom_edge_to_avoid regions (in bp) from edges of chromosomes to avoid sampling from (e.g. avoid acrocentric centromeres)
#' @inheritParams generate_motif_pssm_matrix
#' @return a matrix of peaks (rows) vs. aggregated motif energies (columns)
#' @examples
#' \dontrun{
#' random_genome_motifs <- gen_random_genome_peak_motif_matrix(num_peaks = 5e+4, datasets_of_interest = get_available_pssms(return_datasets_only = TRUE))
#' }
#' @export
gen_random_genome_peak_motif_matrix <- function(num_peaks = 1e+5,
                                                peak_width = 2e+2,
                                                bp_from_chrom_edge_to_avoid = 3e+6,
                                                datasets_of_interest = NULL,
                                                motif_regex = NULL,
                                                motif_tracks = NULL,
                                                parallel = TRUE) {
    ALLGENOME[[1]] <- ALLGENOME[[1]][!grepl("_", ALLGENOME[[1]]$chrom), ]
    chrom_lens <- apply(ALLGENOME[[1]][, 2:3], 1, diff)
    chrom_fracs <- setNames(chrom_lens / sum(chrom_lens), ALLGENOME[[1]][, 1])
    sample_seqs <- mapply(chrom_fracs, names(chrom_fracs), chrom_lens, FUN = function(x, y, z) {
        chrom <- rep(y, round(x * num_peaks))
        start <- sample.int(n = z, size = length(chrom))
        end <- start + peak_width
        return(as.data.frame(rbind(chrom, start, end)))
    })
    sample_seqs <- as.data.frame(do.call("rbind", lapply(sample_seqs, t)))
    sample_seqs[, 2:3] <- apply(sample_seqs[, 2:3], 2, as.numeric)
    sample_seqs <- sample_seqs[with(sample_seqs, order(chrom, start, end)), ]
    end_shift <- ALLGENOME[[1]][match(sample_seqs$chrom, ALLGENOME[[1]][, 1]), 3] - bp_from_chrom_edge_to_avoid
    sample_seqs <- sample_seqs[sample_seqs$start >= bp_from_chrom_edge_to_avoid & sample_seqs$end <= end_shift, ]
    sample_seqs <- PeakIntervals(sample_seqs)
    rg_peak_motif_mat <- generate_motif_pssm_matrix(
        atac = sample_seqs,
        peak_width = peak_width,
        motif_tracks = motif_tracks,
        datasets_of_interest = datasets_of_interest,
        motif_regex = motif_regex,
        parallel = parallel
    )
    return(rg_peak_motif_mat)
}


#' Calculate Kolmogorov-Smirnov D statistics between two interval sets with motif energies
#'
#' This function does a one-sided KS test between a foreground set of peaks (\code{pssm_fg}) and a background set \code{pssm_bg}.
#' The option \code{alternative == "less"}, checks the null hypothesis that the foreground distribution is not less than the
#' background distribution (applicable when looking for motif enrichment; for anti-enrichment, \code{alternative == 'greater'},
#' see ks.test documentation for further details)
#' @param pssm_fg motif energies calculated for a certain set of motifs on a PeakIntervals/ScATAC/McATAC object
#' @param pssm_bg a background set of intervals (e.g. random genome, all ENCODE enhancers etc.) that include all/subset of the motifs (columns) in pssm_fg
#' @param fg_clustering a vector of cluster assignments for the foreground peaks (e.g. from \code{gen_atac_peak_clust})
#' @param parallel (optional) - whether to use parallelize computations
#' @param nc (optional) - number of cores for parallel computations
#' @inheritParams stats::ks.test
#' @return if \code{fg_clustering == TRUE}, returns a matrix of clusters x motifs (rows x columns) with the D-statistic for each combination
#' @examples
#' \dontrun{
#' pssm_fg <- generate_motif_pssm_matrix(my_atac_mc, datasets_of_interest = "jaspar")
#' pssm_bg <- gen_random_genome_peak_motif_matrix(num_peaks = nrow(my_atac_mc@peaks), datasets_of_interest = "jaspar")
#' d_vs_rg <- calculate_d_stats(pssm_fg, pssm_bg)
#' peak_clust <- gen_atac_peak_clust(my_atac_mc, k = 12)
#' d_vs_rg_cl <- calculate_d_stats(pssm_fg, pssm_bg, fg_clustering = peak_clust)
#' }
#' @export
calculate_d_stats <- function(pssm_fg, pssm_bg, fg_clustering = NULL, parallel = getOption("mcatac.parallel"), alternative = "less", nc = getOption("mcatac.parallel.nc")) {
    cols_fg <- grep("chrom|start|end$|interval", colnames(pssm_fg), ignore.case = T, invert = T, value = T)
    cols_bg <- grep("chrom|start|end$|interval", colnames(pssm_bg), ignore.case = T, invert = T, value = T)
    cols_both <- intersect(cols_fg, cols_bg)
    if (!parallel) {
        nc <- pmax(2, round(0.1 * nc))
    }
    if (!is.null(fg_clustering)) {
        ks_test_results <- parallel::mclapply(cols_both, FUN = function(x, i) {
            return(tapply(x[, i], fg_clustering, function(y) suppressWarnings(ks.test(y, pssm_bg[, i], alternative = alternative))))
        }, mc.cores = nc, x = pssm_fg)
        ks_d <- sapply(ks_test_results, function(x) sapply(x, function(y) y$statistic))
        colnames(ks_d) <- cols_both
    } else {
        ks_test_results <- parallel::mclapply(cols_both, function(x) {
            ks.test(x = pssm_fg[, x], y = pssm_bg[, x], alternative = alternative)
        },
        mc.cores = nc
        )
        ks_d <- setNames(sapply(ks_test_results, function(x) x$statistic), cols_both)
    }
    return(ks_d)
}

#' Utility function for validating parameter combinations for PSSM extraction
#'
#' @param atac - an ScATAC/McATAC or PeakIntervals object
#' @return the relevant peak set (from ATAC object if specified, else \code{peak_set})
#' @noRd
get_peaks_for_pssm <- function(atac) {
    cl <- class(atac)
    if (any(c("ScATAC", "McATAC") %in% cl[[1]])) {
        peaks <- atac@peaks
    } else if ("PeakIntervals" %in% cl[[1]]) {
        peaks <- atac
    } else {
        cli_abort("Class of {.var atac} is not recognized (should be either ScATAC, McATAC or PeakIntervals object).")
    }
    return(peaks)
}

#' Function for validating \code{pssm_path} and getting all/selected PSSM datasets and tracks from path/gdb
#'
#' @param pssm_path (optional) - path to misha pssm datasets (*.key-*.data file combinations), no need to specify if misha gdb is set up correctly
#' @param datasets_of_interest (optional) - pssm datasets to get available pssms from. Available datasets can be obtained by calling this function with \code{return_datasets_only = TRUE}
#' @param return_datasets_only (optional) - whether to return only available pssm datasets
#' @return either: 1) (default) a list containing a list of the names and a list of corresponding dataframes of the motif datasets available, OR 2) if \code{return_datasets_only == TRUE}, just the names of the available motif datasets
#' @examples
#' \dontrun{
#' pssm_datasets <- get_available_pssms(return_datasets_only = TRUE)
#' all_pssms <- get_available_pssms()
#' all_pssms <- get_available_pssms(datasets_of_interest = pssm_datasets)
#' all_pssms_jaspar_homer <- get_available_pssms(datasets_of_interest = c("jaspar", "homer"))
#' }
#' @export
get_available_pssms <- function(pssm_path = NULL, datasets_of_interest = NULL, return_datasets_only = FALSE) {
    if (is.null(pssm_path)) {
        pssm_path <- file.path(GROOT, "pssms")
        if (!dir.exists(pssm_path)) {
            cli_abort("The default {.var pssm_path} directory ({.val pssm_path}) does not exist. Make sure you have set up misha gdb correctly.")
        }
        if (is.null(datasets_of_interest) && !return_datasets_only) {
            cli_abort("Must specify either {.var pssm_path} or {.var datasets_of_interest} or set {.var return_datasets_only = TRUE}.")
        }
    }
    if (!dir.exists(pssm_path)) {
        cli_abort("The specified {.var pssm_path} directory ({.val pssm_path}) does not exist.")
    }
    pssm_file_list <- list.files(pssm_path)
    motif_keys <- grep(".key", pssm_file_list, value = TRUE)
    motif_datasets <- grep(".data", pssm_file_list, value = TRUE)
    unique_keys <- sort(unique(gsub(".key", "", motif_keys)))
    unique_datasets <- sort(unique(gsub(".data", "", motif_datasets)))
    good_keys <- unique_keys[!is.na(match(unique_keys, unique_datasets))]
    if (length(good_keys) == 0) {
        cli_abort("There are no matching *.key-*.data file name combinations")
    }
    if (return_datasets_only) {
        return(good_keys)
    }
    if (!is.null(datasets_of_interest)) {
        pssm_datasets_in <- good_keys[good_keys %in% datasets_of_interest]
        pssm_datasets_out <- good_keys[good_keys %!in% datasets_of_interest]
    } else {
        pssm_datasets_in <- good_keys
        pssm_datasets_out <- c()
    }
    if (length(pssm_datasets_out) > 0) {
        cli_alert_info("Ignoring available pssm datasets {.val {pssm_datasets_out}}")
    }
    key_df_list <- lapply(pssm_datasets_in, function(k) read.delim(file.path(pssm_path, paste0(k, ".key")), header = F, col.names = c("key", "track", "bid")))
    data_df_list <- lapply(pssm_datasets_in, function(k) read.delim(file.path(pssm_path, paste0(k, ".data")), header = F, col.names = c("key", "pos", "A", "C", "G", "T")))
    is_non_negative <- sapply(data_df_list, function(ddf) all(apply(ddf[, c("A", "C", "G", "T")], 2, function(x) all(x >= 0))))
    if (length(which(!is_non_negative)) > 0) {
        cli_alert_warning("Dataset(s) {.val {pssm_datasets_in[!is_non_negative]}} had negative probabilities in pssms and were discarded from the analysis")
    }
    data_df_list <- data_df_list[is_non_negative]
    key_df_list <- key_df_list[is_non_negative]
    names(key_df_list) <- pssm_datasets_in[is_non_negative]
    names(data_df_list) <- names(key_df_list)
    return(list("keys" = key_df_list, "datasets" = data_df_list))
}

#' Utility function that makes fake intervals for chromosome/contig not present in an interval set
#' @description Explanation: as a validation step, gextract requires that every chromosome/contig in the ALLGENOME variable also be present in the track directory, so we add fake sequences in all missing chromosomes/contigs
#' @param peaks (optional) - a PeakIntervals or regular misha intervals object
#' @return the same object with fake peaks for missing chromosomes
#' @examples
#' \dontrun{
#' peaks_fix <- fix_missing_chroms_in_peaks(atac_sc@peaks)
#' }
#' @noRd
fix_missing_chroms_in_peaks <- function(peaks) {
    chroms_missing <- ALLGENOME[[1]]$chrom[ALLGENOME[[1]]$chrom %!in% unique(peaks$chrom)]
    if (length(chroms_missing) > 0) {
        fake_seqs <- gintervals.all() %>%
            filter(chrom %!in% peaks$chrom) %>%
            misha.ext::gintervals.centers() %>%
            mutate(start = start - 33, end = end + 34) %>%
            PeakIntervals() %>%
            mutate(peak_name = peak_names(.))
        peaks <- bind_rows(peaks, fake_seqs) %>% arrange(chrom, start)
        peaks$peak_name <- make.unique(peaks$peak_name)
        peaks <- PeakIntervals(peaks)
    }
    return(peaks)
}
