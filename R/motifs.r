#' Generate peak-motif matrix
#'
#' This function calculates aggregated log-PSSM energies (as derived by the misha package's \code{create_pssm_energy} module) for a set of peaks 
#' (default - all peaks in the dataset) and for a set of motifs (default - all available)

#' @param atac (optional) - an ScATAC/McATAC object
#' @param peak_set (optional) - a PeakIntervals object
#' @param peak_width (optional) - size of region around peak centers to extract motif energies for
#' @param pssm_path (optional) - path to directory containing misha-formatted pssm files (e.g. \code{motifs.key} and \code{motifs.data})
#' @param datasets_of_interest (optional) - names of pssm datasets ([name].key-[name].data file combinations) to calculate PSSM values for
#' @param motif_tracks (optional) - misha track names for which to extract motif PSSMs
#' @param motif_regex (optional) - a vector of regular expressions for which to match motif track names and extract motif PSSMs
#' @param parallel (optional) - whether to use parallelize computations
#' @return a matrix of peaks (rows) vs. aggregated motif energies (columns)
#' @examples
#' \dontrun{
#'      peak_motif_mat <- generate_motif_pssm_matrix(peak_set=head(scatac@peaks), 
#'                                motif_regex = c("Bcl", "Atf"), 
#'                                datasets_of_interest = c("homer", "jaspar", "jolma"),
#'                                parallel = F)
#' }
#' @export
generate_motif_pssm_matrix <- function(atac = NULL, 
                                        peak_set = NULL, 
                                        peak_width = 200, 
                                        pssm_path = NULL, 
                                        datasets_of_interest = NULL, 
                                        motif_tracks = NULL, 
                                        motif_regex = NULL, 
                                        parallel = FALSE) {
    set.seed(1337)
    suffix <- stringi::stri_rand_strings(1,5)
    peaks <- get_peaks_for_pssm(atac, peak_set)
    if (!is.null(motif_tracks)) {
        tracks_exist <- sapply(motif_tracks, gtrack.exists())
        bad_tracks <- motif_tracks[!tracks_exist]
        cli_alert_warning("Tracks {.val {paste(head(bad_tracks))}} and {.val {pmax(0, length(bad_tracks) - 6)}} more were not found in misha gdb. Check that you are querying appropriate motif tracks.")
        tracks <- setdiff(motif_tracks, bad_tracks)
    }
    else if (!is.null(motif_regex)) {
        all_pssms <- get_available_pssms(pssm_path, datasets_of_interest)
        pssm_filt_inds <- lapply(all_pssms[["keys"]], function(kdf) {
            sapply(motif_regex, function(mre) {
                kdf$key[grep(mre, kdf$track, ignore.case = TRUE)]
            }) %>% unlist
        })
        keys_filt <- lapply(seq_along(all_pssms[["keys"]]), function (i) {
            filter(all_pssms[["keys"]][[i]], key %in% pssm_filt_inds[[i]])
        })
        datasets_filt <- lapply(seq_along(all_pssms[["datasets"]]), function (i) {
            filter(all_pssms[["datasets"]][[i]], key %in% pssm_filt_inds[[i]])
        })
        names(keys_filt) <- names(all_pssms[["keys"]])
        names(datasets_filt) <- names(all_pssms[["datasets"]])
        pssms <- list("keys" = keys_filt, "datasets" = datasets_filt)
    }
    else {
        cli_alert_info("No motif tracks or regular expressions were specified. Extracting PSSMs for all available motifs.")
        pssms <- get_available_pssms(pssm_path, datasets_of_interest)
    }
    peak_mids <- round((peaks$end + peaks$start)/2)
    peaks <- mutate(peaks, start = round(peak_mids-peak_width/2), end = round(peak_mids+peak_width/2))
    peaks <- fix_missing_chroms_in_peaks(peaks)
    res <- purrr::map2(pssms[["keys"]], names(pssms[["keys"]]), function(kdf, n) {
                        purrr::map2(.x = kdf$key, .y = kdf$track, .f = function(.x, .y) {
                            if (!dir.exists(file.path(GWD, n))) {dir.create(file.path(GWD, n))}
                            track_name <- paste0(n, ".", gsub("\\||\\.", "_", .y), "_", suffix)
                            if (!gtrack.exists(track_name)) {
                                if (parallel) {
                                    cmnd <- glue::glue("gtrack.create_pwm_energy(track = track_name, description = 'Temporary track created for track_name', pssmset = n, pssmid = .x, prior = 0.01, iterator = peaks")
                                    return(cmnd)
                                }
                                else {
                                    gtrack.create_pwm_energy(track = track_name, 
                                                        description = glue::glue("Temporary track created for {track_name}"), 
                                                        pssmset = n, 
                                                        pssmid = .x, 
                                                        prior = 0.01, 
                                                        iterator = peaks)
                                }
                            }
                            return(track_name)
                        })
                    }) %>% unlist
    if (parallel) {
        error_log <- gcluster.run(res)
        print(error_log)
    }
    trks_motifs <- gtrack.ls(paste0("_", suffix))
    if (length(trks_motifs) >= 1e+3) {
        trk_inds <- unlist(do.call('c', lapply(1:4, function(i) 
                    ifelse(i != 4, 
                            rep(i, ceiling(length(trks_motifs)/4)), 
                            rep(i, floor(length(trks_motifs)/4)))
                                    )
                    ))
        df_list <- lapply(sort(unique(trk_inds)), function(u) gextract(res[trk_inds == u], 
                                                                        intervals=peaks, 
                                                                        iterator=peaks, 
                                                                        colnames = gsub(paste0("_", suffix), "", trks_motifs[trk_inds == u])))
        peak_motif_matrix <- plyr::join_all(aa_lst, by=c('chrom', 'start', 'end'), type='left')
    }
    else {
        peak_motif_matrix <- gextract(res, intervals=peaks, iterator=peaks, colnames = gsub(paste0("_", suffix), "", trks_motifs))
    }
    dummy <- sapply(trks_motifs, gtrack.rm, force=TRUE)
    peak_motif_matrix <- filter(peak_motif_matrix, end - start >= peak_width)
    return(peak_motif_matrix)
}


#' Utility function for validating parameter combinations for PSSM extraction
#'
#' @param atac (optional) - an ScATAC/McATAC object
#' @param peak_set (optional) - a PeakIntervals object
get_peaks_for_pssm <- function(atac, peak_set) {
    if (is.null(atac) && is.null(peak_set)) {
        cli_abort("Must specify either {.var atac} or {.var peak_set} to calculate motif PSSMs")
    }
    else if (is.null(atac) && !is.null(peak_set)) {
        peaks <- peak_set
    }
    else if (!is.null(atac) && is.null(peak_set)) {
        peaks <- atac@peaks
    }
    else {
        cli_alert_warning("Both {.var atac} and {.var peak_set} were specified. Using only {.var atac@peaks} to calculate motif PSSMs.")
        peaks <- atac@peaks
    }
    return(peaks)
}

#' Utility function for validating \code{pssm_path} and getting all PSSM tracks from path/gdb
#'
get_available_pssms <- function(pssm_path = NULL, datasets_of_interest = NULL) {
    if (is.null(pssm_path)) {
        pssm_path <- file.path(GROOT, "pssms")
        if (!dir.exists(pssm_path)) {
            cli_abort("The default {.var pssm_path} directory ({.val pssm_path}) does not exist. Make sure you have set up misha gdb correctly.")
        }
    }
    if (!dir.exists(pssm_path)) {
        cli_abort("The specified {.var pssm_path} directory ({.val pssm_path}) does not exist.")
    }
    pssm_file_list <- list.files(pssm_path)
    motif_keys <- grep(".key", pssm_file_list, value=TRUE)
    motif_datasets <- grep(".data", pssm_file_list, value=TRUE)
    unique_keys <- sort(unique(gsub(".key", "", motif_keys)))
    unique_datasets <- sort(unique(gsub(".data", "", motif_datasets)))
    good_keys <- unique_keys[!is.na(match(unique_keys, unique_datasets))]
    if (length(good_keys) == 0) {
        cli_abort("There are no matching *.key-*.data file name combinations")
    }
    if (!is.null(datasets_of_interest)) {
        pssm_datasets_in <- good_keys[good_keys %in% datasets_of_interest]
        pssm_datasets_out <- good_keys[good_keys %!in% datasets_of_interest]
    }
    else {
        pssm_datasets_in <- good_keys
    }
    if (length(pssm_datasets_out) > 0) {cli_alert_info("Ignoring available pssm datasets {.val {pssm_datasets_out}}")}
    key_df_list <- lapply(pssm_datasets_in, function(k) read.delim(file.path(pssm_path, paste0(k, ".key")), header = F, col.names = c('key', 'track', 'bid')))
    data_df_list <- lapply(pssm_datasets_in, function(k) read.delim(file.path(pssm_path, paste0(k, ".data")), header = F, col.names = c('key', 'pos', 'A', 'C', 'G', 'T')))
    is_non_negative <- sapply(data_df_list, function(ddf) all(apply(ddf[,c("A","C","G","T")], 2, function(x) all(x >= 0))))
    if (length(which(!is_non_negative)) > 0) {cli_alert_warning("Dataset(s) {.val {pssm_datasets_in[!is_non_negative]}} had negative probabilities in pssms and are discarded from the analysis")}
    data_df_list <- data_df_list[is_non_negative]
    key_df_list <- key_df_list[is_non_negative]
    names(key_df_list) <- pssm_datasets_in[is_non_negative]
    names(data_df_list) <- names(key_df_list)
    return(list("keys" = key_df_list, "datasets" = data_df_list))
}

#' Utility function that makes fake intervals for chromosome/contig not present in an interval set
#' Explanation: as a validation step, gextract requires that every chromosome/contig in the ALLGENOME
#'              variable also be present in the track directory, so we add fake sequences in all 
#'              missing chromosomes/contigs
fix_missing_chroms_in_peaks <- function(peaks) {
    chroms_missing = ALLGENOME[[1]]$chrom[ALLGENOME[[1]]$chrom %!in% unique(peaks$chrom)]
    fake_seqs = do.call('rbind', lapply(chroms_missing, function(chrom) {
        coord = round(ALLGENOME[[1]]$end[match(chrom, ALLGENOME[[1]]$chrom)]/2)
        return(as.data.frame(list('chrom' = chrom, 'start' = coord - 34, 'end' = coord + 34)))
    }))
    fake_seqs$peak_name <- misha.ext::convert_misha_intervals_to_10x_peak_names(fake_seqs)
    peaks = arrange(bind_rows(peaks, fake_seqs), chrom, start)
    return(peaks)
}