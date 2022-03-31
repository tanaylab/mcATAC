#' Write a McATAC or ScATAC object to an h5ad file
#'
#'
#' @param object McATAC or ScATAC object
#' @param out_file name of the output file
#'
#' @return None.
#'
#' @examples
#' \dontrun{
#' atac_sc <- import_from_10x("pbmc_data", genome = "hg38")
#' export_to_h5ad(atac_sc, "pbmc_data/atac_sc.h5ad")
#' }
#'
#' @inheritDotParams anndata::write_h5ad
#' @export
export_to_h5ad <- function(object, out_file, ...) {
    validate_atac_object(object)

    mat <- object@mat
    mat <- t(mat)

    peaks <- as.data.frame(object@peaks)

    if (nrow(object@ignore_peaks) > 0) {
        peaks <- bind_rows(
            peaks %>% mutate(ignore = FALSE),
            object@ignore_peaks %>% mutate(ignore = TRUE)
        ) %>% as.data.frame()

        mat <- cbind(mat, t(object@ignore_pmat))
    }

    rownames(peaks) <- peak_names(peaks)

    if (!is.null(object@metadata)) {
        metadata <- data.frame(rowname = colnames(object@mat)) %>%
            bind_cols(object@metadata) %>%
            column_to_rownames("rowname")
    } else {
        metadata <- NULL
    }

    uns <- list(class = class(object)[1], genome = object@genome)

    # add all other slots
    slots <- methods::slotNames(object)
    slots <- slots[slots %!in% c(
        "egc", "fp", "mat", "peaks", "genome", "metadata",
        "ignore_peaks", "ignore_pmat", "rna_egc"
    )]
    for (s in slots) {
        uns[[s]] <- slot(object, s)
    }

    if (has_rna(object)) {
        uns[["rna_egc"]] <- object@rna_egc
        uns[["rna_mcs"]] <- colnames(object@rna_egc)
        uns[["rna_gene_names"]] <- rownames(object@rna_egc)
    }

    cli_ul("Creating an AnnData object")
    adata <- anndata::AnnData(
        X = mat,
        var = peaks,
        obs = metadata,
        uns = uns
    )

    cli_ul("Writing to file")
    anndata::write_h5ad(
        adata,
        out_file,
        ...
    )

    cli_alert_success("Successfully exported to {.file {out_file}}")
}

#' Generate a ucsc genome browser file for each metacell cluster
#'
#' @param mc_atac McATAC object
#' @param track_prefix prefix for generated misha tracks.
#' @param output_dir (optional) name of the directory to write files to
#' @param clust_vec (optional) a vector of length #metacells representing an annotation/clustering (can be output of \code{gen_atac_mc_clust})
#' @param normalization (optional) normalization method, either 'none', 'lfcom' (log2 fold-change over median), 'zs' (Z-scores)
#'
#' @export
export_atac_clust_ucsc <- function(mc_atac, track_prefix, output_dir = getwd(), clust_vec = NULL, normalization = "none") {
    res_lst <- prepare_clusters(mc_atac, clust_vec, normalization)
    fns <- purrr::map_chr(res_lst$clusts, function(cl) {
        atac_vec <- res_lst$atac_mc_mat_clust[, cl]
        fn <- file.path(output_dir, paste0(track_prefix, "_", gsub("\\/", "_", cl), ".ucsc"))
        misha.ext::fwrite_ucsc(
            intervals = mutate(mc_atac@peaks, "score" = atac_vec),
            file = fn,
            name = paste0(track_prefix, "_", cl),
            color = res_lst$col_key[[as.character(cl)]],
            type = "bedGraph",
            description = glue::glue("UCSC track for cluster {cl} of dataset {track_prefix}")
        )
        return(fn)
    })
    cli_alert_success("Successfully exported to ucsc. Files generated:")
    purrr::walk(fns, cli_ul)
}



#' Generate a misha track for each atac metacell cluster
#'
#' @description generate a track for each metacell cluster, of the form \code{track_prefix.name}, where names
#' are given at \code{clust_names}
#'
#' @param mc_atac_mat metacell ATAC matrix (like mc_atac@mat)
#' @param clust_vec (optional) a vector of length #metacells representing an annotation/clustering (can be output of \code{gen_atac_mc_clust})
#' @param track_prefix prefix for generated misha tracks.
#' @param description (optional) description for tracks (can be a glue-formatted expression)
#' @param parallel (optional) run function with parallel processing
#' @param num_cores (required if parallel == TRUE) number of cores to use for parallel processing

#' @return track_names names of generated misha tracks
#' @export
export_atac_clust_misha <- function(mc_atac, track_prefix, description = NULL, clust_vec = NULL, normalization = "none", parallel = FALSE, num_cores = NULL) {
    res_lst <- prepare_clusters(mc_atac, clust_vec, normalization)
    if (!is.null(description)) {
        if (!is.null(clust_vec) && length(description) == 1) {
            description <- rep(description, length(unique(clust_vec)))
        } else {
            assert_that(!is.null(clust_vec), length(description) == length(unique(clust_vec)))
        }
    } else {
        description <- paste0(track_prefix, res_lst$clusts)
    }
    if (parallel) {
        track_names <- parallel::mclapply(seq_along(res_lst$clusts), function(cl, i) {
            write_cluster_misha_track(cl[[i]], res_lst$atac_mc_mat_clust, track_prefix, description = description[[i]])
        }, cl = res_lst$clusts, mc.cores = num_cores)
    } else {
        track_names <- sapply(seq_along(res_lst$clusts), function(cl, i) {
            write_cluster_misha_track(cl[[i]], res_lst$atac_mc_mat_clust, track_prefix, description = description[[i]])
        }, cl = res_lst$clusts)
    }
    gdb.reload()
    return(track_names)
}

#' Generate a misha track for a metacell cluster (backend function)
#'
#' @description generates a track for a metacell cluster
#'
#' @param cl ATAC metacell cluster
#' @param atac_mc_mat_clust a matrix with ATAC signal averaged across metacells in each cluster
#' @param track_prefix prefix for track name
#' @param description (optional) description for track (can be a glue-formatted expression)
#' @param override (optional) whether to override an existing track; warning: may slow runtime substantially (due to reloading misha db)
#' @return trknm - the name of the generated track in the misha database
#' @export
write_cluster_misha_track <- function(cl, atac_mc_mat_clust, track_prefix, description = NULL, override = FALSE) {
    atac_vec <- atac_mc_mat_clust[, cl]
    cl <- gsub("[\\/\\.-]", "_", cl)
    if (!is.null(track_prefix)) {
        trknm <- paste0(track_prefix, "_", cl)
    }
    if (is.null(description)) {
        description <- trknm
    } else {
        trknm <- cl
    }
    if (!gtrack.exists(trknm)) {
        gtrack.create_sparse(track = trknm, description = glue::glue(description), intervals = mc_atac@peaks, values = atac_vec)
    } else {
        if (override) {
            gtrack.rm(trknm, force = TRUE)
            gdb.reload()
            gtrack.create_sparse(track = trknm, description = glue::glue(description), intervals = mc_atac@peaks, values = atac_vec)
        } else {
            cli_alert_warning("Track exists and {.code override = FALSE}, no track written")
        }
    }
    return(trknm)
}

#' Prepare peak clusters for export (backend function)
#'
#' @param mc_atac McATAC object
#' @param clust_vec (optional) a vector of length #metacells representing an annotation/clustering (can be output of \code{gen_atac_mc_clust})
#' @param normalization (optional) normalization method, either 'none', 'lfcom' (log2 fold-change over median), 'zs' (Z-scores)
#' @param eps_q if \code{normalization == 'lfcom'}, use quantile \code{eps_q} (default = 0.05) for regularizing log expression
#' @return a list of:
#' \itemize{
#' \item{atac_mc_mat_clust: }{ATAC signal per peak averaged over clusters}
#' \item{clusts: }{names of clusters}
#' \item{col_key: }{mapping of cluster names to colors}
#' }
#'
#' @export
prepare_clusters <- function(mc_atac, clust_vec = NULL, normalization = "none", eps_q = 0.05) {
    if (is.null(clust_vec)) {
        if (!any(has_name(mc_atac@metadata, c("cell_type", "cluster")))) {
            cli_abort('There is no "cell_type" or "cluster" field in metadata and no clustering vector was supplied')
        } else if (has_name(mc_atac@metadata, "cell_type")) {
            clust_vec <- unlist(mc_atac@metadata$cell_type)
        } else if (has_name(mc_atac@metadata, "cluster_")) {
            clust_vec <- unlist(mc_atac@metadata[, grep("cluster_", colnames(mc_atac@metadata))[[1]]])
        } else {
            print(round(ncol(mc_atac@mat) / 10))
            cli_warn(glue::glue("No clustering vector identified. Clustering with k == {round(ncol(mc_atac@mat)/10)}"))
            clust_vec <- gen_atac_mc_clust(mc_atac, k = round(ncol(mc_atac@mat) / 10))
        }
    }
    if (!has_name(mc_atac@metadata, "color")) {
        num_clrs <- length(unique(clust_vec))
        col_key <- setNames(chameleon::distinct_colors(num_clrs)$name, sort(unique(clust_vec)))
    } else {
        col_key <- tibble::deframe(unique(mc_atac@metadata[, c("cell_type", "color")]))
    }
    eps <- quantile(Matrix::rowMeans(mc_atac@mat), eps_q)
    if (normalization == "lfcom") {
        atac_mc_mat <- t(apply(mc_atac@mat, 1, function(x) log2((x + eps) / median(x + eps))))
        cli_alert_info("Using eps_q={eps_q} and eps = {eps} for regularization")
    } else if (normalization == "zs") {
        atac_mc_mat <- t(apply(mc_atac@mat, 1, function(x) (x - mean(x)) / sd(x)))
    } else {
        atac_mc_mat <- mc_atac@mat
    }
    atac_mc_mat_clust <- t(tgs_matrix_tapply(atac_mc_mat, clust_vec, mean))
    clusts <- sort(unique(clust_vec))
    return(list("atac_mc_mat_clust" = atac_mc_mat_clust, "clusts" = clusts, "col_key" = col_key))
}
