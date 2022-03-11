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
#' atac_sc <- import_from_10x("pbmc_data")
#' export_to_h5ad(atac_sc, "pbmc_data/atac_sc.h5ad")
#' }
#'
#' @inheritDotParams anndata::write_h5ad
#' @export
export_to_h5ad <- function(object, out_file, ...) {
    validate_atac_object(object)

    mat <- object$mat
    mat <- t(mat)

    peaks <- as.data.frame(object$peaks)
    rownames(peaks) <- peak_names(peaks)

    if (!is.null(object$metadata)) {
        metadata <- data.frame(rowname = colnames(object$mat)) %>%
            left_join(metadata) %>%
            column_to_rownames("rowname")
    } else {
        metadata <- NULL
    }

    cli_ul("Creating an AnnData object")
    adata <- anndata::AnnData(
        X = mat,
        var = peaks,
        obs = metadata,
        uns = list(class = class(object)[1])
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

export_atac_clust_ucsc <- function(mc_atac, track_prefix, output_dir = getwd(), clust_vec = NULL, normalization = 'none') {
    res_lst = prepare_clusters(mc_atac, clust_vec, normalization)
    purrr::walk(res_lst$clusts, function(cl) {
        atac_vec = res_lst$atac_mc_mat_clust[,cl]
        misha.ext::fwrite_ucsc(intervals = dplyr::mutate(mc_atac$peaks, 'score' = atac_vec), 
                            file = paste0(output_dir, '/', track_prefix, '_', gsub('\\/', '_', cl), '.ucsc'),
                            name = paste0(track_prefix, '_', cl),
                            color = res_lst$col_key[[as.character(cl)]],
                            type = 'bedGraph',
                            description = glue::glue('UCSC track for cluster {cl} of dataset {track_prefix}')
        )
    })
}



#' Generate a misha track for each atac metacell cluster
#'
#' @description generate a track for each metacell cluster, of the form \code{track_prefix.name}, where names
#' are given at \code{clust_names}
#'
#' @param mc_atac_mat metacell ATAC matrix (like mc_atac$mat)
#' @param clust_vec (optional) a vector of length #metacells representing an annotation/clustering (can be output of \code{gen_atac_mc_clust})
#' @param track_prefix prefix for generated misha tracks.
#' @param parallel (optional) run function with parallel processing
#' @param num_cores (required if parallel == TRUE) number of cores to use for parallel processing

#' @return track_names names of generated misha tracks
#' @export
export_atac_clust_misha <- function(mc_atac, track_prefix, clust_vec = NULL, normalization = 'none', parallel = FALSE, num_cores = NULL) {
    res_lst = prepare_clusters(mc_atac, clust_vec, normalization)
    if (parallel) {
        new_tracks = parallel::mclapply(res_lst$clusts, function(cl) write_cluster_misha_track(cl, res_lst$atac_mc_mat_clust, track_prefix), mc.cores = num_cores)
    }
    else  {
        new_tracks = sapply(res_lst$clusts, function(cl) write_cluster_misha_track(cl, res_lst$atac_mc_mat_clust, track_prefix))
    }
    gdb.reload()
}

#' Generate a misha track for a metacell cluster (backend function)
#'
#' @description generates a track for a metacell cluster
#'
#' @param cl ATAC metacell cluster
#' @param atac_mc_mat_clust a matrix with ATAC signal averaged across metacells in each cluster
#' @param track_prefix prefix for track name

#' @return trknm - the name of the generated track in the misha database
write_cluster_misha_track <- function(cl, atac_mc_mat_clust, track_prefix) {
    atac_vec = atac_mc_mat_clust[,cl]
    cl = gsub('[\\/\\.-]', '_', cl)
    if (!is.null(track_prefix)) {trknm = paste0(track_prefix, '_', cl)}
    else {trknm = cl}
    if (!gtrack.exists(trknm)) {
        gtrack.create_sparse(trknm, 'ATAC signal',intervals = mc_atac$peaks, values = atac_vec)
    }
    return(trknm)
}

#' Prepare peak clusters for export (backend function)
#'
#' @param mc_atac McATAC object
#' @param clust_vec (optional) a vector of length #metacells representing an annotation/clustering (can be output of \code{gen_atac_mc_clust})
#' @param normalization (optional) normalization method, either 'none', 'lfcom' (log2 fold-change over median), 'zs' (Z-scores)

#' @return a list of:
#' @return atac_mc_mat_clust - ATAC signal per peak averaged over clusters
#' @return clusts - names of clusters
#' @return col_key - mapping of cluster names to colors
#'

prepare_clusters <- function(mc_atac, clust_vec = NULL, normalization = 'none') {
    print(normalization)
    if (is.null(clust_vec)) {
        if (all(!grepl('cell_type', colnames(mc_atac$metadata))) && all(!grepl('cluster_k_', colnames(mc_atac$peaks)))) {
            cli_abort('There is no "cell_type" or "cluster" field in metadata and no clustering vector was supplied')
        }
        else if (any(grep('^cell_type$', colnames(mc_atac$metadata)))) {
            clust_vec = unlist(mc_atac$metadata$cell_type)
        }
        else if (any(grep('^cluster_k_', colnames(mc_atac$metadata)))) {
            clust_vec = unlist(mc_atac$metadata[,grep('^cluster_k_', colnames(mc_atac$metadata))[[1]]])
        }
        else {
            warning(glue::glue('No clustering vector identified. Clustering with k == {round(ncol(atac_mc)/10)'))
            clust_vec = gen_atac_mc_clust(atac_mc, k = round(ncol(atac_mc)/10))
        }
    }
    if (all(!grepl('color', colnames(mc_atac$metadata)))) {
        num_clrs = length(unique(clust_vec))
        col_key = setNames(sample(grep('gray|white|grey', colors(), v=T), num_clrs), sort(unique(clust_vec)))
    }
    else {
        col_key = unique(mc_atac$metadata[,c('cell_type', 'color')])
        col_key = setNames(col_key$color, col_key$cell_type)
    }
    # if (!is.null(clust_names)) {names(clusts) = clust_names}
    eps = quantile(apply(mc_atac$mat, 1, mean), 0.05)
    if (normalization == 'lfcom') {
        atac_mc_mat = t(apply(mc_atac$mat, 1, function(x) log2((x + eps)/median(x + eps))))
    }
    else if (normalization == 'zs') {
        atac_mc_mat = t(apply(mc_atac$mat, 1, function(x) (x - mean(x))/sd(x)))
    }
    else {atac_mc_mat = mc_atac$mat}
    atac_mc_mat_clust = t(tgs_matrix_tapply(atac_mc_mat, clust_vec, mean))
    # clusts = gsub('\\/', '_', sort(unique(clust_vec)))
    clusts = sort(unique(clust_vec))
    return(list('atac_mc_mat_clust' = atac_mc_mat_clust, 'clusts' = clusts, 'col_key' = col_key))
}