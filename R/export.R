#' Write a McATAC or ScATAC object to an h5ad file
#'
#'
#' @param object McATAC or ScATAC object
#' @param out_file name of the output file
#'
#' @return None.
#'
#' @export
export_to_h5ad <- function(object, out_file) {

}

#' Generate a ucsc genome browser file for each metacell cluster
#'
#' @param mc_atac McATAC object
#' @param track_prefix prefix for generated misha tracks.
#' @param output_dir (optional) name of the directory to write files to
#' @param clust_vec (optional) a vector of length #metacells representing an annotation/clustering 
#'   (can be output of \code{gen_atac_mc_clust})
# @param clust_names (optional) a vector of length(unique(clust_vec)) with names for each metacell cluster
#' @param normalization (optional) normalization method, either 'none', 'lfcom' (log2 fold-change over median), 'zs' (Z-scores)

#'
#' @export

export_atac_clust_ucsc <- function(mc_atac, track_prefix, output_dir = '.', clust_vec = NULL, normalization = 'none') {
    res_lst = prepare_clusters(mc_atac, clust_vec, normalization)
    blabla = parallel::mclapply(res_lst$clusts, function(cl) {
    # blabla = purrr::walk(res_lst$clusts, function(cl) {
        atac_vec = res_lst$atac_mc_mat_clust[,cl]
        misha.ext::fwrite_ucsc(intervals = dplyr::mutate(mc_atac$peaks, 'score' = atac_vec), 
                               file = paste0(output_dir, '/', track_prefix, '_', gsub('\\/', '_', cl), '.ucsc'),
                               name = paste0(track_prefix, '_', cl),
                               color = res_lst$col_key[[as.character(cl)]],
                               type = 'bedGraph',
                               description = glue::glue('UCSC track for cluster {cl} of dataset {track_prefix}')
        )
    # })
    }, mc.cores = 20)
}

#' Generate a misha track for each atac metacell cluster
#'
#' @description generate a track for each metacell cluster, of the form \code{track_prefix.name}, where names
#' are given at \code{clust_names}
#'
#' @param mc_atac_mat metacell ATAC matrix (like mc_atac$mat)
#' @param clust_vec (optional) a vector of length #metacells representing an annotation/clustering 
#'   (can be output of \code{gen_atac_mc_clust})
# @param clust_names (optional) a vector of length(unique(clust_vec)) with names for each metacell cluster
#' @param track_prefix prefix for generated misha tracks.
#'
#' @export
export_atac_clust_misha <- function(mc_atac, track_prefix, clust_vec = NULL, normalization = 'none') {
    gdb.reload()
    res_lst = prepare_clusters(mc_atac, clust_vec, normalization)
    # blabla = parallel::mclapply(res_lst$clusts, function(cl) {
    blabla = purrr::walk(res_lst$clusts, function(cl) {
        atac_vec = res_lst$atac_mc_mat_clust[,cl]
        cl = gsub('[\\/\\.-]', '_', cl)
        if (!is.null(track_prefix)) {trknm = paste0(track_prefix, '_', cl)}
        else {trknm = cl}
        if (!gtrack.exists(trknm)) {
            gtrack.create_sparse(trknm, 'ATAC signal',intervals = mc_atac$peaks, values = atac_vec)
        }
        return(trknm)
    # }, mc.cores = 20)
    })
}

#' Generate a ucsc genome browser file for each metacell cluster
#'
#' @param mc_atac McATAC object
#' @param track_prefix prefix for generated misha tracks.
#' @param output_dir (optional) name of the directory to write files to
#' @param clust_vec (optional) a vector of length #metacells representing an annotation/clustering 
#'   (can be output of \code{gen_atac_mc_clust})
# @param clust_names (optional) a vector of length(unique(clust_vec)) with names for each metacell cluster
#' @param normalization (optional) normalization method, either 'none', 'lfcom' (log2 fold-change over median), 'zs' (Z-scores)

#'
#' @export
prepare_clusters <- function(mc_atac, clust_vec = NULL, normalization = 'none') {
    print(normalization)
    if (is.null(clust_vec)) {
        if (all(!grepl('cell_type', colnames(mc_atac$metadata))) && all(!grepl('cluster_', colnames(mc_atac$peaks)))) {
            stop('There is no "cell_type" or "cluster" field in metadata and no clustering vector was supplied')
        }
        else if (any(grep('^cell_type$', colnames(mc_atac$metadata)))) {
            clust_vec = unlist(mc_atac$metadata$cell_type)
        }
        else if (any(grep('^cluster_', colnames(mc_atac$metadata)))) {
            clust_vec = unlist(mc_atac$metadata[,grep('^cluster_', colnames(mc_atac$metadata))[[1]]])
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