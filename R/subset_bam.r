#' Subset 10X atac_possorted.bam file into per-metacell BAM files
#' 
#' @param bam_path path to the 10X atac_possorted.bam file
#' @param out_dir (output) directory to output per-metacell bam files
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
    error_log <- withr::with_path(new = "/home/feshap/src/parallel-20220322/src/parallel", 
                                    code = system2(command = fp_sb2mc,
                                                    args = c(bam_path, 
                                                             c2mc_path, 
                                                             out_dir,
                                                             fp_rgp))
                                )
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
