library(mcATAC)
library(tidyverse)
library(misha)
library(regioneR)
library(rtracklayer)
devtools::load_all(export_all = FALSE)
# devtools::load_all(export_all = TRUE)
library(tidyverse)
library(tgutil)
library(here)


load("ab_ct_mm_egc.Rda")
load("ct_oc_egc.Rda")

mct <- mct_create("mm10", tracks = gtrack.ls("csf_reik_small.c"), metacells = gsub("csf_reik_small\\.", "", gtrack.ls("csf_reik_small.c")), id = "Reik Gastrulation multiome (cell types)")
mct_oc <- mct_create("GCF_009806435.1_UM_NZW_1.0_renamed", tracks = gtrack.ls("GD8E.c"), metacells = gsub("GD8E\\.", "", gtrack.ls("GD8E.c")), id = "Rabbit Gastrulation multiome (cell types)")


ct_colors_mm <- fread(here("work/color_types_mm.csv"))
mct <- add_mc_metadata(mct, ct_colors_mm %>% mutate(metacell = gsub("#", "c", color)))
ct_colors_oc <- fread(here("work/color_types_oc.tsv"))
mct_oc <- add_mc_metadata(mct_oc, ct_colors_oc %>% mutate(metacell = gsub("#", "c", color)))
mm2oc <- rtracklayer::import.chain("/home/ofirr/amos_ofirr/flo/fresh_chains/mm10.um_nzw_1.all/mm10.um_nzw_1.allfilled.chain")
oc2mm <- rtracklayer::import.chain("/home/ofirr/amos_ofirr/flo/fresh_chains/um_nzw_1.mm10.all/um_nzw_1.mm10.allfilled.chain")
mm2oc@metadata = list(genome1="mm10", genome2="GCF_009806435.1_UM_NZW_1.0_renamed")
oc2mm@metadata = list(genome1="GCF_009806435.1_UM_NZW_1.0_renamed", genome2="mm10")
# oc_peaks <- fread("/home/ofirr/amos_ofirr/cambridge_rabbit/umnzw_peaks_is.csv") %>% as_tibble()
oc_peaks <- fread("/home/ofirr/amos_ofirr/cambridge_rabbit/umnzw_peaks_is_p.csv") %>% as_tibble()
mm_peaks <- fread("/home/ofirr/amos_ofirr/cambridge_rabbit/mm10_peaks_is.csv") %>% as_tibble()
# mm_peaks <- fread("/home/aviezerl/proj/motif_reg/enhflow/data/all_peaks.csv") %>% as_tibble()

colnames(ct_mm) = gsub("#", "c", colnames(ct_mm))
colnames(ct_oc_egc) = gsub("#", "c", colnames(ct_oc_egc))

rna_mc_not_in_atac <- colnames(ct_mm)[colnames(ct_mm) %!in% mct@metacells]
atac_mc_not_in_rna <- mct@metacells[mct@metacells %!in% colnames(ct_mm)]
both_mcs <- intersect(mct@metacells, colnames(ct_mm))
mct@rna_egc <- ct_mm[, both_mcs]
rna_mc_not_in_atac <- colnames(ct_oc_egc)[colnames(ct_oc_egc) %!in% mct_oc@metacells]
atac_mc_not_in_rna <- mct_oc@metacells[mct_oc@metacells %!in% colnames(ct_oc_egc)]
both_mcs <- intersect(mct_oc@metacells, colnames(ct_oc_egc))
mct_oc@rna_egc <- ct_oc_egc[, both_mcs]
hc_mm <- group_to_hclust(ct_colors_mm %>% select(cell_type = color, group, order) %>% mutate(cell_type = gsub("^#", "c", cell_type), order = -order))
hc_oc <- group_to_hclust(ct_colors_oc %>% select(cell_type = color, group, order) %>% mutate(cell_type = gsub("^#", "c", cell_type), order = -order))
get_peaks_data = function(mct, peaks){
    gset_genome(mct@genome)
    for(mcc in mct@tracks){
        gvtrack.create(vtrack = paste0(mcc, "_sum"), src = mcc, func = "sum")
    }
    mat_peaks <- gextract(
        paste0(mct@tracks, "_sum"),
        intervals=peaks,
        iterator=peaks,
        colnames = mct@metacells)
    mat_peaks# %>% misha.ext::intervs_to_mat(remove_intervalID = TRUE)
}

gset_genome(mct@genome)
mm_mat_peaks = get_peaks_data(mct, mm_peaks)

gset_genome(mct_oc@genome)
oc_mat_peaks = get_peaks_data(mct_oc, gintervals.force_range(as.data.frame(oc_peaks)))

intersecting_mcs = setdiff(intersect(colnames(mm_mat_peaks), colnames(oc_mat_peaks)), c("chrom", "start", "end", "intervalID"))
mm_mat_peaks[is.na(mm_mat_peaks)] <- 0
oc_mat_peaks[is.na(oc_mat_peaks)] <- 0

run_app(mct, mct2 = mct_oc, chain = mm2oc, chain2 = oc2mm, port = 6663, annotations1 = mm_peaks, annotations2 = oc_peaks, all_peaks = mm_mat_peaks, all_peaks2 = oc_mat_peaks, hc = hc_mm, hc2 = hc_oc)