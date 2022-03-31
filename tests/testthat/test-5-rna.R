test_that("add_mc_rna works", {
    atac_mc <- add_mc_rna(atac_mc, rna_mc_mat)
    egc <- t(t(rna_mc_mat) / colSums(rna_mc_mat))
    both_mcs <- intersect(colnames(atac_mc@mat), colnames(egc))
    expect_equal(atac_mc@rna_egc, egc[, both_mcs])
})


test_that("add_mc_rna from metacell1 works", {
    metacell::scdb_init(fs::path(raw_dir, "scdb"), force_reinit = TRUE)
    mc_rna <- metacell::scdb_mc("rna")
    atac_mc <- add_mc_rna(atac_mc, mc_rna)
    egc <- mc_rna@e_gc
    both_mcs <- intersect(colnames(atac_mc@mat), colnames(egc))
    expect_equal(atac_mc@rna_egc[, both_mcs], egc[, both_mcs])
})

test_that("plot_atac_rna works", {
    atac_mc <- add_mc_rna(atac_mc, rna_mc_mat)
    plot_atac_rna(atac_mc, "CD4", "CD8A", plot_object_id = FALSE) %>% expect_ggplot_ok()
    plot_atac_rna(atac_mc, "CD4", "CD8A") %>% expect_ggplot_ok()
    plot_atac_rna(atac_mc, "CD4") %>% expect_ggplot_ok()
    plot_atac_rna(atac_mc, "CD4", peak = atac_mc@peaks$peak_name[1], plot_object_id = TRUE, normalize_atac = FALSE) %>% expect_ggplot_ok()
})

test_that("export of rna egc works", {
    atac_mc <- add_mc_rna(atac_mc, rna_mc_mat)
    export_to_h5ad(atac_mc, fs::path(raw_dir, "atac_mc.h5ad"), compression = "gzip")
    atac_mc1 <- import_from_h5ad(fs::path(raw_dir, "atac_mc.h5ad"))
    expect_true(mean(abs(atac_mc@rna_egc - atac_mc1@rna_egc)) <= 1e-9)
    expect_equal(rownames(atac_mc@rna_egc), rownames(atac_mc1@rna_egc))
    expect_equal(colnames(atac_mc@rna_egc), colnames(atac_mc1@rna_egc))
})
