sc_counts <- scc_read(file.path(raw_dir, "reads"))
mc_counts <- scc_project_on_mc(sc_counts, cell_to_metacell_pbmc_example)

test_that("write_sparse_matrix_from_fragments works", {
    gset_genome("hg38")
    reg <- giterator.intervals(intervals = gintervals.all() %>%
        filter(chrom == "chr22"), iterator = 5e7) %>%
        slice(2)

    mm_file <- tempfile()
    fragments_file <- file.path(raw_dir, "fragments.tsv.gz")
    write_sparse_matrix_from_fragments(
        fragments_file,
        out_file = mm_file,
        cell_names = atac_sc,
        region = reg,
        use_tabix = TRUE,
        overwrite = TRUE
    )

    expect_true(file.exists(paste0(mm_file, ".gz")))
    mat <- tgutil::fread_mm(paste0(mm_file, ".gz"))
    colnames(mat) <- colnames(atac_sc@mat)
    expect_equal(dim(mat), c(818468, 11909))
    expect_equal(sum(mat), 501934)
    withr::with_options(list(scipen = 1e5), {
        df <- tgutil::fread(cmd = glue("tabix {fragments_file} {reg$chrom}:{reg$start}-{reg$end}"), header = FALSE, col.names = c("chrom", "start", "end", "id", "reads")) %>%
            as_tibble()
    })

    df_reg <- df %>% filter(id %in% colnames(atac_sc@mat))
    df_reg <- df_reg %>%
        tidyr::pivot_longer(cols = c("start", "end"), values_to = "start") %>%
        mutate(end = start + 1) %>%
        select(chrom, start, end, id) %>%
        misha.ext::gintervals.neighbors1(reg) %>%
        filter(dist == 0) %>%
        select(chrom, start, end, id)
    expect_equal(nrow(df_reg), sum(mat))

    df_bin <- extract_bin(mat, reg, colnames(atac_sc@mat), colnames(atac_sc@mat))
    expect_equal(
        df_reg %>%
            rename(metacell = id) %>%
            mutate(start = start - 1, end = end - 1) %>%
            count(chrom, start, end, metacell, name = "value") %>%
            anti_join(df_bin, by = c("chrom", "start", "end", "metacell", "value")) %>%
            nrow(),
        0
    )
})

test_that("write_sparse_matrix_from_bam works", {
    gset_genome("hg38")
    reg <- giterator.intervals(intervals = gintervals.all() %>%
        filter(chrom == "chr22"), iterator = 5e7) %>%
        slice(2)
    bam_file <- file.path(raw_dir, "possorted_bam.bam")
    mm_file <- tempfile()
    write_sparse_matrix_from_bam(bam_file,
        out_file = mm_file, cell_names = atac_sc,
        region = reg
    )
    expect_true(file.exists(paste0(mm_file, ".gz")))
    mat <- tgutil::fread_mm(paste0(mm_file, ".gz"))
    expect_equal(dim(mat), c(818468, 11909))
    expect_equal(sum(mat), 571580)

    # Extract the reads directly from the BAM file and compare
    cell_names <- colnames(atac_sc@mat)
    region_str <- paste0(reg$chrom, ":", reg$start, "-", reg$end)
    cmd <- glue("{samtools_bin} view --keep-tag CB --exclude-flags 1024 {bam_file} {region_str} | grep CB | {awk_cmd}", awk_cmd = "awk '{print $4,substr($12, 6)}'")
    df <- tgutil::fread(cmd = cmd, col.names = c("start", "cell_name"), header = FALSE) %>% as_tibble()
    df_reg <- df %>%
        mutate(chrom = reg$chrom, end = start + 1) %>%
        select(chrom, start, end, cell_name) %>%
        mutate(start = start - 1, end = end - 1) %>%
        misha.ext::gintervals.neighbors1(reg) %>%
        filter(dist == 0, start >= reg$start) %>%
        select(chrom, start, end, cell_name)

    df_count <- df_reg %>%
        filter(cell_name %in% cell_names) %>%
        count(chrom, start, end, cell_name)

    colnames(mat) <- colnames(atac_sc@mat)
    df_bin <- extract_bin(mat, reg, colnames(atac_sc@mat), colnames(atac_sc@mat)) %>%
        rename(cell_name = metacell, n = value)

    expect_equal(
        df_count %>%
            anti_join(df_bin, by = c("chrom", "start", "end", "cell_name", "n")) %>%
            nrow(),
        0
    )

    expect_equal(
        df_bin %>%
            anti_join(df_count, by = c("chrom", "start", "end", "cell_name", "n")) %>%
            nrow(),
        0
    )
})

test_that("mcc_to_mcatac works", {
    atac_mc_new <- mcc_to_mcatac(mc_counts, atac_sc@peaks)
    expect_equal(colnames(atac_mc@mat), colnames(atac_mc_new@mat))
    expect_true(atac_mc@peaks$peak_name %in% rownames(atac_mc_new@mat))
    intervs <- intersect(atac_mc@peaks$peak_name, atac_mc_new@peaks$peak_name)

    expect_equal(sum(abs(atac_mc@mat[intervs, ] - atac_mc_new@mat[intervs, ])), 0)
})

# skipped for now since it is heavy
# test_that("mcc_to_tracks works", {
#     gset_genome("hg38")
#     mcc_to_tracks(mc_counts, "temp.pbmc_mc_norm", overwrite = TRUE, normalize = TRUE)
#     mcc_to_tracks(mc_counts, "temp.pbmc_mc_new", overwrite = TRUE, normalize = FALSE)
#     a <- gextract(c("temp.pbmc_mc_norm.mc1", "temp.pbmc_mc_new.mc1"),
#         iterator = "temp.pbmc_mc_new.mc1",
#         intervals = "temp.pbmc_mc_new.mc1", colnames = c("norm", "raw")
#     ) %>% as_tibble()
#     expect_equal(sum(a$norm), 0)
#     expect_equal(
#         a %>%
#             mutate(norm1 = raw / sum(raw)) %>%
#             filter(abs(norm - norm1) >= 1e-9) %>%
#             nrow(),
#         0
#     )
# })
