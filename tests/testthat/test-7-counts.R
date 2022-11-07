sc_counts <- scc_read(file.path(raw_dir, "reads"))
mc_counts <- scc_to_mcc(sc_counts, cell_to_metacell_pbmc_example)

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
    expect_true(all(atac_mc@peaks$peak_name %in% rownames(atac_mc_new@mat)))
    intervs <- intersect(atac_mc@peaks$peak_name, atac_mc_new@peaks$peak_name)

    expect_equal(sum(abs(atac_mc@mat[intervs, ] - atac_mc_new@mat[intervs, ])), 0)
})

test_that("create_smoothed_track_from_dataframe works", {
    prev_groot <- GROOT
    withr::defer(gsetroot(GROOT))
    misha.ext::gset_genome("hg38")
    temp_track <- temp_track_name("temp.")
    df <- giterator.intervals(iterator = 1, intervals = gintervals(1, 0, 100)) %>%
        mutate(value = ifelse(start %in% 40:60, 1, 0)) %>%
        as_tibble()
    create_smoothed_track_from_dataframe(df, track_prefix = "temp", track = temp_track, description = "", window_size = 5, resolution = 1, overwrite = TRUE)
    a <- gextract(temp_track, gintervals(1, 0, 100), colnames = "v2")
    df1 <- df %>%
        mutate(v1 = zoo::rollsum(value, k = 11, na.pad = TRUE)) %>%
        left_join(a)
    expect_equal(df1 %>% na.omit() %>% filter(v1 != v2) %>% nrow(), 0)
})

# A heavy test - uncomment to run
# test_that("mcc_to_tracks works", {
#     mct <- mcc_to_tracks(mc_counts, "pbmc_mc", overwrite = TRUE, resolution = 20, window_size = NULL)
#     expect_equal(mc_counts@cell_names,mct@metacells)
#     expect_true(all(gtrack.exists(mct@tracks)))
#     expect_equal(mc_counts@genome, mct@genome)
#     expect_equal(gsummary(mct@tracks[6])[[5]], mct@total_cov[6])
# })
